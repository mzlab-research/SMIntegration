dataidenti2 <- function(databasepath,mode,data0){

  database <- paste0(databasepath,"/animal_",mode,".txt")
  globel_data <- read_delim(database,delim="\t") %>%
    select(`Molecular Weight`,Name,Formula,`KEGG ID`,`HMDB ID`,`Final Class`,Pathway,ID) %>%
    as.data.frame()

  HMDB <- read_delim(file.path(databasepath,"HMDB/end-hmdb5.0.gff"),delim="\t")
  hmdb_class <- read_csv(file.path(databasepath,"HMDB/hmdb_class.csv"))
  hmdb_path <- read_csv(file.path(databasepath,"HMDB/hmdb_info_cas_pathway.csv"))
  HMDBdata <- HMDB %>%
    left_join(hmdb_class[,c(1,5)],c("hmdbid"="hmdbid")) %>%
    left_join(hmdb_path[,c(1,5)],c("hmdbid"="hmdbid")) %>%
    select(mono_mw,name,formula,kegg_id,hmdbid,family,PATHWAY) %>%
    dplyr::rename(`Molecular Weight`=mono_mw,Name=name,Formula=formula,`KEGG ID`=kegg_id,`HMDB ID`=hmdbid,`Final Class`=family,Pathway=PATHWAY) %>%
    mutate(ID=`HMDB ID`)

  KEGG <- read_delim(file.path(databasepath,"KEGG/KEGG_106.txt"),delim="\t")
  KEGG$`MONO_MW` <- gsub(";","",KEGG$`MONO_MW`)
  KEGG$`MONO_MW` <- gsub(" ","",KEGG$`MONO_MW`)
  KEGG$`MONO_MW` <- as.numeric(KEGG$`MONO_MW`)
  kegg_class <- read_csv(file.path(databasepath,"KEGG/kegg_class.csv"))
  KEGGdata <- KEGG %>%
    left_join(kegg_class[,c(1,5)],c("KEGG_ID"="KEGGID")) %>%
    select(MONO_MW,NAME,FORMULA,KEGG_ID,hmdbid,family,PATHWAY) %>%
    dplyr::rename(`Molecular Weight`=MONO_MW,Name=NAME,Formula=FORMULA,`KEGG ID`=KEGG_ID,`HMDB ID`=hmdbid,`Final Class`=family,Pathway=PATHWAY) %>%
    mutate(ID=`KEGG ID`)

  globel_data %<>%
    left_join(KEGG %>% select(KEGG_ID,MONO_MW) %>% as.data.frame(),c("ID"="KEGG_ID")) %>%
    left_join(HMDB %>% select(hmdbid,mono_mw),c("ID"="hmdbid")) %>%
    mutate(`Molecular Weight` = ifelse(!is.na(MONO_MW),MONO_MW,ifelse(is.na(MONO_MW)&!is.na(mono_mw),mono_mw,`Molecular Weight`))) %>%
    select(-MONO_MW,-mono_mw)
  if (!"mz" %in% colnames(data0)) {

    globel_data$Database <- "Project"
    KEGGdata$Database <- "KEGG"
    HMDBdata$Database <- "HMDB"
    combined_db <- bind_rows(globel_data, KEGGdata, HMDBdata)

    clean_name <- function(name) {
      name %>% 
        tolower() %>% 
        stringr::str_remove_all("[^[:alnum:]]") %>% 
        stringr::str_trim()
    }
    

    combined_db <- combined_db %>%
      mutate(CleanName = clean_name(Name))
    

    data0 <- data0 %>%
      mutate(CleanName = clean_name(Name))
    

    result <- data0 %>%
      left_join(combined_db, by = "CleanName") %>%
      arrange(Database) %>%  
      group_by(CleanName) %>%
      slice(1) %>%  
      ungroup() %>%
      mutate(
        DatabaseLevel = case_when(
          Database == "Project" ~ "1",
          Database == "KEGG" ~ "2",
          Database == "HMDB" ~ "3",
          TRUE ~ "Unknown"
        ),
        Adduct = "Direct Name Match",
        massError = 0  
      ) %>%
      dplyr::select(-CleanName, -Database)
    
    return(result)
  }

  MSidenti=function(id,Data,database){

    dat <- Data %>%
      filter(mz==id)

    DataIdenti <- data.frame()
    for(i in 2:ncol(dat)){
      mw <- dat[,i] %>% as.numeric()
      print(names(dat)[i])
      Adducts <- names(dat)[i]

      identi <- database %>%
        mutate(massError = abs(`Molecular Weight`-mw)) %>%
        filter(massError == min(massError)) %>%
        mutate(Adduct=Adducts) %>%
        filter(massError <0.003)
      if(nrow(identi)>0){  
        DataIdenti <- bind_rows(DataIdenti,identi)
      }
    }

    if(nrow(DataIdenti)>0){
      result <- cbind(dat,DataIdenti) 
    }else{
      result <- dat
    }
    return(result)
  }

  cpu = 4
  cl <- makeCluster(getOption("cl.cores", cpu))
  clusterExport(cl, c("MSidenti"),envir=environment())
  clusterEvalQ(cl,library("dplyr"))
  Res <- parLapply(cl,unique(data0$mz),MSidenti,Data=data0,database=globel_data)
  stopCluster(cl)
  endidenti0 <- as.data.frame(rbindlist(Res,fill=TRUE))
  if(length(grep("Name",names(endidenti0)))>0){
    endidenti0 %<>%
      filter(!is.na(Adduct)) %>%
      mutate(DatabaseLevel="1")
  }

  data1 <- data0 %>%
    filter(!mz %in% unique(endidenti0$mz))
  cpu = 4
  cl <- makeCluster(getOption("cl.cores", cpu))
  clusterExport(cl, c("MSidenti"),envir=environment())
  clusterEvalQ(cl,library("dplyr"))
  Res <- parLapply(cl,unique(data1$mz),MSidenti,Data=data1,database=KEGGdata)
  stopCluster(cl)
  endidenti1 <- as.data.frame(rbindlist(Res,fill=TRUE))
  if(length(grep("Name",names(endidenti1)))>0){
    endidenti1 %<>%
      filter(!is.na(Adduct)) %>%
      mutate(DatabaseLevel="2")
  }

  data2 <- data1 %>%
    filter(!mz %in% unique(endidenti1$mz))
  cpu = 4
  cl <- makeCluster(getOption("cl.cores", cpu))
  clusterExport(cl, c("MSidenti"),envir=environment())
  clusterEvalQ(cl,library("dplyr"))
  Res <- parLapply(cl,unique(data2$mz),MSidenti,Data=data2,database=HMDBdata)
  stopCluster(cl)
  endidenti2 <- as.data.frame(rbindlist(Res,fill=TRUE))
  if(length(grep("Name",names(endidenti2)))>0){
    endidenti2 %<>%
      filter(!is.na(Adduct)) %>%
      mutate(DatabaseLevel="3")
  }
  endidenti <- bind_rows(endidenti0,endidenti1,endidenti2)
  return(endidenti)
}
