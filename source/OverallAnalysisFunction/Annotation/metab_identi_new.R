source("./source/OverallAnalysisFunction/Annotation/MSIdenti_DESI.R")
Dirdatabase <- "./source/Database"
run_sp<-function(data){
  sp_df <- as.data.frame(t(as.matrix(data@assays$Spatial$counts)))
  col<-colnames(sp_df)
  sp_df$x<-data$x
  sp_df$y<-data$y
  sp <- sp_df %>%
    dplyr::select(x,y,all_of(col))
  return(sp)
}
detect_input_type <- function(sp) {
  mz_cols <- colnames(sp)[!colnames(sp) %in% c("x", "y")]
  suppressWarnings(num_test <- as.numeric(mz_cols))
  valid_num_ratio <- sum(!is.na(num_test)) / length(mz_cols)
  if(valid_num_ratio > 0.8) {
    mz_values <- num_test[!is.na(num_test)]
    mz_range_ok <- mean(mz_values >= 50 & mz_values <= 2000) > 0.8
    
    if(mz_range_ok) {
      return("mz") 
    }
  }
  name_features <- grepl("[[:alpha:]]", mz_cols)  
  if(mean(name_features) > 0.8) {
    return("metab_name")
  }
  return("mz")
}
run_metab_name_processing <- function(sp,mode) {
  metab_info <- data.frame(
    Name = colnames(sp)[-c(1, 2)],  
    check.names = FALSE
  )
  iso_add_identi <- dataidenti2(databasepath=Dirdatabase,mode=mode,data0=metab_info) 
  iso_add_identi=iso_add_identi|>
    dplyr::select(-Name.y) |>
    dplyr::rename("metabolite"="Name.x")
  spotNb <- nrow(sp)
  stats_data <- apply(sp[, -c(1, 2)], 2, function(x) {
    non_miss_nb <- sum(!is.na(x))
    miss_ratio <- 1 - non_miss_nb / spotNb
    median_intensity <- median(x, na.rm = TRUE)
    return(c(miss_ratio = miss_ratio, median_intensity = median_intensity))
  })
  
  stats_data <- as.data.frame(t(stats_data))
  stats_data$metabolite <- rownames(stats_data)

  iso_add_identi <- merge(iso_add_identi, stats_data, by = "metabolite")

  final_result <- iso_add_identi %>%
    dplyr::select(
      metabolite, `Molecular Weight`, Formula, `KEGG ID`,
      `HMDB ID`, `Final Class`, Pathway, ID, 
      miss_ratio, median_intensity
    )
  
  return(final_result)
}

run_isotope<-function(sp,mode){

  sp6 <- sp
  spotNb <- nrow(sp6)
  mz_view_data <- apply(sp6[,-c(1,2)],2,function(x){
    non_miss_nb <- sum(!is.na(x))
    miss_ratio <- 1-non_miss_nb/spotNb
    total_intensity <- sum(x,na.rm = T)
    median_intensity <- median(x,na.rm = T)
    mz_wise <- c(non_miss_nb,miss_ratio,total_intensity,median_intensity)
    names(mz_wise) <- c("non_miss_nb","miss_ratio","total_intensity","median_intensity")
    return(mz_wise)
  })
  mz_view_data <- as.data.frame(t(mz_view_data))
  mz_view_data$mz <- as.numeric(row.names(mz_view_data))
  y <- mz_view_data[,c(5,4)]
  names(y) <- c("mz","intensity")
  isos <- isotopologues(y, ppm = 5)
  sp6[is.na(sp6)] <- 0.1
  for (i in seq_along(isos)) {
    z <- isos[[i]]
    m <- sp6[,c(z+2)] 
    co <- cor(m)
    co[co==1] <- 0
    high_co <- unique(as.vector(which(co > 0.5, arr.ind = TRUE)))
    if(length(high_co) == 0){
      isos[[i]] <- NA
    }else{
      isos[[i]] <- z[high_co]
    }
  }
  cleaned_list <- isos[!is.na(isos)]  
  sorted_list <- lapply(cleaned_list, sort)
  all_iso_mark <- data.frame()
  for(i in 1:length(sorted_list)){
    temp <- sp6[,sorted_list[[i]]+2] 
    temp_data <- data.frame(mz = names(temp),single_ios = min(as.numeric(names(temp))))
    all_iso_mark <- bind_rows(all_iso_mark,temp_data)
  }
  dim(all_iso_mark)

  mz_view_data$mz <- as.character(mz_view_data$mz)
  deal_data <- mz_view_data %>%
    dplyr::select(mz,miss_ratio,median_intensity) %>%
    left_join(all_iso_mark,c("mz"="mz")) %>%
    mutate(single_ios = ifelse(is.na(single_ios),"Unknown",single_ios)) 

  iso_filter_data <- deal_data %>%
    filter(single_ios!="Unknown") %>%
    arrange(single_ios) %>%
    group_by(single_ios) %>%
    filter(median_intensity == max(median_intensity)) %>%
    filter(mz!=single_ios) 

  deal_data_new <- deal_data %>%
    mutate(single_ios = ifelse(single_ios %in% iso_filter_data$single_ios,"Unknown",single_ios))

  Non_monoisotope <- deal_data_new %>%
    filter(single_ios != "Unknown") %>%
    filter(single_ios != mz)

  clean_data <- sp6 %>%
    dplyr::select(-Non_monoisotope$mz)
  combine<-list(deal_data_new,clean_data)
  return(combine)
}


run_add<-function(combine,mode){

  deal_data_iso <- combine[[1]] 
  sp5 <- combine[[2]] 
  mz <- as.numeric(names(sp5)[-c(1,2)])
  names(mz) <- mz
  if(mode == "pos"){
    add <- adductNames(polarity = "positive")[c(13,16,18,20,36)]
  }else{
    add <- adductNames(polarity = "negative")[c(3,5)]
  }

  combinations <- combn(add, 2, simplify = FALSE)
  add_all <- NULL
  for(x in c(1:length(combinations))){
    i <- combinations[[x]]
    mz_add <- as.data.frame(mz2mass(mz, i))
    mz_add$mz <- row.names(mz_add)
    mz_add <- reshape2::melt(mz_add,ic=c("mz"),variable.name = "adduct",value.name = "netural_mass")
    neutral_mass <- mz_add$netural_mass
    names(neutral_mass) <- paste(mz_add$mz,mz_add$adduct,sep = "_")

    diff_mz <- outer(neutral_mass, neutral_mass, function(a, b) abs((a - b) / b*1e6))
    diag(diff_mz) <- 10
    diff_mz_compare <- diff_mz<5

    true_indices1 <- which(diff_mz_compare, arr.ind = TRUE)
    if(nrow(true_indices1)>0){

      adduct_cor <- apply(true_indices1,1,function(x){
        a <- mz_add[x[1],1]
        b <- mz_add[x[2],1]
        k <- which(names(sp5) %in% c(a,b))
        if(length(k)==2){
          return(cor(sp5[,k[1]],sp5[,k[2]]))
        }else{
          return(1)
        }
      })

      true_indices2 <- as.data.frame(cbind(true_indices1,adduct_cor))
      true_indices2 <- true_indices2[true_indices2$adduct_cor>0.6,]

      if(nrow(true_indices2)>0){
        mz_add$mark <- "unpair"
        for(j in c(1:nrow(true_indices2))){
          a <- true_indices2[j,1]
          b <- true_indices2[j,2]
          mz_add$mark[c(a,b)] <- paste("Two",x,j,sep = "_")
        }
        mz_add2 <- mz_add[mz_add$mark != "unpair",]
        mz_add2 <- mz_add2[order(mz_add2$mark),]
        add_all <- rbind(add_all,mz_add2)
      }
    }
  }

  add_all$test <- paste(add_all$adduct,add_all$mark)
  d <- which(duplicated(add_all$test)==TRUE)
  add_all2 <- add_all[!(add_all$test %in% add_all$test[d]),]

  mz_add_nu <- aggregate(add_all2$adduct,by=list(add_all2$mz),function(x){
    length(unique(x))
  })
  unknow_mz <- mz_add_nu$Group.1[which(mz_add_nu$x>1)]
  unknown_mark <- add_all2$mark[add_all2$mz %in% unknow_mz]
  add_all2 <- add_all2[!(add_all2$mark %in% unknown_mark),]
  names(add_all2)[4] <- "adduct_mark"

  add_all2$mz_adduct <- paste(add_all2$mz,add_all2$adduct)
  add_all2$duplicated <- duplicated(add_all2$mz_adduct) | duplicated(add_all2$mz_adduct, fromLast = TRUE)
  d2 <- which(add_all2$duplicated ==TRUE)
  du_mz <- unique(add_all2$mz_adduct[add_all2$duplicated == TRUE])
  du_remove <- NULL
  for(j in c(1:length(du_mz))){
    i <- du_mz[j]
    row_nu <- which(add_all2$mz_adduct == i)
    marks <- add_all2$adduct_mark[row_nu]
    add_all2$adduct_mark[which(add_all2$adduct_mark %in% marks)] <- paste("Multiple_3_",j,sep = "")
    du_remove <- c(du_remove,row_nu[-1])
  }
  if(length(du_remove)>0){
    add_all3 <- add_all2[-du_remove,]
  }else{
    add_all3 <- add_all2
  }

  adduct_mark <- unique(add_all3$adduct_mark)
  for(i in adduct_mark){
    k <- which(add_all3$adduct_mark %in% i)
    sub_add <- add_all3[k,]
    netural <- mean(sub_add$netural_mass)
    add_all3$netural_mass[k] <- netural
  }
  add_all3$adduct <- as.character(add_all3$adduct)

  
  end_deal_data <- deal_data_iso %>%
    left_join(add_all3 %>% dplyr::select(mz,adduct,netural_mass,adduct_mark),c("mz"="mz")) %>%
    mutate(adduct = ifelse(is.na(adduct),"Unknown",adduct)) %>%
    mutate(netural_mass = ifelse(is.na(netural_mass),"Unknown",netural_mass)) %>%
    mutate(adduct_mark = ifelse(is.na(adduct_mark),"Unknown",adduct_mark))

  return(end_deal_data)
}

run_metab_identi<-function(end_deal_data,mode){

  end_deal_data0 <- end_deal_data 
  
  monoisotope_and_unknown <- end_deal_data0 %>% 
    filter(single_ios == "Unknown"|single_ios == mz)

  unadduct_data <- monoisotope_and_unknown %>%
    filter(adduct=="Unknown")

  adduct_data <- monoisotope_and_unknown %>%
    filter(adduct!="Unknown") %>%
    arrange(adduct_mark) %>%
    group_by(adduct_mark) %>%
    mutate(iosmark = ifelse(single_ios!="Unknown",2,1)) %>%
    filter(iosmark == max(iosmark)) %>%
    ungroup()

  adduct_and_iso_data <- adduct_data %>%
    filter(single_ios != "Unknown")

  adduct_and_not_iso_data <- adduct_data %>%
    filter(single_ios == "Unknown") %>%
    arrange(adduct_mark) %>%
    group_by(adduct_mark) %>%
    filter(median_intensity == max(median_intensity))

  adduct_to_identi_data <- bind_rows(adduct_and_iso_data,adduct_and_not_iso_data)
  iso_add_data0 <- tibble(mz=as.numeric(adduct_to_identi_data$mz),mw.all=adduct_to_identi_data$netural_mass)
  
  iso_add_identi <- dataidenti2(databasepath=Dirdatabase,mode=mode,data0=iso_add_data0)

  iso_add_identi$mz <- as.character(iso_add_identi$mz)
  if(length(grep("Formula",(names(iso_add_identi)))) >0){
    iso_add_identi$Formula <- gsub(" ","",iso_add_identi$Formula)
    iso_add_identi_deal <- iso_add_identi %>%
      left_join(adduct_to_identi_data,c("mz"="mz"))
  }else{
    iso_add_identi_deal <- iso_add_identi %>%
      left_join(adduct_to_identi_data,c("mz"="mz"))
  }
  if(mode == "pos"){
    iso_add_identi_deal$adduct <- gsub("\\[M\\+Na\\]\\+","Na1",iso_add_identi_deal$adduct)
    iso_add_identi_deal$adduct <- gsub("\\[M\\+K\\]\\+","K1",iso_add_identi_deal$adduct)
    iso_add_identi_deal$adduct <- gsub("\\[M\\+H\\]\\+","H1",iso_add_identi_deal$adduct)
    iso_add_identi_deal$adduct <- gsub("\\[M\\+NH4\\]\\+","N1H4",iso_add_identi_deal$adduct)
    iso_add_identi_deal$adduct <- gsub("\\[M\\+H-H2O\\]\\+","H1O1",iso_add_identi_deal$adduct)
  }else{
    iso_add_identi_deal$adduct <- gsub("\\[M-H\\]\\+","H1",iso_add_identi_deal$adduct)
    iso_add_identi_deal$adduct <- gsub("\\[M\\+Cl\\]\\+","Cl1",iso_add_identi_deal$adduct)
  }

  iso_and_not_add <- unadduct_data 

  if(mode == "neg"){
    not_add_data0 <- tibble(mz=as.numeric(iso_and_not_add$mz)) %>%
      mutate(H1=mz+1.0073,Cl1=mz-34.9694,H3O1=mz+19.0178)
  }else{
    not_add_data0 <- tibble(mz=as.numeric(iso_and_not_add$mz)) %>%
      mutate(H1=mz-1.0073,Na1=mz-22.9892,K1=mz-38.9632,N1H4=mz-18.0338,H1O1=mz+17.0033)
  }
  not_add_identi <- dataidenti2(databasepath=Dirdatabase,mode=mode,data0=not_add_data0) 

  not_add_identi$mz <- as.character(not_add_identi$mz)
  not_add_identi$Formula <- gsub(" ","",not_add_identi$Formula)
  not_add_identi_deal <- not_add_identi %>%
    left_join(iso_and_not_add,c("mz"="mz")) %>%
    mutate(adduct=Adduct)


  all_identi_mz <- bind_rows(iso_add_identi_deal,not_add_identi_deal) %>%
    mutate(addLevels=ifelse(adduct=="H1",1,
                            ifelse(adduct=="K1",2,
                                   ifelse(adduct=="Na1",3,
                                          ifelse(adduct=="H1O1",4,5))))) %>%
    dplyr::select(mz,`Molecular Weight`,Name,Formula,`KEGG ID`,`HMDB ID`,`Final Class`,Pathway,ID,massError,DatabaseLevel,miss_ratio,median_intensity,single_ios,adduct,netural_mass,adduct_mark,iosmark,addLevels)
  
  add_ones <- function(molecule){

    result <- gsub('([A-Z][a-z]?)(?![0-9])', '\\11', molecule, perl = TRUE)
    return(result)
  }

  all_identi_mz$Formula <- sapply(all_identi_mz$Formula, add_ones)
  all_identi_mz$Formula <- gsub("C1l","Cl",all_identi_mz$Formula)
  all_identi_mz$Formula <- gsub("N1a","Na",all_identi_mz$Formula)
  all_identi_mz$Formula <- gsub("S1i","Si",all_identi_mz$Formula)
  all_identi_mz$Formula <- gsub("B1r","Br",all_identi_mz$Formula)
  all_identi_mz$Formula_add <- NA

  for(i in 1:nrow(all_identi_mz)){
    chemforms <- all_identi_mz$Formula[i]
    adds <- all_identi_mz$adduct[i]
    if(adds %in% c("Na1","K1","N1H4","Cl1")){
      all_identi_mz$Formula_add[i] <- mergeform(chemforms,adds)
    }else if(adds %in% c("H3O1","H1O1")){
      if(check_ded(chemforms, adds) == "FALSE"){
        all_identi_mz$Formula_add[i] <- subform(chemforms,adds)
      }else{
        print(i)
        print("Erro:adduct")
      }
    }
    if(adds == "H1"){
      if(mode == "neg"){
        if(check_ded(chemforms, adds) == "FALSE"){
          all_identi_mz$Formula_add[i] <- subform(chemforms,adds)
        }else{
          print("Erro:adduct")
        }
      }else{
        all_identi_mz$Formula_add[i] <- mergeform(chemforms,adds)
      }
    }
  }
  all_identi_mz %<>% filter(!is.na(Formula_add))

  not_iso_data <- all_identi_mz %>%
    filter(single_ios == "Unknown")%>%
    arrange(mz) %>%
    group_by(mz) %>%
    filter(addLevels == min(addLevels)) %>%
    dplyr::slice(1)%>%
    ungroup()

  has_iso_data <- all_identi_mz %>%
    filter(single_ios != "Unknown") %>%
    arrange(mz)  %>%
    mutate(mz_ID=paste0(mz,"_",ID))
  data(isotopes)
  similarity_table <- data.frame()
  for(i in 1:nrow(has_iso_data)){

    print(i)
    chemforms <- has_iso_data$Formula_add[i]
    single_iso <- has_iso_data$mz[i]
    ID <- has_iso_data$ID[i]
    pattern<-isopattern(
      isotopes,
      chemforms,
      threshold=0.1,
      plotit=FALSE,
      charge=FALSE,
      emass=0.00054858,
      algo=1
    )

    pattern %<>% as.data.frame()
    names(pattern)[1] <- "mz"
    names(pattern)[2] <- "relative_intensity"
    pattern %<>% dplyr::select(mz,relative_intensity) %>%
      arrange(desc(relative_intensity)) 
    pattern$mz <- round(as.numeric(pattern$mz),4)
    if(nrow(pattern)>5){
      pattern_top5 <- pattern[1:5,]
    }else{
      pattern_top5 <- pattern
    }

    real_data <- end_deal_data0 %>%
      filter(single_ios == single_iso)
    real_data$minmz_intensity <- max(real_data$median_intensity)
    real_data$mz <- round(as.numeric(real_data$mz),4)

    real_data %<>% mutate(relative_intensity = (median_intensity / minmz_intensity) * 100) %>%
      dplyr::select(mz,relative_intensity) %>%
      arrange(desc(relative_intensity)) 

    peaks_a <- cbind(mz = pattern_top5$mz,intensity = pattern_top5$relative_intensity)
    peaks_b <- cbind(mz = real_data$mz,intensity <- real_data$relative_intensity)
    similarity_value <- round(msentropy_similarity(peaks_a, peaks_b, ms2_tolerance_in_da = 0.5)*100,2)
    temp <- data.frame(single_iso=single_iso,ID=ID,Isotope_similarity=similarity_value) %>%
      mutate(mz_ID=paste0(single_iso,"_",ID)) %>%
      dplyr::select(-ID)
    if(nrow(similarity_table)==0){
      similarity_table <- temp
    }else{
      similarity_table <- bind_rows(similarity_table,temp)
    }

    xmin <- round(min(c(pattern_top5$mz,real_data$mz)),0)-1.5
    xmax <- round(max(c(pattern_top5$mz,real_data$mz)),0)+1.5

  }

  has_iso_data %<>%
    left_join(similarity_table,c("mz_ID"="mz_ID")) %>%
    dplyr::select(-mz_ID,-single_iso) %>%
    arrange(mz) %>%
    group_by(mz) %>%
    filter(Isotope_similarity == max(Isotope_similarity)) %>%
    filter(addLevels == min(addLevels)) %>%
    dplyr::slice(1)%>%
    ungroup()

  more_to_one_identi <- bind_rows(has_iso_data,not_iso_data)
  dim(more_to_one_identi)
  return(more_to_one_identi)
}


run_more_to_one<-function(more_to_one_identi){

repeat_result <- more_to_one_identi

identi_result <- repeat_result %>%
  mutate(Isotope_similarity = ifelse(is.na(Isotope_similarity),0,Isotope_similarity)) %>%
  arrange(Name) %>%
  group_by(Name) %>%
  filter(DatabaseLevel == min(DatabaseLevel)) %>%
  filter(Isotope_similarity == max(Isotope_similarity)) %>%
  filter(addLevels == min(addLevels)) %>%
  filter(median_intensity == max(median_intensity)) %>%
  dplyr::slice(1)%>%
  ungroup()

KEGGnotNA <- identi_result %>%
  filter(is.na(`KEGG ID`))
KEGG <- identi_result %>%
  filter(!is.na(`KEGG ID`)) %>%
  arrange(`KEGG ID`) %>%
  group_by(`KEGG ID`) %>%
  filter(DatabaseLevel == min(DatabaseLevel)) %>%
  filter(Isotope_similarity == max(Isotope_similarity)) %>%
  filter(addLevels == min(addLevels)) %>%
  filter(median_intensity == max(median_intensity)) %>%
  dplyr::slice(1)%>%
  ungroup()
KEGGendidenti <- bind_rows(KEGG,KEGGnotNA)
KEGGendidenti$`Final Class` <- gsub("null",NA,KEGGendidenti$`Final Class`)

HMDBnotNA <- KEGGendidenti %>%
  filter(is.na(`HMDB ID`))
HMDB <- KEGGendidenti %>%
  filter(!is.na(`HMDB ID`)) %>%
  arrange(`HMDB ID`) %>%
  group_by(`HMDB ID`) %>%
  filter(DatabaseLevel == min(DatabaseLevel)) %>%
  filter(Isotope_similarity == max(Isotope_similarity)) %>%
  filter(addLevels == min(addLevels)) %>%
  filter(median_intensity == max(median_intensity)) %>%
  dplyr::slice(1)%>%
  ungroup()


HMDBendidenti <- bind_rows(HMDB,HMDBnotNA) %>% 
  dplyr::select(-DatabaseLevel,-miss_ratio,-median_intensity,-single_ios,-adduct_mark,-iosmark,-addLevels,-Formula_add,-adduct,-netural_mass) 
dim(HMDBendidenti)
return(HMDBendidenti)
}
