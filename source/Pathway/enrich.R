enrichdatabase<-function(omics,species){
  path<-"./source/Database"
  if(grepl("^metab",omics)){
    CPD_KO<-"cpd"
    ko_col<-2
    gff_data_filename<-"KEGG_106.txt"
    gff_data_grep<-"map"
  }else if(grepl("^trans",omics)){
    CPD_KO<-"ko"
    ko_col<-3
    gff_data_filename<-"ko_kegg_106.0_chart.txt"
    gff_data_grep<-"ko"
  }

  sp<-species
  lines_file <-paste0(path,"/komap/", sp,"/",sp, ".list")
  lines <- readLines(lines_file)
  filtered_lines <- lines[grep(CPD_KO, lines, ignore.case = TRUE)]
  data <- strsplit(filtered_lines, "\t")
  result <- data.frame(CPD = character(), KO = character(), stringsAsFactors = FALSE)
  ko_dict <- list()
  for (line in data) {
    cpd <- line[1]
    ko <- line[ko_col]
    if (ko %in% names(ko_dict)) {
      ko_dict[[ko]] <- c(ko_dict[[ko]], cpd)
    } else {
      ko_dict[[ko]] <- cpd
    }
  }
  for (ko in names(ko_dict)) {
    result <- rbind(result, data.frame(CPD = ko, KO = paste(ko_dict[[ko]], collapse = "\t"), stringsAsFactors = FALSE))
  }

  result <- result[order(result$CPD), ]
  result <- unique(result)
  result$CPD <- gsub("cpd:", "",result$CPD)

  result$CPD <- gsub("ko:", "",result$CPD)
  result$CPD <- sub(" .*", "", result$CPD)
  result$KO <- gsub("\t", " ",result$KO)
  tab_data <- result


  gff_data_file<-file.path(path,paste0("KEGG/",gff_data_filename))
  gff_data <- suppressMessages(readr::read_delim(gff_data_file, delim = "\t")) |> 
    as.data.frame()

  CtoMap0 <- gff_data |> dplyr::mutate(Map = ifelse(grepl(gff_data_grep, PATHWAY), PATHWAY, NA))
  CtoMap <- setNames(CtoMap0$Map, CtoMap0$KEGG_ID)
  return( list(tab_data=tab_data, gff_data=gff_data, CtoMap=CtoMap))
 
}

################################
run_enrich <- function(selectCid_data=NULL, database=NULL,omics_type=NULL) {
  selectCid_data=selectCid_data[complete.cases(selectCid_data), ]
    tab_data=database$tab_data
    gff_data=database$gff_data
    CtoMap=database$CtoMap

    if(all(grepl("^metab",omics_type))){
      CK<-"^C"
      map_ko<-"map"
    }else if(all(grepl("^trans",omics_type))){
      CK<-"^K"
      map_ko<-"ko"
    }
    
    # Initialize data structures
    c_ko <- list()
    path <- list()
    sum1 <- nrow(selectCid_data)
    sum2 <- nrow(gff_data)
    
    # Process tab file
    for (i in 1:nrow(tab_data)) {
      row <- tab_data[i, ]
      c_ko[[as.character(row$CPD)]] <- row$KO
    }
    
    # Process selectCid file
    for (i in 1:nrow(selectCid_data)) {
      row <- selectCid_data[i, ]
      if (grepl(CK, row$KEGG.ID)) {#"^C"
        if (length(grep(row$KEGG.ID, gff_data$KEGG_ID)) > 0) {
          if (grepl(map_ko, CtoMap[[row$KEGG.ID]])) {#"map"
            maps <- strsplit(CtoMap[[row$KEGG.ID]], ";")[[1]]
            for (map in maps) {
              map_match <- regmatches(map, regexec(paste0(map_ko,"(\\d+)\\s+(.*)"), map))[[1]]#"map(\\d+)\\s+(.*)"
              if (length(map_match) > 0) {
                mapId <- map_match[2]
                mapName <- map_match[3]
                if (exists(row$KEGG.ID, where = c_ko) && grepl(mapId, c_ko[[row$KEGG.ID]])) {
                  path[[mapId]]$mapName <- mapName
                  path[[mapId]]$cId <- unique(c(path[[mapId]]$cId, row[1]))
                }
              }
            }
          }
        }
      }
    }
    
    # Process all GFF data
    for (i in 1:nrow(gff_data)) {
      row <- gff_data[i, ]
      if (grepl(CK, row$KEGG_ID)) {#"^C"
        if (grepl(map_ko, row$PATHWAY)) {#"map"
          maps <- strsplit(row$PATHWAY, ";")[[1]]
          for (map in maps) {
            map_match <- regmatches(map, regexec(paste0(map_ko,"(\\d+)\\s+(.*)"), map))[[1]]#"map(\\d+)\\s+(.*)"
            if (length(map_match) > 0) {
              mapId <- map_match[2]
              if (exists(mapId, where = path)) {
                path[[mapId]]$cIdAll <- unique(c(path[[mapId]]$cIdAll, row$KEGG_ID))
              }
            }
          }
        }
      }
    }
    
    # Calculate P-values
    p_values <- numeric()
    for (pwID in sort(names(path))) {
      num1 <- length(path[[pwID]]$cId)
      num2 <- length(path[[pwID]]$cIdAll)
      p_value <- phyper(num1 - 1, num2, sum2 - num2, sum1, lower.tail = FALSE)
      p_values <- c(p_values, p_value)
    }
    
    # Initialize content variable
    content <- data.frame()
    
    # Add P-values to path list
    i <- 1
    for (pwID in sort(names(path))) {
      path[[pwID]]$pvalue <- p_values[i]
      i <- i + 1
    }
    
    # Sort by P-value
    sorted_path <- path[order(sapply(path, function(x) x$pvalue))]
    
    # Generate output content
    for (pwID in names(sorted_path)) {
      num1 <- length(path[[pwID]]$cId)
      num2 <- length(path[[pwID]]$cIdAll)
      mapName <- path[[pwID]]$mapName
      pvalue <- path[[pwID]]$pvalue
      cIds <- paste(path[[pwID]]$cId, collapse = "+")
      temp <- data.frame(Pathway = mapName, Count = num1, CountAll = num2, Pvalue = pvalue, PathwayID = paste0("map", pwID), `KEGG IDs` = cIds)
      content <- dplyr::bind_rows(content, temp)
    }
    
    if (any(!is.na(content))) {
      plotdata <- content |> dplyr::mutate(RichFactor = Count / CountAll) |>
        # dplyr::filter(Pvalue < enrichment_p_threshold) |>
        dplyr::arrange(Pvalue) 
      plotdata$Count <- as.numeric(plotdata$Count)
    } else {
      plotdata <- content
    }
    
    return(plotdata)

}