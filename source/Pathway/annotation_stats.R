#' @title Prepare Annotation Database
#' @description Loads KEGG database files and species-specific mapping files for pathway annotation.
#' @param omics Character. 'metab' for metabolomics, 'trans' for transcriptomics.
#' @param species Character. Species code (e.g., 'hsa', 'mmu').
#' @return A list containing:
#'   - tab_data: Data frame mapping Compound/Gene IDs to KEGG Orthology (KO).
#'   - gff_data: GFF data containing KEGG pathway information.
#'   - CtoMap: Named vector mapping KEGG IDs to Pathway strings.
annotationdatabase<-function(omics,species){
  # Define path to database resources
  path<-"./source/Database"
  
  # Configure file names and column indices based on omics type
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
  
  # Initialize result container
  result <- data.frame(CPD = character(), KO = character(), stringsAsFactors = FALSE)
  ko_dict <- list()
  
  # Parse mapping file: Map Compound/Gene to KO
  for (line in data) {
    cpd <- line[1]
    ko <- line[ko_col]
    if (ko %in% names(ko_dict)) {
      ko_dict[[ko]] <- c(ko_dict[[ko]], cpd)
    } else {
      ko_dict[[ko]] <- cpd
    }
  }
  
  # Flatten dictionary to data frame
  for (ko in names(ko_dict)) {
    result <- rbind(result, data.frame(CPD = ko, KO = paste(ko_dict[[ko]], collapse = "\t"), stringsAsFactors = FALSE))
  }

  # Clean up IDs
  result <- result[order(result$CPD), ]
  result <- unique(result)
  result$CPD <- gsub("cpd:", "",result$CPD)

  result$CPD <- gsub("ko:", "",result$CPD)
  result$CPD <- sub(" .*", "", result$CPD)
  result$KO <- gsub("\t", " ",result$KO)
  tab_data <- result

  # Load GFF Pathway Data
  gff_data_file<-file.path(path,paste0("KEGG/",gff_data_filename))
  gff_data <- suppressMessages(readr::read_delim(gff_data_file, delim = "\t")) |> 
    as.data.frame()

  CtoMap0 <- gff_data |> dplyr::mutate(Map = ifelse(grepl(gff_data_grep, PATHWAY), PATHWAY, NA))
  CtoMap <- setNames(CtoMap0$Map, CtoMap0$KEGG_ID)
  
  return( list(tab_data=tab_data, gff_data=gff_data, CtoMap=CtoMap))
 
}

#' @title Run Pathway Annotation Count
#' @description Maps input features to KEGG pathways and counts the number of hits per pathway.
#' @param selectCid_data Data frame containing input features (KEGG IDs).
#' @param database List object returned by annotationdatabase().
#' @param omics_type Character. 'metab' or 'trans'.
#' @param skip Logical. If TRUE (default), excludes global/broad metabolic pathways
#' that are too general for meaningful interpretation. These include:
#' - map01100: Metabolic pathways
#' - map01110: Biosynthesis of secondary metabolites
#' - map01120: Microbial metabolism in diverse environments
#' Setting to FALSE will include these pathways in the output.
#' @return A data frame containing:
#'   - Pathway: Pathway Name
#'   - Count: Number of mapped input features
#'   - PathwayID: KEGG Pathway ID 
#'   - KEGG IDs: Semicolon-separated list of mapped feature IDs
run_annotation_count <- function(selectCid_data=NULL, database=NULL,omics_type=NULL,skip=TRUE) {
  selectCid_data=selectCid_data[complete.cases(selectCid_data), ]
    tab_data=database$tab_data
    gff_data=database$gff_data
    CtoMap=database$CtoMap

    # Define patterns based on omics type
    if(all(grepl("^metab",omics_type))){
      CK<-"^C"
      map_ko<-"map"
    }else if(all(grepl("^trans",omics_type))){
      CK<-"^K"
      map_ko<-"ko"
    }
    
    # Initialize data structures for mapping
    c_ko <- list()
    path <- list()
    sum1 <- nrow(selectCid_data)
    sum2 <- nrow(gff_data)
    
    # Process tab file: Build dictionary of Compound->KO
    for (i in 1:nrow(tab_data)) {
      row <- tab_data[i, ]
      c_ko[[as.character(row$CPD)]] <- row$KO
    }
    
    # Process selectCid file: Map input features to Pathways
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
    
    # Process all GFF data: Count total features in each pathway (Background)
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
    
    # Initialize content variable
    content <- data.frame()
    
    # Sort pathways by Annotation Count (num1) descending
    # First, calculate counts for all paths to enable sorting
    for (pwID in names(path)) {
      path[[pwID]]$count <- length(path[[pwID]]$cId)
    }
    
    # Sort path list by count descending
    sorted_path <- path[order(sapply(path, function(x) x$count), decreasing = TRUE)]
    # Global maps to exclude (Too general)
    # 01100: Metabolic pathways
    # 01110: Biosynthesis of secondary metabolites
    # 01120: Microbial metabolism in diverse environments
    
    global_maps <- c("01100", "01110", "01120")
    # Generate output content
    for (pwID in names(sorted_path)) {
      if(skip){
        if (pwID %in% global_maps) next
      }
      num1 <- length(path[[pwID]]$cId) # Annotated Count
      num2 <- length(path[[pwID]]$cIdAll) # Background Count
      mapName <- path[[pwID]]$mapName
      cIds <- paste(path[[pwID]]$cId, collapse = "+")
      
      # No filters applied to match raw output logic
      temp <- data.frame(Pathway = mapName, Count = num1, 
                         PathwayID = paste0("map", pwID), `KEGG IDs` = cIds, check.names = FALSE)
      content <- dplyr::bind_rows(content, temp)
    }
    
    if (nrow(content) > 0) {
      plotdata <- content 
      # No longer arranging by Pvalue, already sorted by Count descending
      plotdata$Count <- as.numeric(plotdata$Count)
    } else {
      plotdata <- content
    }
    
    return(plotdata)

}