#' @title Identify Differential Features
#' @description Performs differential expression/abundance analysis between two groups using Seurat's FindMarkers.
#' Uses the Wilcoxon rank-sum test by default.
#' @param data Seurat object containing the data.
#' @param group Vector of length 2. group[1] is the Treatment group, group[2] is the Control group.
#' @return A data frame of differential test results (p-values, logFC). Returns empty data frame on error.
find_marker<-function(data,group){
  obj <- data 
  # Set active identity to the provided grouping variable
  obj@meta.data$group<-as.factor(obj@meta.data$group)
  obj@active.ident<-obj$group

  # Try-catch to handle potential errors in FindMarkers (e.g., too few cells)
  markers <- tryCatch({
    # Run Wilcoxon test with minimal filtering (min.pct=0, logfc=0) to capture all results
    FindMarkers(obj,ident.1 = group[1],ident.2 = group[2],group.by = 'group', min.pct = 0, logfc.threshold =0)
  }, error = function(e) {
    return(data.frame())
  })

  return(markers)
}

#' @title Process and Classify Differential Markers
#' @description Post-processes FindMarkers results to classify features as Up-regulated, Down-regulated, or Non-significant
#' based on user-defined Fold Change and P-value thresholds.
#' @param markers Data frame returned by find_marker().
#' @param type Character. Feature type label (e.g., "gene" or "metabolite").
#' @param group Vector of length 2. Comparison groups.
#' @param FC_Threshold Numeric. Fold change threshold (linear scale).
#' @param pvalue Numeric. P-value threshold (usually adjusted).
#' @return A processed data frame with 'State' column (Up/Down/Non-significant) and formatted metrics.
find_marker_status<-function(markers,type,group,FC_Threshold,pvalue){
  # Check for empty markers input to prevent crash
  if(nrow(markers) == 0){
      # Return an empty dataframe with expected columns to prevent downstream errors
      empty_df <- data.frame(
          gene = character(),
          p_val = numeric(),
          avg_log2FC = numeric(),
          pct.1 = numeric(),
          pct.2 = numeric(),
          p_val_adj = numeric(),
          log_test = numeric(),
          `log2(Fold Change)` = numeric(),
          `Fold Change` = numeric(),
          State = factor(levels = c('Down', 'Non-significant', 'Up')),
          Sample = character(),
          check.names = FALSE
      )
      colnames(empty_df)[1] <- type
      return(empty_df)
  }

  Pvalue_type1="p_val_adj"
  # Calculate log2 thresholds
  fc<-abs(log2(as.numeric(FC_Threshold)))
  markers$gene<-rownames(markers)
  
  # Format numeric columns
  markers <- markers[order(abs(markers$avg_log2FC),decreasing=TRUE),]
  markers$avg_log2FC<-as.numeric(markers$avg_log2FC)
  markers$p_val_adj <-as.numeric(markers$p_val_adj )
  
  # Calculate -log10(p) and linear FC
  markers %<>%
    mutate(log_test = ifelse(.[[Pvalue_type1]]==0,0,(-log10(.[[Pvalue_type1]])))) %>%
    mutate(`log2(Fold Change)` = avg_log2FC) %>%
    mutate(`Fold Change` =  2^avg_log2FC)
  
  # Define significance state
  py<-abs(-log10(as.numeric(pvalue)))
  markers$State <- ifelse(markers$log_test >= py & abs(markers$`log2(Fold Change)`) >= fc , 
                          ifelse(markers$`log2(Fold Change)` > fc, "Up", "Down"), "Non-significant")
  markers$State <- factor(markers$State, levels = c('Down', 'Non-significant', 'Up'))
  markers$Sample<- paste0(group[1],":",group[2])
  
  # Select final columns
  markers<-markers %>%
    dplyr::select("gene","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","log_test","log2(Fold Change)","Fold Change","State","Sample") %>%
    filter(State %in% c("Up", "Down", "Non-significant"))
  
  # Rename ID column based on type
  markers$gene<-as.character(markers$gene)
  colnames(markers)[1]<-type
  return(markers)
}
