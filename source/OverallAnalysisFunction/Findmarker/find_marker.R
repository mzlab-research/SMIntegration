find_marker<-function(data,group){
  obj <- data 
  obj@meta.data$group<-as.factor(obj@meta.data$group)
  obj@active.ident<-obj$group

  markers <- FindMarkers(obj,ident.1 = group[1],ident.2 = group[2],group.by = 'group', min.pct = 0, logfc.threshold =0)

  return(markers)
}
find_marker_status<-function(markers,type,group,FC_Threshold,pvalue){
  Pvalue_type1="p_val_adj"
  fc<-abs(log2(as.numeric(FC_Threshold)))
  markers$gene<-rownames(markers)
  
  markers <- markers[order(abs(markers$avg_log2FC),decreasing=TRUE),]
  markers$avg_log2FC<-as.numeric(markers$avg_log2FC)
  markers$p_val_adj <-as.numeric(markers$p_val_adj )
  markers %<>%
    mutate(log_test = ifelse(.[[Pvalue_type1]]==0,0,(-log10(.[[Pvalue_type1]])))) %>%
    mutate(`log2(Fold Change)` = avg_log2FC) %>%
    mutate(`Fold Change` =  2^avg_log2FC)
  py<-abs(-log10(as.numeric(pvalue)))
  markers$State <- ifelse(markers$log_test >= py & abs(markers$`log2(Fold Change)`) >= fc , 
                          ifelse(markers$`log2(Fold Change)` > fc, "Up", "Down"), "Non-significant")
  markers$State <- factor(markers$State, levels = c('Down', 'Non-significant', 'Up'))
  markers$Sample<- paste0(group[1],":",group[2])
  markers<-markers %>%
    dplyr::select("gene","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","log_test","log2(Fold Change)","Fold Change","State","Sample") %>%
    filter(State %in% c("Up", "Down", "Non-significant"))
  markers$gene<-as.character(markers$gene)
  colnames(markers)[1]<-type
  return(markers)
}
