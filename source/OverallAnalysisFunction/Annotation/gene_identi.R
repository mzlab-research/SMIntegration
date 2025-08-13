trans_identi_processing<-function(id,speciesname){
  root_path<-"./source/OverallAnalysisFunction/Annotation/"
  ko_genespath<-paste0(root_path,"/","ko_genes_definiton_new.list")#HSA,MMU

  KEGG<-read_delim(ko_genespath,delim="\t", col_names = FALSE)

  KEGG$X1 <- gsub("[ko:]", "", KEGG$X1)

  genedata<-data.frame(ENTREZID=paste0(speciesname,":",id$ENTREZID),gene=id$SYMBOL,id=id$ENTREZID)

  GENEMERGE<-genedata %>%
    left_join(KEGG,by=c("ENTREZID"="X2"))
 
  trans_identi<-GENEMERGE %>%
    dplyr::select(gene,id,X1) %>%
    distinct(gene, .keep_all = TRUE) 
  colnames(trans_identi)[3]<-"KEGG.ID"
return(trans_identi)
 
}
