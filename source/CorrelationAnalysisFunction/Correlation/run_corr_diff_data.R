
run_cordata<-function(diff_m_cor,diff_t_cor){

  metab_diff1<- diff_m_cor[,-1]
  metab_diff2 <- as.data.frame(apply(metab_diff1, 2,as.numeric))
  rownames(metab_diff2)=as.vector(as.matrix(diff_m_cor[,1]))
  ko_diff1<- diff_t_cor[,-1]
  ko_diff2 <- as.data.frame(apply(ko_diff1, 2,as.numeric))
  rownames(ko_diff2)=as.vector(as.matrix(diff_t_cor[,1]))
  data<-list(metab_diff2,ko_diff2)

  return(data)

  
}

run_corr_diff_data <- function(omics_corr_data_save){
  metabdata=omics_corr_data_save[[1]] 
  transdata=omics_corr_data_save[[2]]

  metabdata$rowname <- rownames(metabdata)
  transdata$rowname <- rownames(transdata)

  merged_df <- merge(metabdata, transdata, by = "rowname", all = FALSE) 

  rownames(merged_df) <- merged_df$rowname
  merged_df$rowname <- NULL
  cor_table<-as.matrix(merged_df)
  t_its_cor <- Hmisc::rcorr(cor_table,type="spearman")
  CorrDF <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      from = rownames(cormat)[col(cormat)[ut]],
      to = rownames(cormat)[row(cormat)[ut]],
      cor =(cormat)[ut],
      p = pmat[ut]
    )
  }
  t_its_cor_df <- CorrDF(t_its_cor$r,t_its_cor$P)
  t_its_cor_df$padj <- p.adjust(t_its_cor_df$p, method = "bonferroni") #method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  
  return(t_its_cor_df)
}