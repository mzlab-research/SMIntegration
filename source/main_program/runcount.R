runcount=function(data,sample){
  count=as.data.frame(t(data@assays$Spatial$counts))
  id<-colnames(count)
  count$x_y=rownames(count)
  samplename<-paste0("^",sample,":")
  count$x_y <- gsub(samplename, "", count$x_y)#1_1
  count %<>% separate(x_y, into = c("x", "y"),sep = "_",remove = FALSE)
  count$x_y <- str_replace_all(count$x_y, pattern = "(\\d+)_(\\d+)", replacement = "x\\1_y\\2")#x1_y1

  count %<>% 
    dplyr::select(x_y,x,y,all_of(as.character(id)))

  return(count)

}


omics_corr_processing1<-function(m_diff,t_diff,combine_count,samplelist,groupname){
  combine_count<-combine_count %>%
    dplyr::select(-x,-y)
    rownames(combine_count)<-combine_count$x_y
    count=combine_count
    samplelist_diff<-samplelist %>%
      filter(group %in% groupname)
    treat_samplelist_diff<-samplelist %>%
      filter(group %in% groupname[1])
    con_samplelist_diff<-samplelist %>%
      filter(group %in% groupname[2])
    #300
    if(nrow(m_diff) > 300){
      m_diff %<>% arrange(p_val_adj)
      m_diff <- m_diff[1:300,]
    }
    if(nrow(t_diff) > 300){
      t_diff %<>% arrange(p_val_adj)
      t_diff <- t_diff[1:300,]
    }
    #################runmean
    count_m_diff<- count %>% dplyr::filter(x_y %in% samplelist_diff$metabolomics) %>%
      dplyr::select(x_y,all_of(as.character(m_diff$metabolite)))
    rownames(count_m_diff)<- count_m_diff$x_y
    count_m_diff$x_y=NULL
    count_t_diff<- count %>% dplyr::filter(x_y %in% samplelist_diff$gene) %>%
      dplyr::select(x_y,all_of(as.character(t_diff$gene)))
    rownames(count_t_diff)<- count_t_diff$x_y
    count_t_diff$x_y=NULL
    metab_mean<-run_mean(core_table=count_m_diff,samplelist=samplelist,groupname=groupname)
    trans_mean<-run_mean(core_table=count_t_diff,samplelist=samplelist,groupname=groupname)
    #################treatment
    tr_count_m_diff<- count %>% dplyr::filter(x_y %in% treat_samplelist_diff$metabolomics) %>%
      dplyr::select(x_y,all_of(as.character(m_diff$metabolite)))

    rownames(tr_count_m_diff)<- tr_count_m_diff$x_y

    tr_count_t_diff<- count %>% dplyr::filter(x_y %in% treat_samplelist_diff$gene) %>%
      dplyr::select(x_y,all_of(as.character(t_diff$gene)))
    rownames(tr_count_t_diff)<- tr_count_t_diff$x_y

    ######################################control
    con_count_m_diff<- count %>% dplyr::filter(x_y %in% con_samplelist_diff$metabolomics) %>%
      dplyr::select(x_y,all_of(as.character(m_diff$metabolite)))
    rownames(con_count_m_diff)<- con_count_m_diff$x_y

    con_count_t_diff<- count %>% dplyr::filter(x_y %in% con_samplelist_diff$gene) %>%
      dplyr::select(x_y,all_of(as.character(t_diff$gene)))
    rownames(con_count_t_diff)<- con_count_t_diff$x_y

    tr_spatial_sampling <- spatial_sampling_correlation(data=tr_count_t_diff, sample_size = 1000)
    tr_count_m_diff <- tr_count_m_diff[tr_spatial_sampling, ]
    tr_count_t_diff <- tr_count_t_diff[tr_spatial_sampling, ]
    
    con_spatial_sampling <- spatial_sampling_correlation(data=con_count_t_diff, sample_size = 1000)
    con_count_m_diff <- con_count_m_diff[con_spatial_sampling, ]
    con_count_t_diff <- con_count_t_diff[con_spatial_sampling, ]
    #############
    data<-list(tr_count_m_diff,tr_count_t_diff,con_count_m_diff,con_count_t_diff,
               metab_mean,trans_mean)
  return(data)

}

spatial_sampling_correlation <- function(data, sample_size = 1000, method = "pearson") {

  if(sample_size > nrow(data)) {
    sample_size <- nrow(data)
  }

  sampled_indices <- sample(nrow(data), size = sample_size, replace = FALSE)

  return(sampled_indices)
}
run_mean<-function(core_table=NULL,samplelist=NULL,groupname=NULL){
  ###mean_table
  mean_dat<-as.data.frame(core_table)
  mean_table<- apply(mean_dat,2,function(x){scales::rescale(x, to = c(-1, 1))})
  mean_table<-mean_table |>
    t() |>
    as.data.frame(rownames.force = TRUE,colnames.force = TRUE)
  mean_table$feature_ID<-rownames(mean_table)
  mean_table<-mean_table |>
    dplyr::select(feature_ID, everything())     
  samplelist$sample<-samplelist$gene
  ##samplelist
  compare_table <- tibble::as_tibble_col(groupname) |>
    dplyr::group_by(.data$value) |> 
    dplyr::mutate(
      sample_list = list(
        samplelist[samplelist$group %in% unlist(stringr::str_split(.data$value, "\\+")), ]$sample
      )
    ) |>
    dplyr::mutate(sample_count = purrr::map_int(.data$sample_list, length)) 
  
  result_mean<-mean_table |>
    dplyr::mutate(
      case_norm_mean=rowMeans(pick(compare_table$sample_list[[1]])) ,
      control_norm_mean=rowMeans(pick(compare_table$sample_list[[2]]))
    )  |>
    dplyr::select(feature_ID,case_norm_mean,control_norm_mean)
  return(result_mean)
}