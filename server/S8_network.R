########download

output$download_cornetwork_data1 <- downloadHandler(
  filename = function() {

    paste0("treatment_cornetwork_data.zip")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
      top_corrResultfilter=top_corrResultfilter()
      omics_corrResult<-omics_corrResult1()
      cornetwork_data<-list(omics_corrResult[[1]],top_corrResultfilter[[3]],top_corrResultfilter[[1]])
      tempdir <- setwd(tempdir())
      on.exit(setwd(tempdir))
      fi=c("CorResult.csv","network_nodes.csv","network_edges.csv")
      for (i in 1:length(fi)) {
        data<-cornetwork_data[[i]]
        write.csv(data, fi[i], row.names = FALSE)
      }
      zip(file,fi)
    })
  })
output$download_cornetwork_plot1 <- downloadHandler(
  filename = function() {
    "treatment_cornetwork_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      p <- diffnet_plotsave()[[1]]
      ggsave(file, plot = p, device = "png", width = 6, height = 6, dpi = 300,bg = "#FFFFFF")#
    })
  })


output$download_cornetwork_data2 <- downloadHandler(
  filename = function() {

    paste0("control_cornetwork_data.zip")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
      top_corrResultfilter=top_corrResultfilter()
      omics_corrResult<-omics_corrResult1()
      cornetwork_data<-list(omics_corrResult[[1]],top_corrResultfilter[[3]],top_corrResultfilter[[2]])
      tempdir <- setwd(tempdir())
      on.exit(setwd(tempdir))
      fi=c("CorResult.csv","network_nodes.csv","network_edges.csv")
      for (i in 1:length(fi)) {
        data<-cornetwork_data[[i]]
        write.csv(data, fi[i], row.names = FALSE)
      }
      zip(file,fi)
    })
  })
output$download_cornetwork_plot2 <- downloadHandler(
  filename = function() {
    "control_cornetwork_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      p <- diffnet_plotsave()[[2]]
      ggsave(file, plot = p, device = "png", width = 6, height = 6, dpi = 300,bg = "#FFFFFF")
    })
  })

###################test
# observe({
#   file="celldiff"
#   x<-readRDS(file.path(file,"diff_omics.rds"))
#   diff_omics(x)
# })
# samplelist<- reactive({
#   file="celldiff"
#   samplelist<-readRDS(file.path(file,"samplelist.rds"))
# return(samplelist)
# })
# data_rds_group<- reactive({
#   
#   data_rds_group<-readRDS(file.path("data_rds_group.rds"))
# return(data_rds_group)
# })
# x<-readRDS(file.path("diff_omics.rds"))
# diff_omics(x)
# 
# omics_count<- reactive({
#   
#   omics_count<-readRDS(file.path("omics_count.rds"))
#   return(omics_count)
# })
#################
omics_count<- eventReactive(c(input$start_diff_analysis), {
  req(data_rds_group())
  withProgress(message = "Processing data...",value=0.8,{
    data_rds_group <- data_rds_group()
    data_mrds_group<-data_rds_group[[1]]
    data_trds_group<-data_rds_group[[2]]
    source("./source/main_program/runcount.R")
    omics_count_m<-runcount(data_mrds_group,"sample")
    omics_count_t<-runcount(data_trds_group,"sample")
    count<-omics_count_m %>% 
      left_join(omics_count_t[, setdiff(colnames(omics_count_t), c("x", "y"))],by=c("x_y"="x_y"))
    countt<-count %>%
      dplyr::select(-x,-y)
    rownames(countt)<-countt$x_y
    #
    countt<-countt %>%
      dplyr::select(-x_y)
    
    count1<-as.data.frame(t(countt))
    count1$ID=rownames(count1)
    omics_count<-list(count,count1)
    return(omics_count)
  })
})

omics_corr_data1<-  reactive({
  req(diff_omics(), omics_count())
  withProgress(message = "Processing data...",value=0.8,{
    diff_omics<-diff_omics()
    omics_count<-omics_count()
    samplelist<-samplelist()
    combine_count<-omics_count[[1]]
    diff_m<-diff_omics[[1]]
    diff_t<-diff_omics[[2]]
    t_diff<- diff_t %>%
      filter(State!="Non-significant")
    req(nrow(t_diff) > 0)  
    m_diff<-diff_m %>%
      filter(State!="Non-significant")
    req(nrow(m_diff) > 0)  
    if(length(grep("Level",names(m_diff)))>0){
      m_diff<-subset(m_diff,m_diff$Level!="level5")
    }
    req(nrow(m_diff) > 0)  
    groupname<-c("treatment","control")
    source("./source/main_program/runcount.R")
    omics_corr_data<-omics_corr_processing1(m_diff,t_diff,combine_count,samplelist,groupname)
    return(omics_corr_data)
  })
})

omics_corr_data_save1<-  reactive({
  req(omics_corr_data1())
  withProgress(message = "Processing data...",value=0.8,{
    omics_corr_data<-omics_corr_data1()
    diff_m_cor<-omics_corr_data[[1]]
    diff_t_cor<-omics_corr_data[[2]]
    con_diff_m_cor<-omics_corr_data[[3]]
    con_diff_t_cor<-omics_corr_data[[4]]
    source("./source/CorrelationAnalysisFunction/Correlation/run_corr_diff_data.R")
    treat_omics_corr_data_save<-run_cordata(diff_m_cor,diff_t_cor)
    con_omics_corr_data_save<-run_cordata(con_diff_m_cor,con_diff_t_cor)
    omics_corr_data_save=list(treat_omics_corr_data_save,con_omics_corr_data_save)
    return(omics_corr_data_save)
  })
})


omics_corrResult1<- reactive({
  req(omics_corr_data_save1())
  withProgress(message = "Processing data...",value=0.8,{
    omics_corr_data_save<-omics_corr_data_save1()
    casedata<-omics_corr_data_save[[1]]
    controldata<-omics_corr_data_save[[2]]
    options(warn = -1)
    source("./source/CorrelationAnalysisFunction/Correlation/run_corr_diff_data.R")
    
    case_omics_corrResult<-run_corr_diff_data(casedata)
    con_omics_corrResult<-run_corr_diff_data(controldata) 
    
    omics_corrResult<-list(case_omics_corrResult,con_omics_corrResult)

    options(warn = 0)

    return(omics_corrResult)
  })
})


top_corrResultfilter<- reactive({
  omics_corrResult<-omics_corrResult1()
  req(omics_corrResult)
  withProgress(message = "Processing data...",value=0.8,{
    caseedge=omics_corrResult[[1]]
    controledge=omics_corrResult[[2]]
    
    omics_corr_data<-omics_corr_data1()
    m_mean=omics_corr_data[[5]]
    m_mean$Class="metabolite"
    t_mean=omics_corr_data[[6]]
    t_mean$Class="gene"
    c_mean=rbind(m_mean,t_mean)
    topnum=input$topdiffnum
    diff_omics<-diff_omics()
    
    diff_m<-diff_omics[[1]] %>%
      filter(State!="Non-significant") %>%
      arrange(`p_val_adj`)

    diff_t<-diff_omics[[2]]  %>%
      filter(State!="Non-significant") %>%
      arrange(`p_val_adj`)
    if(length(diff_m$metabolite)<topnum || length(diff_t$gene)<topnum){
      topnum=min(c(length(diff_m$metabolite),length(diff_t$gene)))
    }
    topm=as.character(diff_m$metabolite[1:topnum])
    topt=as.character(diff_t$gene[1:topnum])
    Rthreshold<-as.numeric(input$cor_Coefficient1)
    Pthreshold<-as.numeric(input$cor_p_value1)
    caseedge<-caseedge |>
      dplyr::filter(from %in% c(topm,topt)) |>
      dplyr::filter(to %in% c(topm,topt)) |>
      dplyr::filter(abs(cor) > Rthreshold & padj < Pthreshold) |>
      dplyr::arrange(abs(cor),padj)
    
    
    controledge<-controledge |>
      dplyr::filter(from %in% c(topm,topt)) |>
      dplyr::filter(to %in% c(topm,topt)) |>
      dplyr::filter(abs(cor) > Rthreshold & padj < Pthreshold) |>
      dplyr::arrange(abs(cor),padj)
    

    source("./source/CorrelationAnalysisFunction/Correlation/net_show.R")
    if(any(!is.na(caseedge))){
      caseedge <-run_compare_edge(caseedge)
    }else{
      caseedge <-NULL
      showNotification("After filtering, no significant correlations were found in the treatment group.")
    }
    if(any(!is.na(controledge))){
      controledge <-run_compare_edge(controledge)
    }else{
      controledge <-NULL
      showNotification("After filtering, no significant correlations were found in the control group.")
    }

    if(any(!is.na(caseedge)) | any(!is.na(controledge))){
      edgeboth<-data.frame(from_to=union(controledge$from_to,caseedge$from_to)) 
      if(any(!is.na(caseedge))){
        caseedge_U<-edgeboth |>
          dplyr::left_join(caseedge,by=c("from_to"="from_to")) |>
          dplyr::select(-from,-to) |>
          tidyr::separate(col = from_to, into = c("from", "to"), sep = "_",remove = FALSE) |>
          dplyr::select(-from_to)
        nodedata=data.frame(node=union(caseedge_U$from,caseedge_U$to)) |>
          dplyr::left_join(c_mean,by=c("node"="feature_ID"))
      }else{
        caseedge_U=NULL
      }
      if(any(!is.na(controledge))){
        controledge_U<- edgeboth |>
          dplyr::left_join(controledge,by=c("from_to"="from_to")) |>
          dplyr::select(-from,-to) |>
          tidyr::separate(col = from_to, into = c("from", "to"), sep = "_",remove = FALSE) |>
          dplyr::select(-from_to)
        
        nodedata=data.frame(node=union(controledge_U$from,controledge_U$to)) |>
          dplyr::left_join(c_mean,by=c("node"="feature_ID"))
      }else{
        controledge_U=NULL
      }
      top_corrResultfilter=list(caseedge_U,controledge_U,nodedata)
      return(top_corrResultfilter)
    }else{
      return(NULL)
    }
  })
})
diffnet_plotsave <- reactive({
  req(top_corrResultfilter())
  withProgress(message = "Processing data...",value=0.8,{
    top_corrResultfilter<-top_corrResultfilter()
    case_edges <- top_corrResultfilter[[1]]
    control_edges <- top_corrResultfilter[[2]]
    nodes <- top_corrResultfilter[[3]]

    source("./source/CorrelationAnalysisFunction/Correlation/net_show.R")
    layout_type<-input$layout_type
    color_type<-"Normalized mean"
    show_node_name<-input$show_node_name
    node_name_size<-input$node_name_size
    node_size<-input$node_size
    if(any(!is.na(case_edges))){
      case_igraph <- run_igraph(nodes=nodes,edges=case_edges)
      igraph::V(case_igraph)$`Normalized mean`<-igraph::V(case_igraph)$case_norm_mean
      print("case_p START")
      case_p<-net_show(igraph=case_igraph,layout_type=layout_type,
                       color_type=color_type,show_node_name=show_node_name,node_name_size=node_name_size,
                       node_size=node_size)
      print("case_p DONE")
    }else{
      case_p<-NULL
    }
    if(any(!is.na(control_edges))){
      control_igraph <- run_igraph(nodes=nodes,edges=control_edges)
      print("con_p START")
      igraph::V(control_igraph)$`Normalized mean`<-igraph::V(control_igraph)$control_norm_mean
      con_p<-net_show(igraph=control_igraph,layout_type=layout_type,
                      color_type=color_type,show_node_name=show_node_name,node_name_size=node_name_size,
                      node_size=node_size)
      print("con_p DONE")
    }else{
      con_p<-NULL
    }

    return(list(case_p,con_p))
  })
})
output$cornetwork_plot1 <- renderPlot({
  req( diffnet_plotsave()[[1]])
  withProgress(message = "Plotting...",value=0.8,{
    diffnet_plotsave()[[1]]
  })
})
output$cornetwork_plot2 <- renderPlot({
  req( diffnet_plotsave()[[2]])
  withProgress(message = "Plotting...",value=0.8,{
    diffnet_plotsave()[[2]]
  })
})






