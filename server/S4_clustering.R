
output$download_cluster_plot  <- downloadHandler(
  filename = function() {
    "cluster_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- get("cluster_plot", envir = .GlobalEnv)

    ggsave(file, plot = p, device = "png", width = 24, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })
output$download_cluster_data <- downloadHandler(
  filename = function() {
    "cluster_data.zip"
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
    cluster_plotlist<-cluster_plotlist()
    tempdir <- setwd(tempdir())
    on.exit(setwd(tempdir))
    fi=c("Metabolite_cluster_data.csv","Gene_cluster_data.csv","Merge_cluster_data.csv")
    for (i in 1:length(fi)) {
      data<-cluster_plotlist[[i]]
      write.csv(data[[2]], fi[i], row.names = FALSE)
    }

zip(file,fi)
    })
  })

output$download_cluster_sanky_plot <- downloadHandler(
  filename = function() {
    "cluster_sanky_plot.zip"
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
      p <- cluster_sanky_plot_save()
      tempdir <- setwd(tempdir())
      on.exit(setwd(tempdir))
      fi=c("sankey_plot.html")
      saveNetwork(p,fi[1])
      zip(file,fi)
      
    })
  })
output$download_cluster_sanky_plot_data <- downloadHandler(
  filename = function() {
    paste0("cluster_sanky_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData <- cluster_sanky()
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })



rds_norm<- eventReactive(c(input$Submit), {
  req(data_rds())#data_crds()
  withProgress(message = "Processing data...",value=0.8,{
    data_rds<-data_rds()
    data_mrds <- RUNSCT(data_rds[[1]])
    data_trds <- RUNSCT(data_rds[[2]])

    decon_mtrx = data_mrds@assays$Spatial$counts
    decon_ttrx = data_trds@assays$Spatial$counts
    source("./source/preprocessing/pre_merge.R")
    data_crds<-run_merge_rds(decon_mtrx,decon_ttrx)
    source("./source/preprocessing/run_prerds.R")
    data_crds<-run_prerds(data=data_crds)
    data_crds=RUNSCT(data_crds)
    data_crds@assays$SCT$data=rbind(
      data_mrds@assays$SCT$data,
      data_trds@assays$SCT$data
    )
    data_crds@assays$SCT$scale.data=rbind(
      data_mrds@assays$SCT$scale.data,
      data_trds@assays$SCT$scale.data
    )
    # 提取高变特征
    hv_metab <- VariableFeatures(data_mrds)
    hv_trans <- VariableFeatures(data_trds)
    hv_combined <- c(hv_metab, hv_trans)
    VariableFeatures(data_crds) <- hv_combined

    rds_norm=list(data_mrds=data_mrds,data_trds=data_trds,combine_rds=data_crds)

    return(rds_norm)
  })
})


cluster_plotlist <-eventReactive(c(input$start_cluster), {
  req(rds_norm())
  withProgress(message = "Processing data...",value=0.8,{
    if(!is.null(input$cluster_resolution)){
      rds_norm<-rds_norm()
      source("./source/OverallAnalysisFunction/Clustering/clusterplot.R")
      data_mrds=rds_norm[["data_mrds"]]
      data_trds=rds_norm[["data_trds"]]
      combine_rds=rds_norm[["combine_rds"]]
      if(input$cluster_select=="UMAP-kmeans"){
        clustertype=0
        k=7
      }else if(input$cluster_select=="LV"){
        clustertype=1
      }else if(input$cluster_select=="LM"){
        clustertype=2
      }else{
        clustertype=3
      }
      resolution=input$cluster_resolution
      mplot<-run_clusterplot(data_mrds,type="Metabolite",pointSize=1,breakseq=50,resolution=resolution,clustertype)
      tplot<-run_clusterplot(data_trds,type="Gene",pointSize=1,breakseq=50,resolution=resolution,clustertype)
      cplot<-run_clusterplot(combine_rds,type="Merge",pointSize=1,breakseq=50,resolution=resolution,clustertype)
      plotlist<-list(mplot,tplot,cplot)
      return(plotlist)
      gc()
    }else{
      return(NULL)
    }
    
  })
})

observeEvent(c(input$start_cluster),{
  req(cluster_plotlist())
  cluster_plotlist<-cluster_plotlist()
  mplot<-cluster_plotlist[[1]][[1]]
  tplot<-cluster_plotlist[[2]][[1]]
  cplot<-cluster_plotlist[[3]][[1]]

  output$cluster_plot <- renderPlot({
    withProgress(message = "Plotting...",value=0.8,{
    p<-gridExtra::grid.arrange(mplot,tplot,cplot,ncol=3)
    p
    assign("cluster_plot", p, envir = .GlobalEnv)
    })
  })
})

cluster_sanky <-eventReactive(c(input$start_cluster), {
  req(cluster_plotlist())
  withProgress(message = "Processing data...",value=0.8,{
    cluster_plotlist<-cluster_plotlist()
    mcluster<-cluster_plotlist[[1]][[3]]
    tcluster<-cluster_plotlist[[2]][[3]]
    ccluster<-cluster_plotlist[[3]][[3]]
    source("./source/main_program/slim_seurat.R")
    mcluster <- slim_seurat(seurat_obj=mcluster)
    tcluster <- slim_seurat(seurat_obj=tcluster)
    ccluster <- slim_seurat(seurat_obj=ccluster)
    gc()
    mcluster_data<-data.frame(x_y=mcluster$x_y,clusters_m=paste0("Metabolite_",mcluster$seurat_clusters))
    tcluster_data<-data.frame(x_y=tcluster$x_y,clusters_t=paste0("Gene_",tcluster$seurat_clusters))
    ccluster_data<-data.frame(x_y=ccluster$x_y,clusters_c=paste0("Merge_",ccluster$seurat_clusters))
    cluster_data<-mcluster_data %>%
      dplyr::left_join(tcluster_data,by=c("x_y"="x_y")) %>%
      dplyr::left_join(ccluster_data,by=c("x_y"="x_y"))
    cluster_data$value<-1
    ##part 1
    cluster_data1<-cluster_data %>%
      dplyr::select(x_y,clusters_m,clusters_c,value)
    colnames(cluster_data1)<-c("x_y","source","target","value")
    cluster_data1$color<-cluster_data1$target
    #part 2
    cluster_data2<-cluster_data %>%
      dplyr::select(x_y,clusters_c,clusters_t,value)
    colnames(cluster_data2)<-c("x_y","source","target","value")
    cluster_data2$color<-cluster_data2$source
    cluster_sanky<-rbind(cluster_data1,cluster_data2)
    return(cluster_sanky)
  })
})
cluster_sanky_plot_save <- eventReactive(c(input$start_cluster), {
  req(cluster_sanky())
  withProgress(message = "Processing data...",value=0.8,{
    source("./source/OverallAnalysisFunction/Sanky/sankyplot.R")
    links<-cluster_sanky()
    p=sankyplot_ID(links)
    
    return(p)
  })
})
observeEvent(input$start_cluster, {
  output$cluster_sanky_plot <- renderUI({
    req( cluster_sanky_plot_save())
    withProgress(message = "Plotting...", value = 0.8, {
      cluster_sanky_plot_save()
    })
  })
})


output$culster_resolution_button_container <- renderUI({
  if (input$cluster_select=="UMAP-kmeans") {
    sliderInput("cluster_resolution","Number of Clusters",value = 7,min=2,max=25,step = 1)
  } else {
    sliderInput("cluster_resolution","Resolution of clusters",value =2,min=0.1,max=2,step = 0.1)
  }
})
