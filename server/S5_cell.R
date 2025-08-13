
output$download_cell_annotation_demo <- downloadHandler(
  filename = function() {
    paste0("cell_annotation_demo.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    DemoData <- fread("./example_data/cell_annotation_demo.txt")
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  }
)
output$download_cell_annotation_data <- downloadHandler(
  filename = function() {
    paste0("cell_annotation.txt")
  },
  content = function(file) {
    cell_annotation_data <- cell_annotation_data()
     withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData<-data.frame(x=cell_annotation_data$x,y=cell_annotation_data$y,celltype=cell_annotation_data$celltype)
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  }
)

output$download_cell_annotation_plot <- downloadHandler(
  filename = function() {
    "cell_annotation.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- cell_annotation_plot_save()

    ggsave(file, plot = p, device = "png", width = 8, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })

submitted_data <- reactiveValues(names = list())


observeEvent(input$start_cell_annotation, {
  t_cluster_rds<-cluster_plotlist()[[2]][[3]]
  t_cluster_txt<-data.frame(seurat_clusters=t_cluster_rds$seurat_clusters)

  submitted_data$names <- lapply(unique(t_cluster_txt$seurat_clusters), function(group) {
    input[[paste0("group_", group)]]
  })
  
})
cell_annotation_data<-eventReactive(input$start_cell_annotation, {
  req(cluster_plotlist())
  withProgress(message = "Processing data...",value=0.8,{
    t_cluster_rds<-cluster_plotlist()[[2]][[3]]

    if (input$demo_select == "Use demo data"){
      data <- fread("./example_data/cell_annotation_demo.txt")
      data$x_y<-paste0(data$x,"_",data$y)
      t_cluster_rds$x_y<-paste0(t_cluster_rds$x,"_",t_cluster_rds$y)
      t_cluster_rds$celltype<-t_cluster_rds$seurat_clusters
      t_cluster_txt<-data.frame(x_y=t_cluster_rds$x_y,celltype=t_cluster_rds$celltype)
      t_cluster_txt$celltype <- ifelse(t_cluster_txt$x_y %in% data$x_y, 
                                       sapply(t_cluster_txt$x_y, function(x_y) {
                                         data$celltype[match(x_y, data$x_y)]
                                       }), 
                                       t_cluster_txt$celltype)
      t_cluster_rds$celltype<-t_cluster_txt$celltype
      cell_annotation_data<-t_cluster_rds
      return(cell_annotation_data)
  }
    if (input$cell_annotation_select=="Upload_annotation_table"){
    if(!is.null(input$cellannotationfile) && input$cellannotationfile$name != ""){
      data <- fread(input$cellannotationfile$datapath) 
      if("x" %in% colnames(data) && "y" %in% colnames(data) && "celltype" %in% colnames(data)){
           data$x_y<-paste0(data$x,"_",data$y)
           t_cluster_rds$x_y<-paste0(t_cluster_rds$x,"_",t_cluster_rds$y)
           t_cluster_rds$celltype<-t_cluster_rds$seurat_clusters
           t_cluster_txt<-data.frame(x_y=t_cluster_rds$x_y,celltype=t_cluster_rds$celltype)
           t_cluster_txt$celltype <- ifelse(t_cluster_txt$x_y %in% data$x_y, 
                                            sapply(t_cluster_txt$x_y, function(x_y) {
                                              data$celltype[match(x_y, data$x_y)]
                                            }), 
                                            t_cluster_txt$celltype)
           t_cluster_rds$celltype<-t_cluster_txt$celltype
           cell_annotation_data<-t_cluster_rds
           return(cell_annotation_data)
      }else if("index" %in% colnames(data) && "celltype" %in% colnames(data)){
        data$index<-as.character(data$index)
        t_cluster_txt<-data.frame(index=as.character(t_cluster_rds$index),celltype=t_cluster_rds$celltype)
        t_cluster_txt$celltype <- ifelse(t_cluster_txt$index %in% data$index, 
                                         sapply(t_cluster_txt$index, function(index) {
                                           data$celltype[match(index, data$index)]
                                         }), 
                                         t_cluster_txt$celltype)
        t_cluster_rds$celltype<-t_cluster_txt$celltype
        cell_annotation_data<-t_cluster_rds
        return(cell_annotation_data)
      }else{
        showNotification("Please check if the column names meet the requirements.")
        return(NULL)
      }
    }
    }else{
      nomaldata=GetAssayData(t_cluster_rds@assays$SCT,slot = "scale.data")
      clusters=t_cluster_rds@meta.data$seurat_clusters

      if(input$speciesname_select=="mmu"){
        load("./source/Database/SingleR/mouse_all.rData")
      }else if(input$speciesname_select=="hsa"){
        load("./source/Database/SingleR/human_all.rData")
        }

      cellpred=SingleR(test = nomaldata,ref = refdata,labels = refdata$label.main,clusters = clusters,assay.type.test = "logcounts",assay.type.ref = "logcounts")
      celltype=data.frame(ClusterID=rownames(cellpred),celltype=cellpred$labels,stringsAsFactors = F)

      t_cluster_rds@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(t_cluster_rds@meta.data$seurat_clusters)), from =celltype$ClusterID , to = celltype$celltype)

      cell_annotation_data<-t_cluster_rds
      return(cell_annotation_data)
    }
   
  })
})
cell_annotation_rds<-eventReactive(input$start_cell_annotation, {
  req(cluster_plotlist(),cell_annotation_data())
  withProgress(message = "Processing data...",value=0.8,{
    cell_annotation_data<-cell_annotation_data()
    cluster_rds <- cluster_plotlist()
    cluster_rds_m<-cluster_rds[[1]][[3]]
    cluster_rds_c<-cluster_rds[[3]][[3]]

    celllist<-data.frame(x_y=paste0(cell_annotation_data$x,"_",cell_annotation_data$y),
                         celltype=cell_annotation_data$celltype)
    celllist_m<-data.frame(x_y=paste0(cluster_rds_m$x,"_",cluster_rds_m$y)) %>%
      left_join(celllist,by=c("x_y"="x_y"))
    cluster_rds_m$celltype<-celllist_m$celltype
    
    celllist_c<-data.frame(x_y=paste0(cluster_rds_c$x,"_",cluster_rds_c$y)) %>%
      left_join(celllist,by=c("x_y"="x_y"))
    cluster_rds_c$celltype<-celllist_c$celltype

    source("./source/main_program/slim_seurat.R")
    cluster_rds_m <- slim_seurat(seurat_obj=cluster_rds_m)
    cell_annotation_data <- slim_seurat(seurat_obj=cell_annotation_data)
    cluster_rds_c <- slim_seurat(seurat_obj=cluster_rds_c)
    gc()
    cell_annotation_rds<-list(cluster_rds_m,cell_annotation_data,cluster_rds_c)

    return(cell_annotation_rds)

  })
})

cell_annotation_plot_save<-eventReactive(input$start_cell_annotation, {
  req(cell_annotation_data())
  cell_annotation_data<-cell_annotation_data()
  source("./source/OverallAnalysisFunction/Clustering/clusterplot.R")
  overallplot <- function(obj,pointSize,breakseq) {
     plot1 <- Preview(obj, 'celltype',pointSize,breakseq)+
       labs(title="")+
       theme(plot.title = element_text(hjust = 0.5))
    return(plot1)
  }
  if(length(unique(cell_annotation_data$celltype))>30){
    pointSize=0.5
  }else{
    pointSize=1
  }
  cell_annotation_plot_save<-overallplot(cell_annotation_data,pointSize=pointSize,breakseq=50)
  return(cell_annotation_plot_save)
  })


output$cell_annotation_plot <- renderPlot({
req(cell_annotation_plot_save())
  withProgress(message = "Plotting...",value=0.8,{
    
cell_annotation_plot_save()
    

  })
})




output$cell_annotation_button_container <- renderUI({
  
   if (input$cell_annotation_select=="Upload_annotation_table"){
    tagList(
      p("The cell annotation file (.txt) must contain the following columns: x, y, celltype"),
      fileInput("cellannotationfile", "Upload the cell annotation file",
                accept = ".txt"),
      downloadButton("download_cell_annotation_demo", "Download demo cell annotation file")
      
    )
    
  }else{
    NULL
  }
})

