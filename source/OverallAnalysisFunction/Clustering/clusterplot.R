

Preview <- function(object, feature,pointSize=0.2,breakseq=50){
  patch <- iPlot(object, features = feature, pt.size = pointSize,breakseq=breakseq)
  return(patch)
}
iPlot <- function(object, features, pt.size,breakseq){
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  if(features %in% c('seurat_clusters')){
    object@meta.data$seurat_clusters <- factor(object@meta.data$seurat_clusters,
                                               levels = as.character(sort(unique(as.integer(object@meta.data$seurat_clusters)))))
  }
  plot <- ggplot(object@meta.data, aes(x = x, y = y, color =  !!sym(features))) +
    geom_point(shape = 19, size = pt.size) +
    xlim(min(object@meta.data$x),max(object@meta.data$x))+
    ylim(min(object@meta.data$y),max(object@meta.data$y))+
    coord_equal()+
    xlab("")+ylab("")+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank()) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "plain"))

  if (features %in% c('nCount_Spatial','nCount_SCT')){
    plot <- ggplot(object@meta.data, aes(x = x, y = y, color =  !!sym(features))) +
      geom_point(shape = 19, size = pt.size) +
      scale_color_gradientn(colours = heatmap_Palette(100)) +
      guides(colour = guide_colorbar(title = "") )+
      xlim(min(object@meta.data$x),max(object@meta.data$x))+
      ylim(min(object@meta.data$y),max(object@meta.data$y))+
      coord_equal()+
      xlab("")+ylab("")+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank()) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
  }else if(features %in% c('nFeature_Spatial','nFeature_SCT')){
    min_val <- min(object@meta.data[features], na.rm = TRUE)
    max_val <- max(object@meta.data[features], na.rm = TRUE)
    interval <- round((max_val - min_val) / 5)
    interval <- ifelse(interval %% 1 == 0, interval, floor(interval))
    breaks <- seq(floor(min_val), ceiling(max_val), by = interval)
    
    plot <- ggplot(object@meta.data, aes(x = x, y = y, color =  !!sym(features))) +# xlab("y")+ ylab("x")+
      geom_point(shape = 19, size = pt.size) +
      scale_color_gradientn(colours = heatmap_Palette(100), breaks = breaks) +
      guides(colour = guide_colorbar(title = "") )+
      xlim(min(object@meta.data$x),max(object@meta.data$x))+
      ylim(min(object@meta.data$y),max(object@meta.data$y))+
      coord_equal()+
      xlab("")+ylab("")+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank()) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
  }else if(features %in% c('seurat_clusters')){
    
    k <- length(unique(object$seurat_clusters))
    plot <- plot + scale_color_manual(values = cluster_Palette[1:k]) +
      guides(colour = guide_legend(title = "Group",override.aes = list(size=3), nrow = 10))
  }else if(features %in% c('celltype')){
    k <- length(unique(object$celltype))
    plot <- plot + scale_color_manual(values = cluster_Palette[1:k]) +
      guides(colour = guide_legend(title = "celltype",override.aes = list(size=3), nrow = 10))
  }
  return(plot)
}

Clustering <- function(object,resolution,clustertype, dims = 30){
 options(warn = -1)
 object <- RunPCA(object)
 options(warn = 0)
  object <- FindNeighbors(object, dims = 1:dims)
  object <- FindClusters(object, algorithm=clustertype,verbose = FALSE, resolution = resolution)
  object$index<-as.character(object$x* 2^32 + object$y)
  
  return(object)
}
umap_plot<-function(obj,pt.size = 0.2){
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

  p<-DimPlot(object = obj,group.by="groups",reduction="umap",
          cols = cluster_Palette, pt.size=pt.size,ncol=1, label=F,
          label.size=2) +
    theme(legend.position="right",
          axis.title.x=element_text(size=8),
          axis.title.y=element_text(size=8),
          axis.text.x=element_text(vjust=1,size=8),
          axis.text.y=element_text(vjust=1,size=8),
          plot.title = element_text(hjust = 0.5))
  return(p)
}


run_clusterplot <- function(obj,type,pointSize,breakseq,resolution,clustertype,umapk=7) {
  if(clustertype==0){
    clustertype=1
    umapk=resolution
    resolution=0.5
    obj <- Clustering(object=obj,resolution=resolution,clustertype=clustertype,dims = 30)
    obj <- RunUMAP(obj, dims = 1:30)
    umap = obj@reductions$umap@cell.embeddings %>%
      as.data.frame()
    kmeans_clusters1<- kmeans(umap, centers = umapk)
    Group_cluster<-as.character(kmeans_clusters1$cluster)
    clusterdata<-umap %>% cbind(Group_cluster,y=obj@meta.data[["y"]],x=obj@meta.data[["x"]])
  }else{
    
    obj <- Clustering(object=obj,resolution=resolution,clustertype=clustertype,dims = 30)
    Group_cluster<-as.character(obj@active.ident)
    clusterdata<-data.frame(Group_cluster,y=obj@meta.data[["y"]],x=obj@meta.data[["x"]])
  }
  obj@meta.data$seurat_clusters<-Group_cluster
  plot1 <- iPlot(obj, feature = 'seurat_clusters', pt.size = pointSize,breakseq) +
    labs(title=paste0(type))
  plotlist<-list(plot1,clusterdata,obj)
  return(plotlist)
}
