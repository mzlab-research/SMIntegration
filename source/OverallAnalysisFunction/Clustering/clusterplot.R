

#' @title Preview Spatial Plot
#' @description Wrapper for iPlot to generate a spatial plot of a specific feature.
#' @param object Seurat object.
#' @param feature Character. Name of the feature or metadata column to plot.
#' @param pointSize Numeric. Size of the points.
#' @param breakseq Numeric. Break sequence for color scale (unused in current iPlot).
#' @return A ggplot object.
Preview <- function(object, feature,pointSize=0.2,breakseq=50){
  patch <- iPlot(object, features = feature, pt.size = pointSize,breakseq=breakseq)
  return(patch)
}

#' @title Internal Spatial Plotting Function
#' @description Generates spatial scatter plots for features, clusters, or metadata.
#' Handles color palettes for discrete (clusters) and continuous (counts) data.
#' @param object Seurat object.
#' @param features Character. Feature to plot.
#' @param pt.size Numeric. Point size.
#' @param breakseq Numeric. Unused parameter.
#' @return A ggplot object.
iPlot <- function(object, features, pt.size,breakseq){
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  
  # Ensure cluster labels are sorted factors
  if(features %in% c('seurat_clusters')){
    object@meta.data$seurat_clusters <- factor(object@meta.data$seurat_clusters,
                                               levels = as.character(sort(unique(as.integer(object@meta.data$seurat_clusters)))))
  }
  
  # Base plot structure
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

  # Customize based on feature type (Continuous vs Discrete)
  if (features %in% c('nCount_Spatial','nCount_SCT')){
    # Continuous scale for counts
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
    # Continuous scale with breaks for features
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
    # Discrete scale for clusters
    k <- length(unique(object$seurat_clusters))
    plot <- plot + scale_color_manual(values = cluster_Palette[1:k]) +
      guides(colour = guide_legend(title = "Group",override.aes = list(size=3), nrow = 10))
  }else if(features %in% c('celltype')){
    # Discrete scale for cell types
    k <- length(unique(object$celltype))
    plot <- plot + scale_color_manual(values = cluster_Palette[1:k]) +
      guides(colour = guide_legend(title = "celltype",override.aes = list(size=3), nrow = 10))
  }
  return(plot)
}

#' @title Perform Graph-based Clustering
#' @description Runs PCA, constructs SNN graph, and performs modularity optimization (Louvain/SLM).
#' @param object Seurat object.
#' @param resolution Numeric. Clustering resolution parameter.
#' @param clustertype Integer. Algorithm type (1=Louvain, 2=Louvain with multilevel refinement, 3=SLM).
#' @param dims Integer. Number of PCs to use.
#' @return Seurat object with cluster assignments.
Clustering <- function(object,resolution,clustertype, dims = 50){
 options(warn = -1)
 object <- RunPCA(object)
 options(warn = 0)
  object <- FindNeighbors(object, dims = 1:dims)
  object <- FindClusters(object, algorithm=clustertype,verbose = FALSE, resolution = resolution)
  object$index<-as.character(object$x* 2^32 + object$y)
  
  return(object)
}

#' @title UMAP Plot Wrapper
#' @description Generates a UMAP plot grouped by 'groups'.
#' @param obj Seurat object.
#' @param pt.size Numeric. Point size.
#' @return A ggplot object.
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


#' @title Run Clustering and Generate Plots
#' @description Main controller for clustering analysis. Supports multiple algorithms:
#' PCA-Kmeans (Type 4), UMAP-Kmeans (Type 0), and Graph-based (Louvain/SLM).
#' @param obj Seurat object.
#' @param type Character. Data type label (e.g., "Metabolite").
#' @param pointSize Numeric. Point size for spatial plot.
#' @param breakseq Numeric. Unused.
#' @param resolution Numeric. Resolution for graph-based or K for K-means.
#' @param clustertype Integer. Algorithm code.
#' @param umapk Integer. Redundant parameter (uses resolution as K).
#' @return A list containing: [1] Spatial Plot, [2] Cluster Data (coordinates + ID), [3] Updated Seurat Object.
run_clusterplot <- function(obj,type,pointSize,breakseq,resolution,clustertype,umapk=7) {
  if (clustertype == 4) {
    # PCA-kmeans logic: K-means clustering on PCA embeddings
    pcak = resolution
    # Run PCA if not already present or ensure enough dims
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    obj$index <- as.character(obj$x * 2^32 + obj$y)
    
    pca_embed <- Embeddings(obj, reduction = "pca")[, 1:30] %>% as.data.frame()
    
    kmeans_clusters1 <- kmeans(pca_embed, centers = pcak)
    Group_cluster <- as.character(kmeans_clusters1$cluster)
    clusterdata <- pca_embed %>% cbind(Group_cluster, y = obj@meta.data[["y"]], x = obj@meta.data[["x"]])
    
  } else if(clustertype==0){
    # UMAP-kmeans logic: K-means clustering on UMAP embeddings
    clustertype=1
    umapk=resolution # Treat resolution input as K for K-means
    resolution=0.5 # Default resolution for intermediate graph clustering (if needed)
    obj <- Clustering(object=obj,resolution=resolution,clustertype=clustertype,dims = 30)
    obj <- RunUMAP(obj, dims = 1:30)
    umap = obj@reductions$umap@cell.embeddings %>%
      as.data.frame()
    kmeans_clusters1<- kmeans(umap, centers = umapk)
    Group_cluster<-as.character(kmeans_clusters1$cluster)
    clusterdata<-umap %>% cbind(Group_cluster,y=obj@meta.data[["y"]],x=obj@meta.data[["x"]])
  }else{
    # Standard Graph-based Clustering (Louvain, SLM, etc.)
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
