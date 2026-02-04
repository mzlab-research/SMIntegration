suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(showtext))
suppressMessages(library(Cairo))
suppressMessages(library(patchwork))
# font_add("sans", regular = "arial.ttf", italic = "ariali.ttf")
# showtext_auto()
#2b+2c+2d---------------------------------------------------------------------------------------
#' @title Figure 2b, 2c, 2d: Spatial Clustering
#' @description Visualizes spatial domains identified by clustering.
#' 2b: Metabolomics-derived clusters.
#' 2c: Transcriptomics-derived clusters.
#' 2d: Integrated (Multi-omics) clusters.
#' @return Generates PDFs '2b.pdf', '2c.pdf', '2d.pdf'.
data<-readRDS("./Figures/rds/cluster.rds")
m=data[[1]]
t=data[[2]]
c=data[[3]]
#' @title Generate Spatial Cluster Plot
#' @description Creates a spatial scatter plot colored by cluster ID or cell type.
#' Uses a discrete color palette.
#' @param object Seurat object or data frame with coordinates.
#' @param features Column name to use for coloring (e.g., 'seurat_clusters').
#' @param pt.size Point size.
#' @return A ggplot object.
iPlot <- function(object, features, pt.size){
heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
if(features %in% c('seurat_clusters')){
  object$seurat_clusters <- factor(object$seurat_clusters,
                                             levels = as.character(sort(unique(as.integer(object$seurat_clusters)))))
}
plot <- ggplot(object, aes(x = x, y = y, color = !!sym(features))) +
  geom_point(shape = 19, size = pt.size) +
  xlim(min(object$x), max(object$x)) +
  ylim(min(object$y), max(object$y)) +
  coord_equal() +
  xlab("") + 
  ylab("") +
  theme_minimal(base_family = "sans") +  # 设置基础主题字体
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "plain"),
    text = element_text(family = "sans"),  # 全局字体设置
    legend.text = element_text(family = "sans"),
    legend.title = element_text(family = "sans")
  )

if(features %in% c('seurat_clusters')){
  k <- length(unique(object$seurat_clusters))
  plot <- plot + 
    scale_color_manual(values = cluster_Palette[1:k]) +
    guides(colour = guide_legend(
      title = "Group",
      override.aes = list(size = 3),
      nrow = 10
    ))
} else if(features %in% c('celltype')){
  k <- length(unique(object$celltype))
  plot <- plot + 
    scale_color_manual(values = cluster_Palette[1:k]) +
    guides(colour = guide_legend(
      title = "celltype",
      override.aes = list(size = 3),
      nrow = 10
    ))
}
return(plot)
}
mplot <- iPlot(m, feature = 'seurat_clusters', pt.size = 1) 
tplot <- iPlot(t, feature = 'seurat_clusters', pt.size = 1) 
cplot <- iPlot(c, feature = 'seurat_clusters', pt.size = 1) 
pdf("./Figures/2b.pdf", width = 6, height = 6)
print(mplot)
dev.off()
pdf("./Figures/2c.pdf", width = 6, height = 6)
print(tplot)
dev.off()
pdf("./Figures/2d.pdf", width = 6, height = 6)
print(cplot)
dev.off()
#2e+2f---------------------------------------------------------------------------------------
#' @title Figure 2e & 2f: Spatial Patterns
#' @description Visualizes the spatial patterns (modules) identified by SpaGene.
#' 2e: Metabolomics spatial patterns.
#' 2f: Transcriptomics spatial patterns.
#' @return Generates PDFs '2e.pdf' and '2f.pdf'.
find_spatial_pattern<-readRDS("./Figures/rds/find_spatial_pattern.rds")
data_m_parttern<-find_spatial_pattern[[1]]
data_t_parttern<-find_spatial_pattern[[2]]

#' @title Generate Spatial Pattern Grid
#' @description Creates a grid of spatial plots, one for each identified pattern.
#' Colors points by pattern weight/expression with an alpha transparency gradient.
#' @param pattern List containing pattern weights matrix.
#' @param location Data frame of coordinates.
#' @return A patchwork object combining individual pattern plots.
PlotPattern_c<-function (pattern, location, max.cutoff = 0.9, pt.size = 0.1, alpha.min = 0.1) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer")
  }
  colnames(location) <- c("x", "y")
  npattern <- dim(pattern$pattern)[1]
  plist <- list()
  for (i in 1:npattern) {
    feature = pattern$pattern[i, ]
    max.use <- quantile(feature, max.cutoff)
    feature[feature > max.use] <- max.use
    alpha = (feature - min(feature))/(max(feature) - min(feature)) * 
      (1 - alpha.min) + alpha.min
    tmp <- as.data.frame(cbind(location, exp = feature, alpha = alpha))
    p1 <- ggplot(tmp, aes(x = x, y = y, col = exp, alpha = alpha)) + 
      geom_point(size = pt.size) +
      scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 10, 
                                                                   name = "RdYlBu"))) + xlab("") + ylab("") +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank()) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, family = "sans"))+
      guides(color = "none", alpha = "none") + ggtitle(paste0("Pattern", 
                                                              i))+
      coord_equal()
    plist[[i]] <- p1
  }
  
  if(npattern==8){
    patchwork::wrap_plots(plist,ncol = 4)
  }else{
    patchwork::wrap_plots(plist)
  }
}
mplot2<-PlotPattern_c(data_m_parttern[[1]],data_m_parttern[[2]])
tplot2<-PlotPattern_c(data_t_parttern[[1]],data_t_parttern[[2]])
pdf("./Figures/2e.pdf", width =8, height = 7)
print(mplot2)
dev.off()
pdf("./Figures/2f.pdf", width = 8, height = 7)
print(tplot2)
dev.off()
#2g---------------------------------------------------------------------------------------
#' @title Figure 2g: Pattern Correlation Heatmap
#' @description Displays the correlation matrix between Metabolomics and Transcriptomics spatial patterns.
#' Uses `pheatmap` to cluster and visualize relationships.
#' @return Generates a PDF '2g.pdf'.
data_list<-readRDS("./Figures/rds/spatial_heatmap_data.rds")
cmt <- data_list[[1]]
annotation_col <- data_list[[2]]
devs <- dev.list()
current_dev <- dev.cur()
if (!is.null(devs) && length(devs) > 1) {
  for (d in names(devs)[names(devs) != names(current_dev)]) {
    tryCatch({
      dev.off(dev.list()[[d]])
    }, error = function(e) {})
  }
}

breaks <- seq(-1, 1, length.out = 51)

plot_obj <- pheatmap::pheatmap(
  cmt,
  scale = "none",
  cluster_row = TRUE,
  cluster_col = TRUE,
  border = NA,
  annotation_col = annotation_col,
  fontsize_number = 12,
  number_color = "green",
  cellwidth = 20,
  cellheight = 20,
  color = colorRampPalette(c("#8470FF", "#FFFFFF", "#FF0000"))(50),
  breaks = breaks,  
  fontfamily = "sans",  
  treeheight_row = 50,  
  treeheight_col = 50,  
  silent = TRUE  
)


heatmap_gg <- ggplot() +
  annotation_custom(plot_obj$gtable) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) 

CairoPDF("./Figures/2g.pdf", width = 8, height = 7, family = "sans")
grid::grid.newpage()
grid::grid.draw(plot_obj$gtable)
dev.off()
