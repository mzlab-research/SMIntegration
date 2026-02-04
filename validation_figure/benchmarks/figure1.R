suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(showtext))
suppressMessages(library(Cairo))
suppressMessages(library(networkD3))
suppressMessages(library(webshot))
suppressMessages(library(htmlwidgets))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(ggtext))
suppressMessages(library(VennDiagram))
suppressMessages(library(patchwork))
suppressMessages(library(gridExtra))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(magick))
suppressMessages(library(stringr))

# font_add("Arial", regular = "arial.ttf", italic = "ariali.ttf")
# showtext_auto()
######################################################1a
#' @title Figure 1a: Total Counts Distribution
#' @description Visualizes the spatial distribution of total ion counts (Metabolomics) and total gene counts (Transcriptomics).
#' @return Generates PNGs '1a-1.png' (Metabolomics) and '1a-2.png' (Transcriptomics).
data<-readRDS("./Figures/rds/cluster.rds")
m=data[[1]]
t=data[[2]]
c=data[[3]]
#' @title Generate Spatial Scatter Plot
#' @description Creates a spatial scatter plot colored by a specific feature (counts, clusters, or cell types).
#' Supports continuous (gradient) and discrete (categorical) color scales.
#' @param object Seurat object or data frame with coordinates.
#' @param features Column name to visualize.
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
    theme_minimal(base_family = "Arial") +  
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "plain"),
      text = element_text(family = "Arial"),  
      legend.text = element_text(family = "Arial"),
      legend.title = element_text(family = "Arial")
    )
  if (features %in% c('nCount_Spatial','nCount_SCT')){
    plot <- ggplot(object, aes(x = x, y = y, color =  !!sym(features))) +
      geom_point(shape = 19, size = pt.size) +
      scale_color_gradientn(colours = heatmap_Palette(100)) +
      guides(colour = guide_colorbar(title = "Total Feature Abundance") )+
      xlim(min(object$x),max(object$x))+
      ylim(min(object$y),max(object$y))+
      coord_equal()+
      xlab("")+ylab("")+
      theme_minimal(base_family = "Arial") +  
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
        title = "Cell type",
        override.aes = list(size = 3),
        nrow = 10
      ))
  }
  return(plot)
}
mplot <-iPlot(m, 'nCount_Spatial',pt.size = 1) 
tplot <-iPlot(t, 'nCount_Spatial',pt.size = 1) 
png("./Figures/1a-1.png", width = 400, height = 400)
print(mplot)
dev.off()
png("./Figures/1a-2.png", width = 400, height = 400)
print(tplot)
dev.off()

#' @title Figure 1a-3: RGB Overlay
#' @description Registration Effectiveness Reflected by Spatial Overlay of Omics Data.
rgb_overlay=readRDS("./Figures/rds/rgb.rds")
png("./Figures/1a-3.png", width = 400, height = 400)
grid.raster(rgb_overlay)
dev.off()
######################################################1c
#' @title Figure 1c: Clustering and Correspondence
#' @description Visualizes spatial clusters for Metabolomics and Transcriptomics,
#' along with a PCA scree plot showing the variance explained by each principal component.
#' Also generates a Sankey diagram ('1c-3.png') showing the relationship between clusters across modalities.
#' @return Generates PNGs '1c-0.png', '1c-1.png' , '1c-2.png', and '1c-3.png'.

PCA=readRDS("./Figures/rds/PCA.rds")
PCAPLOT <- ggplot(PCA, aes(x = PC, y = Stdev)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  theme(panel.grid=element_blank())+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Merge Variance", y = "Standard Deviation", x = "PC") +
  theme(plot.title = element_text(hjust = 0.5))
png("./Figures/1c-0.png", width = 400, height = 400)
print(PCAPLOT)
dev.off()

mplot1 <- iPlot(m, feature = 'seurat_clusters', pt.size = 1.5) 
tplot1 <- iPlot(t, feature = 'seurat_clusters', pt.size = 1.5) 
png("./Figures/1c-1.png", width = 400, height = 400)
print(mplot1)
dev.off()
png("./Figures/1c-2.png", width = 400, height = 400)
print(tplot1)
dev.off()
mcluster_data<-data.frame(x_y=m$x_y,clusters_m=paste0("Metabolite_",m$seurat_clusters))
tcluster_data<-data.frame(x_y=paste0("x",t$x,"_y",t$y),clusters_t=paste0("Gene_",t$seurat_clusters))
ccluster_data<-data.frame(x_y=c$x_y,clusters_c=paste0("Merge_",c$seurat_clusters))
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
#' @title Generate Sankey Diagram
#' @description Creates a Sankey diagram using `networkD3` to visualize flows between source and target nodes.
#' @param links Data frame containing 'source', 'target', 'value', and 'color' columns.
#' @return A networkD3 Sankey plot object.
sankyplot_ID=function(links){
  nodes <- data.frame(
    name=c(as.character(links$source),
           as.character(links$target)) %>% unique()
  )
  links$IDsource <- match(links$source, nodes$name)-1
  links$IDtarget <- match(links$target, nodes$name)-1
  p<-sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",
                   LinkGroup = "color", 
                   sinksRight=FALSE,nodeWidth = 10, 
                   nodePadding = 4) 
  return(p)
}
cplot=sankyplot_ID(cluster_sanky)
saveWidget(cplot, "temp_sankey.html")
webshot("temp_sankey.html", "./Figures/1c-3.png", vwidth = 600, vheight = 400,zoom = 2) 
unlink("temp_sankey.html")
######################################################1b
#' @title Figure 1b: Spatial Patterns and Correlation
#' @description Visualizes spatial patterns identified by SpaGene.
#' Also generates a heatmap ('1b-3.png') of cross-omics pattern correlations.
#' @return Generates PNGs '1b-1.png', '1b-2.png', and '1b-3.png'.
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
            plot.title = element_text(hjust = 0.5, family = "Arial"))+
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

png("./Figures/1b-1.png", width = 800, height = 400)
print(mplot2)
dev.off()
png("./Figures/1b-2.png", width = 800, height = 400)
print(tplot2)
dev.off()

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
  fontsize_number = 6,
  number_color = "green",
  cellwidth = 14,
  cellheight = 14,
  color = colorRampPalette(c("#8470FF", "#FFFFFF", "#FF0000"))(50),
  breaks = breaks,  
  fontfamily = "Arial",  
  treeheight_row = 25,  
  treeheight_col = 25,  
  silent = TRUE  
)


heatmap_gg <- ggplot() +
  annotation_custom(plot_obj$gtable) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) 

png("./Figures/1b-3.png", width = 600, height = 400, family = "Arial")
grid::grid.newpage()
grid::grid.draw(plot_obj$gtable)
dev.off()

#################################################1d
#' @title Figure 1d: Cell Type Annotation
#' @description Visualizes cell type annotations derived from SingleR and manual table mapping.
#' @return Generates PNGs '1d-1.png' (SingleR) and '1d-2.png' (Table).
c_table<-readRDS("./Figures/rds/cell_annotation_rds.rds")
c_singler<-readRDS("./Figures/rds/cell_annotation_rds1.rds")
cellplot1 <-iPlot(c_singler, 'celltype',pt.size = 1.4) 
png("./Figures/1d-1.png", width = 400, height = 400)
print(cellplot1)
dev.off()
cellplot2 <-iPlot(c_table, 'celltype',pt.size = 1.4) 
  
png("./Figures/1d-2.png", width = 800, height = 400)
print(cellplot2)
dev.off()

#################################1e
#' @title Figure 1e: Differential Analysis Summary
#' @description Generates a composite view of the differential analysis workflow:
#' 1. Spatial Group Assignment (Treatment vs Control) - '1e-1.png'.
#' 2. UMAP Visualization of groups - '1e-2.png'.
#' 3. Volcano Plot of differential metabolites - '1e-3.png'.
#' @return Generates PNGs '1e-1.png', '1e-2.png', '1e-3.png'.
data=readRDS("./Figures/rds/NA-vs-MO/plotdata_group.rds")
data <- data |>
  group_by(groups) 
group_color<-c("treatment"="red","control"="blue","0"="grey")
p1 <-  ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = groups), size = 1.5) +
  guides(colour = guide_legend(
    title = "Group",
    override.aes = list(size = 3), 
    nrow = 10
  )) +
  scale_color_manual(
    values = group_color,
    labels = c("0" = "Other", 
               "treatment" = "NA", 
               "control" = "MO")
  ) +
  xlab("") + ylab("") +
  coord_equal() +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_text(family = "Arial"),
    legend.text = element_markdown(family = "Arial"),
    axis.title = element_text(family = "Arial"),
    axis.text = element_text(family = "Arial"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

png("./Figures/1e-1.png", width = 400, height = 400)
print(p1)
dev.off()

diff_rds=readRDS("./Figures/rds/NA-vs-MO/umapdiff_rds.rds")
diff_rds_m<-diff_rds[[1]]
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

diff_rds_m$groups <- factor(diff_rds_m$groups,
                            levels = c("treatment", "control"),
                            labels = c("NA", "MO"))
umap_coords <- Embeddings(diff_rds_m, reduction = "umap")
x_range <- range(umap_coords[, 1])
y_range <- range(umap_coords[, 2])
x_expanded <- x_range + c(-0.05, 0.05) * diff(x_range)
y_expanded <- y_range + c(-0.05, 0.05) * diff(y_range)

plotm <- DimPlot(
  object = diff_rds_m,
  group.by = "groups",
  reduction = "umap",
  cols = cluster_Palette,
  pt.size = 1,
  ncol = 1,
  label = FALSE,
  label.size = 2
) +
  ggtitle("") +
  
  xlim(x_expanded[1], x_expanded[2]) +
  ylim(y_expanded[1], y_expanded[2]) +
  theme(legend.position = "right") +
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_text(family = "Arial", size = 12),
    legend.text = element_markdown(family = "Arial", size = 12),
    axis.title = element_text(family = "Arial", size = 12, face = "bold"),
    axis.text = element_text(family = "Arial"),
    axis.text.y = element_text(size = 10, angle = 0),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  )

png("./Figures/1e-2.png", width = 400, height = 400)
print(plotm)
dev.off()

volcano_data_m=readRDS("./Figures/rds/NA-vs-MO/volcano_data.rds")
group<-c("treatment","control")
FC_Threshold=2^0.26
pvalue<-0.05

#' @title Generate Volcano Plot
#' @description Creates a volcano plot to visualize differential analysis results (Log2FC vs -Log10 P-value).
#' Highlights Up-regulated, Down-regulated, and Non-significant features.
#' @param data Data frame containing differential results.
#' @param type Title/Label for the plot.
#' @param group Group labels.
#' @param pvalue Significance threshold.
#' @param FC_Threshold Fold change threshold.
#' @return A ggplot object.
volcano_plot_processing<-function(data,type,group,pvalue,FC_Threshold){
  Pvalue_type="p_val_adj"
  comparename=paste0(group[1],":",group[2])
  xlims <- ceiling(max(data$absfc))
  State_value <- unique(data$State)
  State_len <- length(State_value)
  ###color args
  if(State_len == 1){
    if(State_value == 'Non-significant'){
      scale_color <- "grey"
    }else if(State_value == "Down"){
      scale_color <- "LightSkyBlue"
    }else if(State_value == "Up"){
      scale_color <- "HotPink"
    }
  }else if(State_len == 2){
    if('Down' %in% State_value & 'Non-significant' %in% State_value){
      scale_color <- c("LightSkyBlue","grey")
    }else if('Non-significant' %in% State_value & 'Up' %in% State_value){
      scale_color <- c("grey","HotPink")
    }else if('Down' %in% State_value & 'Up' %in% State_value){
      scale_color <- c("LightSkyBlue","HotPink")
    }
  }else if(State_len == 3){
    scale_color <- c("LightSkyBlue","grey","HotPink")
  }
    p<- data %>%
      ggplot(aes(`log2(Fold Change)`,log_test))+
      theme_classic()+
      labs(title=paste(type))+
      geom_point(alpha= I(1/2),aes(color = State),size = 2.5)+
      scale_color_manual(values = scale_color)+
      scale_shape_manual(values = c(17, 16))+
      geom_hline(yintercept = -log10(pvalue),linetype=6,size = .3,color = "black")+
      geom_vline(xintercept=c(-log2(FC_Threshold),log2(FC_Threshold)),linetype=6,size = .3,color = "black")+
      xlim(-xlims,xlims)+
      xlab("log2(Fold Change)")+
      ylab(paste("-log10(",Pvalue_type,")",sep = ""))+
      theme(
        text = element_text(family = "Arial"),
        legend.title = element_text(family = "Arial", size = 12),
        legend.text = element_markdown(family = "Arial", size = 12),
        axis.title = element_text(family = "Arial", size = 12, face = "bold"),
        axis.text = element_text(family = "Arial"),
        axis.text.y = element_text(size = 10, angle = 0),
        plot.title = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        panel.border = element_blank()
      )
    
  
  return(p)
}
volcanoplot<-volcano_plot_processing(volcano_data_m,"Metabolite",group,pvalue,FC_Threshold)
png("./Figures/1e-3.png", width = 500, height = 400)
print(volcanoplot)
dev.off()

#1e-4---------------------------------------------------------------------------------------
#' @title Figure 1e-4: Case Group Network
#' @description Visualizes the correlation network for the Case/Treatment group.
#' Nodes colored by normalized mean expression, edges by correlation type.
#' @return Generates PNG '1e-4.png'.
cornetworkplot=readRDS("./Figures/rds/NA-vs-MO/cornetworkplot.rds")
nodes=cornetworkplot[[1]]
case_edges=cornetworkplot[[2]]
control_edges=cornetworkplot[[3]]

layout_type<-"kk"
color_type<-"Normalized mean"
show_node_name<-TRUE
node_name_size<-1.5
node_size<-4

edges=case_edges
case_igraph <- igraph::graph_from_data_frame(d=as.data.frame(edges),vertices=as.data.frame(nodes),directed = F)
igraph::E(case_igraph)$color <- ifelse(is.na(igraph::E(case_igraph)$cor), "other",
                                       ifelse(igraph::E(case_igraph)$cor > 0,"positive","negative"))
igraph::E(case_igraph)$edge_width<-ifelse(is.na(igraph::E(case_igraph)$cor), 0.5, igraph::E(case_igraph)$cor)
igraph::E(case_igraph)$edge_width<-abs(igraph::E(case_igraph)$edge_width)#不取正数的话负相关会透明
igraph::V(case_igraph)$size<-1
igraph::V(case_igraph)$`Normalized mean`<-igraph::V(case_igraph)$case_norm_mean

igraph=case_igraph
set.seed(123)
if(layout_type=="kk"){
  LAYOUT<-igraph::layout_with_kk(igraph)
}else if(layout_type=="nicely"){
  LAYOUT<-igraph::layout_nicely(igraph)
}else if(layout_type=="fr"){
  LAYOUT<-igraph::layout_with_fr(igraph)
}else{
  stop("The layout_type parameter must be either ‘kk’, ‘nicely’, or ‘fr’.")
}
colormap <-colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
edge_status <- c("positive", "negative", "other")
edgecolormapping <- setNames(c("red", "green", "gray"), edge_status)
edgecolormap <- edgecolormapping[unique(igraph::E(igraph)$color)]
xlim_max=max(LAYOUT[,1])+0.3
xlim_min=min(LAYOUT[,1])-0.3
ylim_max=max(LAYOUT[,2])+0.3
ylim_min=min(LAYOUT[,2])-0.3
nodeshape<-c(21,24)
omics_levels <- sort(unique(igraph::V(igraph)$Class))
shape_values <- setNames(nodeshape, omics_levels)
case_p <- ggraph::ggraph(igraph, layout = LAYOUT) + 
  ggraph::geom_edge_link(
    ggplot2::aes(
      edge_color = color, edge_width = edge_width),
    alpha = 1,
    show.legend = c(edge_color = TRUE, edge_width = FALSE, linetype = FALSE)
  ) + 
  ggraph::scale_edge_color_manual(
    values = edgecolormap, 
    name = "Correlation",
    guide = ggplot2::guide_legend(order = 2)
  ) + 
  ggraph::geom_node_point(
    ggplot2::aes(
      x = x, y = y, size = size,
      fill = !!rlang::sym(color_type), shape = Class
    ),
    colour = "black",
    show.legend = c(size = FALSE, fill = TRUE, shape = TRUE)
  ) +
  ggplot2::scale_shape_manual(
    name = "Omics",
    values = shape_values
  ) +
  ggplot2::scale_fill_gradientn(
    colours = colormap(100),
    limits = c(-1, 1),
    name = color_type,
    guide = ggplot2::guide_colorbar(order = 3)
  ) +
  ggplot2::scale_size(range = c(as.numeric(node_size), as.numeric(node_size))) +
  ggraph::scale_edge_width(range = c(0.2, 1), limits = c(0, 1)) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, family = "Arial"),
    text = element_text(family = "Arial"),
    legend.title = element_text(family = "Arial"),
    legend.text = element_text(family = "Arial")
  ) +
  coord_equal() + 
  xlim(xlim_min, xlim_max) + 
  ylim(ylim_min, ylim_max) + 
  ggraph::geom_node_text(
    ggplot2::aes(label = sapply(strsplit(name, ";"), function(x) x[1])),  
    vjust = -0.5,
    size = node_name_size,
    family = "Arial"
  )

png("./Figures/1e-4.png", width = 400, height = 400)
print(case_p)
dev.off()

#1e-5---------------------------------------------------------------------------------------
#' @title Figure 1e-5: Control Group Network
#' @description Visualizes the correlation network for the Control group for comparison.
#' @return Generates PNG '1e-5.png'.
edges=control_edges
control_igraph <- igraph::graph_from_data_frame(d=as.data.frame(edges),vertices=as.data.frame(nodes),directed = F)
igraph::E(control_igraph)$color <- ifelse(is.na(igraph::E(control_igraph)$cor), "other",
                                          ifelse(igraph::E(control_igraph)$cor > 0,"positive","negative"))
igraph::E(control_igraph)$edge_width<-ifelse(is.na(igraph::E(control_igraph)$cor), 0.5, igraph::E(control_igraph)$cor)
igraph::E(control_igraph)$edge_width<-abs(igraph::E(control_igraph)$edge_width)#不取正数的话负相关会透明
igraph::V(control_igraph)$size<-1
igraph::V(control_igraph)$`Normalized mean`<-igraph::V(control_igraph)$control_norm_mean
igraph=control_igraph

set.seed(123)
if(layout_type=="kk"){
  LAYOUT<-igraph::layout_with_kk(igraph)
}else if(layout_type=="nicely"){
  LAYOUT<-igraph::layout_nicely(igraph)
}else if(layout_type=="fr"){
  LAYOUT<-igraph::layout_with_fr(igraph)
}else{
  stop("The layout_type parameter must be either ‘kk’, ‘nicely’, or ‘fr’.")
}
colormap <-colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
edge_status <- c("positive", "negative", "other")
edgecolormapping <- setNames(c("red", "green", "gray"), edge_status)
edgecolormap <- edgecolormapping[unique(igraph::E(igraph)$color)]
xlim_max=max(LAYOUT[,1])+0.3
xlim_min=min(LAYOUT[,1])-0.3
ylim_max=max(LAYOUT[,2])+0.3
ylim_min=min(LAYOUT[,2])-0.3
nodeshape<-c(21,24)
omics_levels <- sort(unique(igraph::V(igraph)$Class))
shape_values <- setNames(nodeshape, omics_levels)
con_p <- ggraph::ggraph(igraph, layout = LAYOUT) + 
  ggraph::geom_edge_link(
    ggplot2::aes(
      edge_color = color, edge_width = edge_width),
    alpha = 1,
    show.legend = c(edge_color = TRUE, edge_width = FALSE, linetype = FALSE)
  ) + 
  ggraph::scale_edge_color_manual(
    values = edgecolormap, 
    name = "Correlation",
    guide = ggplot2::guide_legend(order = 2)
  ) + 
  ggraph::geom_node_point(
    ggplot2::aes(
      x = x, y = y, size = size,
      fill = !!rlang::sym(color_type), shape = Class
    ),
    colour = "black",
    show.legend = c(size = FALSE, fill = TRUE, shape = TRUE)
  ) +
  ggplot2::scale_shape_manual(
    name = "Omics",
    values = shape_values
  ) +
  ggplot2::scale_fill_gradientn(
    colours = colormap(100),
    limits = c(-1, 1),
    name = color_type,
    guide = ggplot2::guide_colorbar(order = 3)
  ) +
  ggplot2::scale_size(range = c(as.numeric(node_size), as.numeric(node_size))) +
  ggraph::scale_edge_width(range = c(0.2, 1), limits = c(0, 1)) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, family = "Arial"),
    text = element_text(family = "Arial"),
    legend.title = element_text(family = "Arial"),
    legend.text = element_text(family = "Arial")
  ) +
  coord_equal() + 
  xlim(xlim_min, xlim_max) + 
  ylim(ylim_min, ylim_max) + 
  ggraph::geom_node_text(
    ggplot2::aes(label = sapply(strsplit(name, ";"), function(x) x[1])),  
    vjust = -0.5,
    size = node_name_size,
    family = "Arial"
  )
png("./Figures/1e-5.png", width = 400, height = 400)
print(con_p)
dev.off()
########################################1f
#' @title Figure 1f: Functional Analysis
#' @description Composite figure for Functional Association Analysis.
#' 1f-2: Venn diagram of pathway overlap between Metabolites and Genes.
#' 1f-3: Bubble plot of annotated pathways.
#' 1f-1: KEGG pathway map with mapped differential features.
#' @return Generates PNGs '1f-2.png', '1f-3.png', '1f-1.png'.
vennplot_data=readRDS("./Figures/rds/CA-vs-PAG/vennplot_datanew.RDS")
p3 <- venn.diagram(vennplot_data,
                   col = "transparent",
                   fill = c('yellow', 'skyblue'),
                   alpha = 0.5, cex = 2,
                   fontfamily = "Arial",
                   cat.cex = 2, cat.pos = 0,
                   cat.fontfamily = "Arial",
                   rotation.degree = 0,
                   resolution = 300, filename = NULL)

Cairo("./Figures/1f-2.png", width = 400, height = 400, family = "Arial",bg = "white")
grid.draw(p3)
dev.off()

df=readRDS("./Figures/rds/CA-vs-PAG/dotplot_save.RDS")
plot_title = ""
plot_title_size = 15
text_size =10 
axis_title_size = 12
legend_title_size = 12
legend_text_size =12 
top_num=20
# Calculate Total Count for each Pathway to determine Top 20
pathway_counts <- df %>%
  dplyr::group_by(Pathway) %>%
  dplyr::summarise(TotalCount = sum(Count)) %>%
  dplyr::arrange(desc(TotalCount)) %>%
  dplyr::slice_head(n = top_num)

top_pathways <- pathway_counts$Pathway

# Filter original df to keep only Top 20 pathways
df <- df %>% 
  dplyr::filter(Pathway %in% top_pathways) %>%
  # Join with total counts for sorting
  dplyr::left_join(pathway_counts, by = "Pathway") %>%
  dplyr::mutate(Pathway = factor(Pathway, levels = rev(top_pathways))) # Sort Y-axis by Total Count
# Plot: Annotation Count Dot Plot
# Color and Shape both map to Types (Gene/Metabolite)
g <- ggplot(df)+
  geom_point(aes(x= Count, y= Pathway, colour= Types, shape = Types), size = 4)+ 
  theme_minimal() +
  labs(title="", x="Annotation Count", y="Pathway", colour="Omics Type", shape="Omics Type")+
  scale_shape_manual(values = c("Gene" = 16, "Metabolite" = 17, "Co-annotation" = 15)) + # Circle, Triangle
  scale_x_continuous(breaks = function(x) unique(floor(pretty(x)))) +
  scale_colour_manual(values = c("Gene" = "#377EB8", "Metabolite" = "#E41A1C", "Co-annotation" = "#4DAF4A")) + # Blue, Red, Green
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 10, angle = 0),
    plot.title = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", size = 1), 
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  )


# Title styling
if(plot_title == ""){
  g <- g + theme(plot.title = element_blank())
}else{
  g <- g + theme(plot.title = element_text(size = plot_title_size, hjust = 0.5, face="bold"))
}

# Axis title styling
if(axis_title_size == 0){
  g <- g + theme(axis.title = element_blank())
}else{
  g <- g + theme(axis.title = element_text(size = axis_title_size, face="bold"))
}
png("./Figures/1f-3.png", width = 600, height = 400, family = "Arial")
print(g)
dev.off()

annotation_plotdata=readRDS("./Figures/rds/CA-vs-PAG/PAG/Retrograde endocannabinoid signaling/annotation_plotdata.RDS")
m_fplot <- annotation_plotdata[[1]]
t_fplot <- annotation_plotdata[[2]]
pathway_id <- annotation_plotdata[[5]]
species <- "mmu"
devs <- dev.list()
if (!is.null(devs)) {
  for (d in devs) {
    tryCatch({
      dev.off(d)
    }, error = function(e) {})
  }
}

mappath <- "./source/Database/KEGG/map"
img <- image_read(file.path(mappath, paste0(pathway_id, ".png")))
img <- image_convert(img, colorspace = 'gray')

lines <- suppressWarnings(readLines(file.path(mappath, paste0(pathway_id, ".conf"))))
dim <- image_info(img)
width <- dim$width
height <- dim$height

canvas <- image_blank(width, height, "none")

img_draw <- image_draw(canvas)
on.exit(dev.off(), add = TRUE)

circle <- function(x, y, r, border = "black", lty = "solid", lwd = 1) {
  angle <- seq(0, 2 * pi, length.out = 100)
  x_circle <- x + r * cos(angle)
  y_circle <- y + r * sin(angle)
  lines(x_circle, y_circle, col = border, lty = lty, lwd = lwd)
}

for (line in lines) {
  line <- trimws(line)
  
  if (grepl("^circ", line) && any(!is.na(m_fplot))) {
    coords <- as.numeric(str_extract_all(line, "\\d+")[[1]])
    id_match <- str_match(line, "(C\\d{5})")
    matched_ids <- intersect(id_match, names(m_fplot))
    if (length(matched_ids) > 0 && length(coords) >= 3) {
      x <- coords[1]
      y <- coords[2]
      r <- coords[3]
      
      values <- unlist(m_fplot[matched_ids])
      median_value <- median(values)
      
      if (!is.na(median_value)) {
        color <- ifelse(median_value > 0, "orange",
                        ifelse(median_value < 0, "blue", "gray"))
        circle(x, y, r, border = color, lty = "solid", lwd = 3)
      }
    }
  } else if (grepl("^rect", line) && any(!is.na(t_fplot))) {      
    coords <- as.numeric(str_extract_all(line, "\\d+")[[1]])
    ids <- str_extract_all(line, "K\\d{5}")[[1]]
    matched_ids <- intersect(ids, names(t_fplot))
    
    if (length(matched_ids) > 0 && length(coords) >= 4) {
      xleft <- coords[1]
      ybottom <- coords[2]
      xright <- coords[3]
      ytop <- coords[4]
      
      values <- unlist(t_fplot[matched_ids])
      median_value <- median(values)
      if (!is.na(median_value)) {
        color <- ifelse(median_value > 0, "red",
                        ifelse(median_value < 0, "green", "gray"))
        rect(xleft, ybottom, xright, ytop, border = color, lty = "solid", lwd = 3)
      }
    }
  }
}


img <- image_composite(img, img_draw, offset = "+0+0")

output_file <-  file.path("./Figures",  "1f-1.png")
target_width <- 4  
target_height <- 4 
dpi <- 300  
width_px <- target_width * dpi
height_px <- target_height * dpi

img_adjusted <- img %>%
  image_scale(geometry = sprintf("%dx%d", width_px, height_px)) %>% 
  image_extent(geometry = sprintf("%dx%d", width_px, height_px),   
               gravity = "center",                                  
               color = "white")                                  

image_write(
  img_adjusted,
  path = output_file,
  format = "png",
  density = dpi  
)
devs <- dev.list()
if (!is.null(devs)) {
  for (d in devs) {
    tryCatch({
      dev.off(d)
    }, error = function(e) {})
  }
}
########################################1g
#' @title Figure 1g: Multi-Modal Visualization
#' @description Composite figure for Data Visualization module.
#' 1g-1: Single ion spatial heatmap.
#' 1g-2: Grid of top positively correlated ions.
#' 1g-3: RGB co-visualization of selected features.
#' @return Generates PNGs '1g-1.png', '1g-2.png', '1g-3.png'.
single_ionsplot_save=readRDS("./Figures/rds/single_ionsplot_save.rds")
data=single_ionsplot_save[[1]]
title_text=single_ionsplot_save[[2]]
heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
p <- ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = norm_intensity), size = 1.5) +
  scale_color_gradientn(name = "Relative Abundance",colours = heatmap_Palette(100)) +
  ggtitle(title_text)+
  coord_equal() +
  xlim(range(data$x)) +
  ylim(range(data$y)) +
  theme_minimal() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 10),
    panel.grid = element_blank()
  )
png("./Figures/1g-1.png", width = 400, height = 400, family = "Arial")
print(p)
dev.off()

spatial_coordlist=readRDS("./Figures/rds/positively_correlated_ions_m.rds")
correlated_ions=spatial_coordlist[[1]]
spatial_coord=spatial_coordlist[[2]]
k=tail(correlated_ions)
plot_list <- list()
for(i in c(1:nrow(k))){
  featurename<- k$feature[i]
  w <- which(names(spatial_coord) %in% c("x","y",featurename))
  data <- spatial_coord %>%
    dplyr::select(all_of(w))
  
  names(data)[3] <- "intensity"
  data$norm_intensity <-100*(data$intensity)/max(data$intensity)
  
  title_text<-unlist(strsplit(as.character(featurename), ";"))[1]
  if(k$p[i]<0.001){
    ptext="p<0.001"
  }else{
    ptext=paste0("p:",k$p[i])
  }
  sub_title<-paste("Cor:",round(k$cor[i],2),"; ",ptext,sep = "")
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  plot_list[[i]] <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = norm_intensity), size = 1) +
    scale_color_gradientn(name = "Relative Abundance",colours = heatmap_Palette(100)) +
    coord_equal() +
    xlim(range(data$x)) +
    ylim(range(data$y)) +
    theme_minimal() +
    ggtitle(title_text,subtitle =sub_title )+
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5,size = 12),
      plot.subtitle = element_text(hjust = 0.5,size = 10),
      panel.grid = element_blank()
    )
  #eval(parse(text = paste("p",i," <- p",sep = "")))
}
png("./Figures/1g-2.png", width = 800, height = 600, family = "Arial")
do.call(gridExtra::grid.arrange, c(plot_list, ncol = 3))
dev.off()

k <- tail(correlated_ions)
spatial_coordlist=readRDS("./Figures/rds/CA-vs-PAG/PAG/Retrograde endocannabinoid signaling/cor_gene_metabolite_Co_visualisation_plot_save.rds")
spatial_coord=spatial_coordlist[[1]]
pair=spatial_coordlist[[2]]
pt.size=0.5
pairname=c("GABA",NA,"Slc32a1")
LRpair =pair[!is.na(pair)]
data = spatial_coord[,c('x','y',LRpair)] 
if(!is.na(pair[1])){
  f1=data[,pair[1]]
}else{
  f1=0
}
if(!is.na(pair[2])){
  f2=data[,pair[2]]
}else{
  f2=0
}
if(!is.na(pair[3])){
  f3=data[,pair[3]]
}else{
  f3=0
}


data$color <- rgb(f1,f2,f3, maxColorValue = 1)


p5<-ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = color), size = pt.size) +
  xlim(min(data$x), max(data$x)*1.3) +
  ylim(min(data$y), max(data$y)) +
  coord_equal() +
  scale_color_identity() +  
  xlab("")+ylab("")+
  theme_minimal() +
  theme(legend.text = element_text(size = 4 ) ,
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) 
if(!is.na(pair[1])){
  
  p5<-p5 + annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)*1.4/2, ymax = max(data$y)*1.4/2+4, 
                    fill = "red", color = "black") +  
    annotate("text", x = max(data$x)+ 5, y = max(data$y)*1.4/2+2, label = pairname[1], 
             hjust = 0, vjust = 0.5, size = 3, color = "black", family = "Arial")   
  
}
if(!is.na(pair[2])){
  
  p5<-p5 +annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)/2, ymax = (max(data$y)/2)+4, 
                   fill = "green", color = "black") +  
    annotate("text", x = max(data$x)+ 5, y = (max(data$y)/2)+2, label = pairname[2], 
             hjust = 0, vjust = 0.5, size = 3, color = "black", family = "Arial")   
}
if(!is.na(pair[3])){
  
  p5<-p5 +annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)*1.2/2, ymax = (max(data$y)*1.2/2)+4, 
                   fill = "blue", color = "black") +  
    annotate("text", x = max(data$x)+ 5, y =  (max(data$y)*1.2/2)+2, label =pairname[3], 
             hjust = 0, vjust = 0.5, size = 3, color = "black", family = "Arial", fontface = "italic")  
}
png("./Figures/1g-3.png", width = 400, height = 400, family = "Arial")
print(p5)
dev.off()
