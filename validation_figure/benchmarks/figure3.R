suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(showtext))
suppressMessages(library(Seurat))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggtext))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
# font_add("sans", regular = "arial.ttf", italic = "ariali.ttf")
# showtext_auto()
#3a---------------------------------------------------------------------------------------
#' @title Figure 3a: Spatial Group Assignment
#' @description Visualizes the spatial distribution of Treatment (NA) and Control (MO) groups.
#' @return Generates a PDF '3a.pdf'.
data=readRDS("./Figures/rds/NA-vs-MO/plotdata_group.rds")
data <- data |>
  group_by(groups) 
group_color<-c("treatment"="red","control"="blue","0"="grey")
p1 <-  ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = groups), size = 1) +
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
    text = element_text(family = "sans"),
    legend.title = element_text(family = "sans"),
    legend.text = element_markdown(family = "sans"),
    axis.title = element_text(family = "sans"),
    axis.text = element_text(family = "sans"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

pdf("./Figures/3a.pdf", width = 6, height = 6)
print(p1)
dev.off()

#3b+3c---------------------------------------------------------------------------------------
#' @title Figure 3b & 3c: UMAP Visualization
#' @description Generates UMAP plots for Metabolomics (3b) and Transcriptomics (3c)
#' comparing Treatment (NA) and Control (MO) groups.
#' @return Generates PDFs '3b.pdf' and '3c.pdf'.
diff_rds=readRDS("./Figures/rds/NA-vs-MO/umapdiff_rds.rds")
diff_rds_m<-diff_rds[[1]]
diff_rds_t<-diff_rds[[2]]
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
    text = element_text(family = "sans"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_markdown(family = "sans", size = 12),
    axis.title = element_text(family = "sans", size = 12, face = "bold"),
    axis.text = element_text(family = "sans"),
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

pdf("./Figures/3b.pdf", width = 6, height = 6)
print(plotm)
dev.off()


diff_rds_t$groups <- factor(diff_rds_t$groups,
                            levels = c("treatment", "control"),
                            labels = c("NA", "MO"))
umap_coords <- Embeddings(diff_rds_t, reduction = "umap")
x_range <- range(umap_coords[, 1])
y_range <- range(umap_coords[, 2])

x_expanded <- x_range + c(-0.05, 0.05) * diff(x_range)
y_expanded <- y_range + c(-0.05, 0.05) * diff(y_range)
x_expanded<-c(-15, 10)
y_expanded<-c(-30,10)
plott <- DimPlot(
  object = diff_rds_t,
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
    text = element_text(family = "sans"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_markdown(family = "sans", size = 12),
    axis.title = element_text(family = "sans", size = 12, face = "bold"),
    axis.text = element_text(family = "sans"),
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

pdf("./Figures/3c.pdf", width = 6, height = 6)
print(plott)
dev.off()

#3d---------------------------------------------------------------------------------------
#' @title Figure 3d: Differential Feature Counts
#' @description Bar plot showing the number of Up- and Down-regulated features for each sample type.
#' @return Generates a PDF '3d.pdf'.
data=readRDS("./Figures/rds/NA-vs-MO/data_count_long_all.rds")
UD_Class <- unique(data$State)
UD_len <- length(UD_Class)

if(UD_len == 1){
  if(UD_Class == 'Up'){
    scale_fill_values <- "HotPink"
    scale_fill_labels <- "Up"
  }else if(UD_Class == 'Down'){
    scale_fill_values <- "LightSkyBlue"
    scale_fill_labels <- "Down"
  }
}else{
  scale_fill_values <- c("HotPink", "LightSkyBlue")
  scale_fill_labels <- c("Up","Down")
}
p2 <- ggplot(data,mapping = aes(x=Sample_type, y=value, fill=State)) +
  theme_bw()+ theme(panel.grid=element_blank())+
  geom_bar(stat="identity",width = 0.5, position = position_dodge(0.7)) +
  geom_text(aes(label=value),hjust = 0.3,vjust = 1,position = position_dodge(0.7))+
  scale_fill_manual(values = scale_fill_values, labels = scale_fill_labels)+
  labs(x="",y="")+
  scale_y_continuous(
    limits = c(0, 1500),
    breaks = seq(0, 1500, by = 500)  
  ) +
  theme(
    text = element_text(family = "sans"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_markdown(family = "sans", size = 12),
    axis.title = element_text(family = "sans", size = 12, face = "bold"),
    axis.text = element_text(family = "sans"),
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

pdf("./Figures/3d.pdf", width = 6, height = 6)
print(p2)
dev.off()
#3e+3f---------------------------------------------------------------------------------------
#' @title Figure 3e & 3f: Differential Correlation Networks
#' @description Constructs and visualizes correlation networks for Case (3e) and Control (3f) groups.
#' Nodes are colored by normalized mean expression, edges by correlation type.
#' @return Generates PDFs '3e.pdf' and '3f.pdf'.
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
    plot.title = ggplot2::element_text(hjust = 0.5, family = "sans"),
    text = element_text(family = "sans"),
    legend.title = element_text(family = "sans"),
    legend.text = element_text(family = "sans")
  ) +
  coord_equal() + 
  xlim(xlim_min, xlim_max) + 
  ylim(ylim_min, ylim_max) + 
  ggraph::geom_node_text(
    ggplot2::aes(label = sapply(strsplit(name, ";"), function(x) x[1])),  
    vjust = -0.5,
    size = node_name_size,
    family = "sans"
  )
pdf("./Figures/3e.pdf", width = 6, height = 6)
print(case_p)
dev.off()


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
    plot.title = ggplot2::element_text(hjust = 0.5, family = "sans"),
    text = element_text(family = "sans"),
    legend.title = element_text(family = "sans"),
    legend.text = element_text(family = "sans")
  ) +
  coord_equal() + 
  xlim(xlim_min, xlim_max) + 
  ylim(ylim_min, ylim_max) + 
  ggraph::geom_node_text(
    ggplot2::aes(label = sapply(strsplit(name, ";"), function(x) x[1])),  
    vjust = -0.5,
    size = node_name_size,
    family = "sans"
  )


pdf("./Figures/3f.pdf", width = 6, height = 6)
print(con_p)
dev.off()
