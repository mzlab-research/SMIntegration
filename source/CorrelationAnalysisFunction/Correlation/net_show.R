net_show<-function(igraph,layout_type="fr",color_type ="Normalized mean",show_node_name=FALSE,
                   node_name_size=2,node_size=6){
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
  if(color_type=="Class"){
  colormap <- setNames(c("#348BBB","#F47346"), c("metabolite","gene"))
  }else if(color_type=="Normalized mean"){
  colormap <-colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
  }
  edge_status <- c("positive", "negative", "other")
  edgecolormapping <- setNames(c("red", "green", "gray"), edge_status)
  edgecolormap <- edgecolormapping[unique(igraph::E(igraph)$color)]
  xlim_max=max(LAYOUT[,1])+0.3
  xlim_min=min(LAYOUT[,1])-0.3
  ylim_max=max(LAYOUT[,2])+0.3
  ylim_min=min(LAYOUT[,2])-0.3
  plot <- run_ggraph_plot(igraph=igraph,colormap=colormap,edgecolormap=edgecolormap, color_type =color_type,
                          plot_layout=LAYOUT,show_node_legend=TRUE,show_edge_legend=TRUE,alpha=1,
                          edge_color_type="color",node_size=node_size) 
  plot <-plot + xlim(xlim_min, xlim_max)+ 
    ylim(ylim_min, ylim_max)
  if (show_node_name){
    plot <- plot + ggraph::geom_node_text(ggplot2::aes(label =name), vjust = -0.5,size = node_name_size)
  }
  return(plot)
} 


run_color<-function(annotation_table=NULL){
  unique_classes<-unique(annotation_table$Class)
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  if(length(unique_classes)>length(cluster_Palette)){
    cluster_Palette <-colorRampPalette(cluster_Palette)(length(unique_classes))#避免颜色不够
  }
  color_mapping <- setNames(cluster_Palette[1:length(unique_classes)], unique_classes)
  return(color_mapping)
}
run_igraph <- function(nodes=NULL,edges=NULL){
  net1 <- igraph::graph_from_data_frame(d=as.data.frame(edges),vertices=as.data.frame(nodes),directed = F)

  igraph::E(net1)$color <- ifelse(is.na(igraph::E(net1)$cor), "other",
                                  ifelse(igraph::E(net1)$cor > 0,"positive","negative"))
  igraph::E(net1)$edge_width<-ifelse(is.na(igraph::E(net1)$cor), 0.5, igraph::E(net1)$cor)
  igraph::E(net1)$edge_width<-abs(igraph::E(net1)$edge_width)#不取正数的话负相关会透明
  igraph::V(net1)$size<-1
  
  igraph::V(net1)$betweenness <- igraph::betweenness(
    graph = net1,
    v = igraph::V(net1),
    directed = FALSE
  )
  igraph::V(net1)$degree <- igraph::degree(
    graph = net1,
    v= igraph::V(net1)
  )
  igraph::V(net1)$eigenvector <- igraph::eigen_centrality(
    graph = net1
    , directed = FALSE)$vector
  return(net1)
}

run_ggraph_plot <- function(igraph=NULL, colormap=NULL,edgecolormap=NULL, color_type = "Normalized mean",plot_layout="fr",
                            edge_color_type="color",
                            alpha = 1,show_edge_legend=FALSE,show_node_legend=FALSE,node_size=6) {
  plot <- ggraph::ggraph(igraph, layout = plot_layout)

    plot <- plot + ggraph::geom_edge_link(
      ggplot2::aes(
        edge_color =!!rlang::sym(edge_color_type), edge_width = edge_width),
      alpha = alpha,
      show.legend = c(edge_color = show_edge_legend, edge_width = FALSE, linetype = FALSE)#, size = FALSE

    )

  if(edge_color_type=="color"){
    plot <- plot + ggraph::scale_edge_color_manual(values = edgecolormap, name = "Correlation",
                                                   guide = ggplot2::guide_legend(order=2))
  }else{
    plot <- plot + ggraph::scale_edge_color_manual(values = edgecolormap, name = "Correlation status",
                                                   guide = ggplot2::guide_legend(order=2))
  }

    nodeshape<-c(21,24)

    omics_levels <- sort(unique(igraph::V(igraph)$Class))

    shape_values <- setNames(nodeshape, omics_levels)
    plot <- plot + ggraph::geom_node_point(
      ggplot2::aes(x=x,y=y,size = size,
                   fill = !!rlang::sym(color_type),shape=Class),#color
             
      colour = "black",
      show.legend = c(size = FALSE, fill = show_node_legend,shape=show_node_legend)#color
    ) +
      ggplot2::scale_shape_manual(name="Omics",values =shape_values) 
 
  if(is.function(colormap)){
    if(color_type=="Normalized mean"){
      plot <- plot + ggplot2::scale_fill_gradientn(colours = colormap(100),limits = c(-1, 1),
                                                   name = color_type,
                                                   guide = ggplot2::guide_colorbar(order=3)) 
    }else{
      plot <- plot + ggplot2::scale_fill_gradientn(colours = colormap(100), name = paste0(color_type,"(case:control)"),
                                                   guide = ggplot2::guide_colorbar(order=3)) 
    }
  }else{
    plot <- plot + ggplot2::scale_fill_manual(values = colormap, name = color_type,
                                              guide = ggplot2::guide_legend(order=3)) #scale_color_manual
  }

  plot <- plot + ggplot2:: scale_size(range = c(as.numeric(node_size), as.numeric(node_size)))


  plot <- plot + ggraph::scale_edge_width(range = c(0.2, 1),limits = c(0, 1))
  

  plot <- plot +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    coord_equal()
  
  return(plot)
}

run_compare_edge<-function(edges){
  edges <-transform(edges, sorted_from = pmin(from, to), sorted_to = pmax(from, to)) 
  edges$from_to<-paste0(edges$sorted_to,"_",edges$sorted_from)
  edges<-edges |>
    dplyr::select(-sorted_from, -sorted_to) 
  return(edges)
}