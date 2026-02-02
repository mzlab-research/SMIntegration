AnnotationDotPlot <- function(
    df,
    plot_title = "Top 20 Annotated Pathways", 
    plot_title_size = 15, 
    font_family='serif', 
    text_size =10 , 
    axis_title_size = 12, 
    legend_title_size = 12, 
    legend_text_size =12 
){
  
  # Calculate Total Count for each Pathway to determine Top 20
  pathway_counts <- df %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarise(TotalCount = sum(Count)) %>%
    dplyr::arrange(desc(TotalCount)) %>%
    dplyr::slice_head(n = 20)
  
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
    theme_bw()+
    labs(title=plot_title, x="Annotation Count", y="Pathway", colour="Omics Type", shape="Omics Type")+
    scale_shape_manual(values = c("Gene" = 16, "Metabolite" = 17, "Co-annotation" = 15)) + # Circle, Triangle, Square
    scale_colour_manual(values = c("Gene" = "#377EB8", "Metabolite" = "#E41A1C", "Co-annotation" = "#4DAF4A")) + # Blue, Red, Green
    scale_x_continuous(breaks = function(x) unique(floor(pretty(x)))) + # Ensure integer breaks on x-axis
    theme(legend.title=element_text(size=legend_title_size,family = font_family),
          legend.text=element_text(size=legend_text_size,family = font_family),
          axis.text.y = element_text(size=text_size,angle=0,family = font_family),
          panel.grid.major.x = element_line(linetype = "dashed"), # Vertical grid lines (dashed)
          panel.grid.major.y = element_line(linetype = "dotted")) # Horizontal grid lines (dotted)
  
  
  # Title styling: hide title if plot_title is empty string
  if(plot_title == ""){
    g <- g + theme(plot.title = element_blank())
  }else{
    g <- g + theme(plot.title = element_text(size = plot_title_size, hjust = 0.5, face="bold", family = font_family))
  }
  
  # Axis title styling: hide axis titles if axis_title_size is 0
  if(axis_title_size == 0){
    g <- g + theme(axis.title = element_blank())
  }else{
    g <- g + theme(axis.title = element_text(size = axis_title_size, face="bold", family = font_family))
  }
  
  return(g)
}