

BubbleDiagram <- function(
    df,
    plot_title = "", 
    plot_title_size = 15, 
    font_family='serif', 
    text_size =10 , 
    axis_title_size = 12, 
    legend_title_size = 12, 
    legend_text_size =12 
){
 
    maxnum <- df %>% dplyr::select(all_of("Count")) %>%
      max() %>%
      unique()
    
    g <- ggplot(df)+
      geom_point(aes(x= RichFactor, y= Pathway, size= Count, colour= Pvalue, shape = Types))+
      theme_bw()+
      labs(title=plot_title)+
      theme(legend.title=element_text(size=legend_title_size,family = font_family),
            legend.text=element_text(size=legend_text_size,family = font_family),
            axis.text.y = element_text(size=text_size,angle=0,family = font_family))
    if(maxnum < 5){
      g <- g + scale_size_continuous(limits=c(1,5),range=c(1,5))
    }else{
      g <- g + scale_size_continuous(range=c(2,8))
    }

    if(sum(grepl("^Pvalue$",names(df))) > 0){
      pvaluenum <- df %>% dplyr::select(all_of("Pvalue")) %>%
        max() %>%
        unique()
      g <- g + scale_colour_continuous(low="red",high="blue",limits=c(0,pvaluenum))
    }else{
      classnum <- df %>% dplyr::select(all_of("Pvalue")) %>%
        unique() %>%
        nrow()
      if(classnum >20){
        g <- g + guides(colour = guide_legend(ncol = 3))
      }else{
        g <- g + guides(colour = guide_legend(ncol = 2))
      }
    }

    if(plot_title == ""){
      g <- g + theme(plot.title = element_blank())
    }else{
      g <- g + theme(plot.title = element_text(size = plot_title_size, hjust = 0.5, face="bold", family = font_family))
    }

    if(axis_title_size == 0){
      g <- g + theme(axis.title = element_blank())
    }else{
      g <- g + theme(axis.title = element_text(size = axis_title_size, face="bold", family = font_family))
    }
    
return(g)
  
}
