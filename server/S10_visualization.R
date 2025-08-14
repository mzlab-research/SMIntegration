
output$download_single_feature_data <- downloadHandler(
  filename = function() {

    paste0("single_feature_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData=single_ionsplotdata()
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
      
    })
  })
output$download_single_feature_plot <- downloadHandler(
  filename = function() {
    "single_feature_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      p <- single_ionsplot_save()

      ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
      
    })
  })



output$download_positively_correlated_ions_m_data <- downloadHandler(
  filename = function() {

    paste0("positively_correlated_metabolites_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData=positively_correlated_ions_m()[[2]]
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
      
    })
  })
output$download_positively_correlated_ions_m_plot <- downloadHandler(
  filename = function() {
    "positively_correlated_metabolites_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      pp<- positively_correlated_ions_m()[[1]]
      p<-gridExtra::grid.arrange(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],ncol=3)

      ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
      
    })
  })

output$download_negatively_correlated_ions_m_data <- downloadHandler(
  filename = function() {

    paste0("negatively_correlated_metabolites_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData=negatively_correlated_ions_m()[[2]]
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
      
    })
  })
output$download_negatively_correlated_ions_m_plot <- downloadHandler(
  filename = function() {
    "negatively_correlated_metabolites_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      pp <- negatively_correlated_ions_m()[[1]]
      p<-gridExtra::grid.arrange(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],ncol=3)

      ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
      
    })
  })

output$download_positively_correlated_ions_t_data <- downloadHandler(
  filename = function() {

    paste0("positively_correlated_genes_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData=positively_correlated_ions_t()[[2]]
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
      
    })
  })
output$download_positively_correlated_ions_t_plot <- downloadHandler(
  filename = function() {
    "positively_correlated_genes_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      pp <- positively_correlated_ions_t()[[1]]
      p<-gridExtra::grid.arrange(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],ncol=3)

      ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
      
    })
  })

output$download_negatively_correlated_ions_t_data <- downloadHandler(
  filename = function() {

    paste0("negatively_correlated_genes_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData=negatively_correlated_ions_t()[[2]]
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
      
    })
  })
output$download_negatively_correlated_ions_t_plot <- downloadHandler(
  filename = function() {
    "negatively_correlated_genes_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      pp <- negatively_correlated_ions_t()[[1]]
      p<-gridExtra::grid.arrange(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],ncol=3)

      ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
      
    })
  })


spatial_coord<-eventReactive(c(input$start_visualisation_analysis), {
  req(data_rds())
  withProgress(message = "Processing data...",value=0.8,{
    data_rds<-data_rds()
    data_mrds<-data_rds[[1]]
    data_trds<-data_rds[[2]]
    decon_mtrx = data_mrds@assays$Spatial$counts
    decon_ttrx = data_trds@assays$Spatial$counts
    stopifnot(identical(colnames(decon_mtrx), colnames(decon_ttrx)))
    combined_matrix <- rbind2(decon_mtrx, decon_ttrx)
    meta.data<-data_mrds@meta.data
    source("./source/CorrelationAnalysisFunction/Co_visualisation/Co_visualisation.R")
    spatial_coord<-spatial_coord_processing(combined_matrix=combined_matrix,
                                            meta.data=meta.data,rescale=FALSE)
    gc()
    return(spatial_coord)
  })
})


observe({
  data_rds<-data_rds()
  req(data_rds)
  data_mrds<-data_rds[[1]]
  data_trds<-data_rds[[2]]

  if(input$single_ions_type=="Gene"){
    g = unique(as.character(rownames(data_trds@assays$Spatial$counts)))
    g  <- sort(g)
  }else{
    g = unique(as.character(rownames(data_mrds@assays$Spatial$counts)))
    g  <- sort(g)
  }
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      inputId = "select_single_ions",
      choices = g,
      server = TRUE,
      selected = if (length(g) > 0) g[1] else NULL,
      options = list(
        maxOptions = 1000,          
        searchConjunction = 'and',   
        render = I("{                
      option: function(item, escape) {
        return '<div>' + escape(item.label) + '</div>';
      }
    }"),
        loadThrottle = 300,          
        placeholder = "Search by feature name...",
        closeAfterSelect = TRUE      
      )
    )

})


single_ionsplotdata <-eventReactive(c(input$start_visualisation_analysis), {
  req(input$select_single_ions)
  withProgress(message = "Processing data...",value=0.8,{
    ion <- tryCatch(
      as.character(input$select_single_ions),
      error = function(e) {
        showNotification("Invalid feature selection!", type = "error")
        return(NULL)
      })
    req(ion, nzchar(ion)) 

    spatial_coord<-spatial_coord()
    plotdata<-spatial_coord |>
      dplyr::select(x,y,any_of(ion))
    colnames(plotdata)[3]<-"intensity"
    plotdata$norm_intensity <-100*(plotdata$intensity)/max(plotdata$intensity)
    req(nrow(plotdata) > 0) 
    return(plotdata)
  })
})
single_ionsplot_save <- reactive({ 
  data <- single_ionsplotdata()
  req(input$select_single_ions)
  req(data)

  title_text<-unlist(strsplit(as.character(input$select_single_ions), ";"))[1]

  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = norm_intensity), size = 1) +
    scale_color_gradientn(colours = heatmap_Palette(100)) +
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
  return(p)
  
})
observeEvent(input$start_visualisation_analysis,{
output$single_ion_visualizaiton <- renderPlot({
  single_ionsplot_save()
})
})

correlated_ions_m <- eventReactive(c(input$start_visualisation_analysis), {
  data_rds<-data_rds()
  spatial_coord <- spatial_coord()
  req(spatial_coord)
  req(data_rds)
  req(input$select_single_ions)
  withProgress(message = "Processing data...",value=0.8,{
    data_mrds<-data_rds[[1]]
    m = as.character(rownames(data_mrds@assays$Spatial$counts))

    ion <- input$select_single_ions
    m <- subset(m, m !=ion)
    spatial_coord_m<-spatial_coord |>
      dplyr::select(x,y,all_of(m))
    
    w <- which(names(spatial_coord) == ion)
    ion_intensity <- spatial_coord[,w]

    corr_datalist <- lapply(colnames(spatial_coord_m)[-c(1:2)], function(x) {
      test_result <- cor.test(spatial_coord_m[,x], ion_intensity, method = "spearman",exact = FALSE)
      c(feature=x,cor = test_result$estimate, p = test_result$p.value)
    })
    corr_data<-data.frame(do.call(rbind,corr_datalist))
    colnames(corr_data)<-c("feature","cor","p")
    correlated_ions_m <- corr_data |>
      arrange(as.numeric(cor)) |>
      mutate(cor=round(as.numeric(cor),3),
             p=round(as.numeric(p),3))

    return(correlated_ions_m)
  })
})
correlated_ions_t <- eventReactive(c(input$start_visualisation_analysis), {
  data_rds<-data_rds()
  spatial_coord <- spatial_coord()
  req(spatial_coord)
  req(data_rds)
  req(input$select_single_ions)
  withProgress(message = "Processing data...",value=0.8,{
    data_trds<-data_rds[[2]]
    g = as.character(rownames(data_trds@assays$Spatial$counts))

    ion <-input$select_single_ions
    g <- subset(g, g !=ion)
    spatial_coord_t<-spatial_coord |>
      dplyr::select(x,y,all_of(g))
    
    w <- which(names(spatial_coord) == ion)
    ion_intensity <- spatial_coord[,w]

    corr_datalist <- lapply(colnames(spatial_coord_t)[-c(1:2)], function(x) {
      test_result <- cor.test(spatial_coord_t[,x], ion_intensity, method = "spearman",exact = FALSE)
      c(feature=x,cor = test_result$estimate, p = test_result$p.value)
    })
    corr_data<-data.frame(do.call(rbind,corr_datalist))
    colnames(corr_data)<-c("feature","cor","p")
    correlated_ions_t <- corr_data |>
      arrange(as.numeric(cor)) |>
      mutate(cor=round(as.numeric(cor),3),
             p=round(as.numeric(p),3))

    return(correlated_ions_t)
  })
})
positively_correlated_ions_m <- eventReactive(c(input$start_visualisation_analysis), {
  correlated_ions <- correlated_ions_m()
  spatial_coord <- spatial_coord()
  req(spatial_coord)
  req(correlated_ions)
  correlated_ions<- correlated_ions |>
    dplyr::filter(cor>0)
  if(nrow(correlated_ions)==0){
    return(NULL)
  }
  k <- tail(correlated_ions)
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
      scale_color_gradientn(colours = heatmap_Palette(100)) +
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
  posdata=spatial_coord %>%
    dplyr::select(x,y,all_of(k$feature))
  #pp<-list(p1,p2,p3,p4,p5,p6)
  positively_correlated_ions_m<-list(plot_list,posdata)
  return(positively_correlated_ions_m)
})
observeEvent(input$start_visualisation_analysis,{
output$positively_correlated_ions_m_plot <- renderPlot({
  req(positively_correlated_ions_m())
  withProgress(message = "Plotting...", value = 0.8, {
    p<-positively_correlated_ions_m()[[1]]
    do.call(gridExtra::grid.arrange, c(p, ncol = 3))
    #gridExtra::grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],ncol=3)
  })
})
})

negatively_correlated_ions_m <-eventReactive(c(input$start_visualisation_analysis), {
  correlated_ions <- correlated_ions_m()
  spatial_coord <- spatial_coord()
  req(spatial_coord)
  req(correlated_ions)
  correlated_ions<- correlated_ions |>
    dplyr::filter(cor<0)
  if(nrow(correlated_ions)==0){
    return(NULL)
  }
  k <- head(correlated_ions)
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
      scale_color_gradientn(colours = heatmap_Palette(100)) +
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

  }
  negdata=spatial_coord %>%
    dplyr::select(x,y,all_of(k$feature))
  negatively_correlated_ions_m<-list(plot_list,negdata)
  return(negatively_correlated_ions_m)
  
})
observeEvent(input$start_visualisation_analysis,{
output$negatively_correlated_ions_m_plot <- renderPlot({
  req( negatively_correlated_ions_m())
  withProgress(message = "Plotting...", value = 0.8, {
    p<-negatively_correlated_ions_m()[[1]]
    do.call(gridExtra::grid.arrange, c(p, ncol = 3))
  })
})
})
positively_correlated_ions_t <-eventReactive(c(input$start_visualisation_analysis), {
  correlated_ions <- correlated_ions_t()
  spatial_coord <- spatial_coord()
  req(spatial_coord)
  req(correlated_ions)
  correlated_ions<- correlated_ions |>
    dplyr::filter(cor>0)
  if(nrow(correlated_ions)==0){
    return(NULL)
  }
  k <- tail(correlated_ions)
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
      scale_color_gradientn(colours = heatmap_Palette(100)) +
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
    
  }
  posdata=spatial_coord %>%
    dplyr::select(x,y,all_of(k$feature))
  positively_correlated_ions_t<-list(plot_list,posdata)
  return(positively_correlated_ions_t)
})
observeEvent(input$start_visualisation_analysis,{
output$positively_correlated_ions_t_plot <- renderPlot({
  req( positively_correlated_ions_t())
  withProgress(message = "Plotting...", value = 0.8, {
    p<-positively_correlated_ions_t()[[1]]
    do.call(gridExtra::grid.arrange, c(p, ncol = 3))
  })
})
})

negatively_correlated_ions_t <- eventReactive(c(input$start_visualisation_analysis), {
  correlated_ions <- correlated_ions_t()
  spatial_coord <- spatial_coord()
  req(spatial_coord)
  req(correlated_ions)
  correlated_ions<- correlated_ions |>
    dplyr::filter(cor<0)
  if(nrow(correlated_ions)==0){
    return(NULL)
  }
  k <- head(correlated_ions)
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


    plot_list[[i]]  <- ggplot(data, aes(x = x, y = y)) +
      geom_point(aes(color = norm_intensity), size = 1) +
      scale_color_gradientn(colours = heatmap_Palette(100)) +
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

  }
  negdata=spatial_coord %>%
    dplyr::select(x,y,all_of(k$feature))
  negatively_correlated_ions_t<-list(plot_list,negdata)
  return(negatively_correlated_ions_t)
})

observeEvent(input$start_visualisation_analysis,{
output$negatively_correlated_ions_t_plot <-renderPlot({
  req( negatively_correlated_ions_t())
  withProgress(message = "Plotting...", value = 0.8, {
    p<- negatively_correlated_ions_t()[[1]]
    do.call(gridExtra::grid.arrange, c(p, ncol = 3))
  })
})
})
