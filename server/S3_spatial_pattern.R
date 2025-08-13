


output$download_spatial_spatial_pattern_polt_m <- downloadHandler(
  filename = function() {
    "Metabolite_spatial_pattern.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    spatial_pattern_plot_save<-spatial_pattern_plot_save()
    p <- spatial_pattern_plot_save[[1]]

    ggsave(file, plot = p, device = "png", width = 12, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })

output$download_spatial_spatial_pattern_polt_t <- downloadHandler(
  filename = function() {
    "Gene_spatial_pattern.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    spatial_pattern_plot_save<-spatial_pattern_plot_save()
    p <- spatial_pattern_plot_save[[2]]

    ggsave(file, plot = p, device = "png", width = 12, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })

##
output$download_spatial_heatmap_data <- downloadHandler(
  filename = function() {

    paste0("spatial_heatmap_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
      spatial_heatmap_data<-spatial_heatmap_data()
           DemoData<-as.data.frame(spatial_heatmap_data[[1]])
           write.table(DemoData,file,row.names = T, quote = F, sep = '\t')
    })
  })

output$download_spatial_heatmap_plot <- downloadHandler(
  filename = function() {
    "spatial_heatmap_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      spatial_heatmap_data=spatial_heatmap_data()
      cmt<-spatial_heatmap_data[[1]]
      annotation_col<-spatial_heatmap_data[[2]]
      p<-spatial_heatmap_plotsave()
      ggsave(file, plot = p, device = "png", width = 10, height = 10,bg = "#FFFFFF", dpi = 300)

    })
  })

###
find_spatial_pattern<- eventReactive(c(input$start_Spatial_Pattern), {
  req(data_rds())
  withProgress(message = "Processing data...",value=0.8,{
    data_rds<-data_rds()
    data_mrds<-data_rds[[1]]
    data_trds<-data_rds[[2]]

    spatialpattern<-function(rds){
      expr=rds@assays$Spatial$counts
      location=data.frame(x=rds@meta.data$x,y=rds@meta.data$y)
      rownames(location)<-paste0("sample:",location$x,"_",location$y)
      k=SpaGene(expr,location)

      parttern=FindPattern(k)
      data<-list(parttern,location,k$spagene_res)
      return(data)
    }
    
    data_m_parttern<-spatialpattern(rds=data_mrds)
    data_t_parttern<-spatialpattern(rds=data_trds)

    find_spatial_pattern<-list(data_m_parttern,data_t_parttern)

    return(find_spatial_pattern)
  })
})
spatial_pattern_plot_save<- eventReactive(c(input$start_Spatial_Pattern), {
  req(find_spatial_pattern())
  withProgress(message = "Plotting...",value=0.8,{
  find_spatial_pattern <- find_spatial_pattern()

  data_m_parttern<-find_spatial_pattern[[1]]
  data_t_parttern<-find_spatial_pattern[[2]]

source("./source/OverallAnalysisFunction/spatialpattern/spatialpattern.R")
  mplot<-PlotPattern_c(data_m_parttern[[1]],data_m_parttern[[2]])
  tplot<-PlotPattern_c(data_t_parttern[[1]],data_t_parttern[[2]])

  spatial_pattern_plot_save<-list(mplot,tplot)
  return(spatial_pattern_plot_save)
  })
  })
    
observeEvent(input$start_Spatial_Pattern,{
  req(spatial_pattern_plot_save())
  spatial_pattern_plot_save<-spatial_pattern_plot_save()
  output$spatial_pattern_plot_m <- renderPlot({
    withProgress(message = "Plotting...",value=0.8,{
      spatial_pattern_plot_save[[1]]
      
    })
  })
  output$spatial_pattern_plot_t <- renderPlot({
    withProgress(message = "Plotting...",value=0.8,{
      spatial_pattern_plot_save[[2]]
      
    })
  })

})
spatial_heatmap_datalist<- eventReactive(c(input$start_Spatial_Pattern), {
  req(find_spatial_pattern())
  withProgress(message = "Processing data...",value=0.8,{
    find_spatial_pattern <- find_spatial_pattern()
    data_m_parttern<-find_spatial_pattern[[1]]
    data_t_parttern<-find_spatial_pattern[[2]]

    source("./source/OverallAnalysisFunction/spatialpattern/spatialpattern.R")
    mdata<-Patterndata_c(data_m_parttern[[1]],data_m_parttern[[2]],colname="Metabolite")
    tdata<-Patterndata_c(data_t_parttern[[1]],data_t_parttern[[2]],colname="Gene")

    merged_df <- Reduce(function(x, y) full_join(x, y, by = "location"), list(mdata,tdata))

 return(list(mdata,tdata,merged_df))
  })
})

spatial_heatmap_data<- eventReactive(c(input$start_Spatial_Pattern), {
  req(spatial_heatmap_datalist())
  withProgress(message = "Processing data...",value=0.8,{
    spatial_heatmap_datalist <- spatial_heatmap_datalist()
    df=spatial_heatmap_datalist[[3]] %>%
      tidyr::separate(col = location,
                      into = c("prefix", "x", "y"),
                      sep = "[:\\.]|_") %>%
      dplyr::select(-prefix) %>%
      dplyr::mutate(x = as.numeric(x), y = as.numeric(y))
    source("./source/CorrelationAnalysisFunction/Correlation/run_moran.R")
    moran_matrix<-run_moran(df)[[1]]

    annotation_col <- data.frame(
      Group = factor(ifelse(grepl("^Metabolite", rownames(moran_matrix)), "Metabolite", "Gene"))
    )
    rownames(annotation_col) <- colnames(moran_matrix)

    return(list(moran_matrix,annotation_col))
  })
})

spatial_heatmap_plotsave<- eventReactive(c(input$start_Spatial_Pattern), {
  req(spatial_heatmap_data())
  withProgress(message = "Processing data...",value=0.8,{
    spatial_heatmap_datalist <- spatial_heatmap_datalist()
    data_list <- spatial_heatmap_data()

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

    plot_obj<- pheatmap::pheatmap(
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
      silent = TRUE  
    )
    heatmap_gg <- ggplot() +
      annotation_custom(plot_obj$gtable) +
      theme_void() 

    return(heatmap_gg)
  })
})
output$spatial_heatmap_plot <- renderPlot({
  req(spatial_heatmap_plotsave())
  
  withProgress(message = "Plotting heatmap...", value = 0.8, {

    spatial_heatmap_plotsave()
  })
})  

output$download_Pattern_specific_data_m <- downloadHandler(
  filename = function() {

    paste0("Metabolite_Pattern_specific.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      Pattern_data<-Pattern_data()
      DemoData<-Pattern_data[[1]]
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_Pattern_specific_polt_m <- downloadHandler(
  filename = function() {
    "Metabolite_Pattern_specific.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      Pattern_data<-Pattern_specific_plotsave()
      p <- Pattern_data[[1]]
  
      ggsave(file, plot = p, device = "png", width = 8, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })
output$download_Pattern_specific_data_t <- downloadHandler(
  filename = function() {

    paste0("Gene_Pattern_specific.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      Pattern_data<-Pattern_data()
      DemoData<-Pattern_data[[2]]
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_Pattern_specific_polt_t <- downloadHandler(
  filename = function() {
    "Gene_Pattern_specific.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      Pattern_specific<-Pattern_specific_plotsave()
      p <- Pattern_specific[[2]]
  
      ggsave(file, plot = p, device = "png", width = 8, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })
##

output$download_Pattern_specific_m_plot <- downloadHandler(
  filename = function() {
    "Pattern_specific_metabolite_distribution_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      p <- patternspecific_distribution_plot_save_m()
  
      ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })
output$download_Pattern_specific_m_data <- downloadHandler(
  filename = function() {
    paste0("Pattern_specific_metabolite_distribution_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData <- patternspecific_distribution_plotdata_m()
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_Pattern_specific_t_plot <- downloadHandler(
  filename = function() {
    "Pattern_specific_gene_distribution_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      p <- patternspecific_distribution_plot_save_t()
  
      ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })
output$download_Pattern_specific_t_data <- downloadHandler(
  filename = function() {
    paste0("Pattern_specific_gene_distribution_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData <- patternspecific_distribution_plotdata_t()
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })

Pattern_data<- reactive({
  req(find_spatial_pattern())
  find_spatial_pattern<-find_spatial_pattern()
  data_m_parttern<-find_spatial_pattern[[1]]
  data_t_parttern<-find_spatial_pattern[[2]]
  

  Pattern_specific_m<-as.data.frame(data_m_parttern[[1]]$genepattern)#[, pattern_select_m])
  Pattern_specific_t<-as.data.frame(data_t_parttern[[1]]$genepattern)#[, pattern_select_t])
  pattern_cutoff <- 0.5
  print(pattern_cutoff)
  run_pattern_class<- function(Pattern_specific,pattern_cutoff){
    max_cols <- apply(Pattern_specific, 1, function(x) {
      max_val <- max(x)          
      if (max_val > pattern_cutoff) {       
        colnames(Pattern_specific)[which.max(x)]  
      } else {
        "unclassified"           
      }
    })

    max_v <- apply(Pattern_specific, 1, function(x) {
      max_val <- max(x)          
      if (max_val > pattern_cutoff) {
        max_val
      } else {
        NA       
      }
    })
    pattern_class <- data.frame(
      feature = rownames(Pattern_specific),
      class = max_cols,
      similarity=max_v,
      stringsAsFactors = FALSE
    )
    return(pattern_class)
  }
  pattern_class_m<- run_pattern_class(Pattern_specific=Pattern_specific_m,pattern_cutoff=pattern_cutoff) %>%
    arrange(desc(similarity))
  pattern_class_t<- run_pattern_class(Pattern_specific=Pattern_specific_t,pattern_cutoff=pattern_cutoff) %>%
    arrange(desc(similarity))
  Pattern_data<-list(pattern_class_m,pattern_class_t)
  return(Pattern_data)
})
Pattern_specific_plotsave<-reactive({
  req(Pattern_data())
  pattern_class <- Pattern_data()

  pattern_class_m<- pattern_class[[1]]
  pattern_count_m <- table(pattern_class_m$class)
  pattern_count_df_m <- as.data.frame(pattern_count_m)
  colnames(pattern_count_df_m) <- c("Pattern", "Count")
  p1 <- ggplot(pattern_count_df_m, aes(x = Pattern, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Pattern Metabolite Count", x = "Pattern", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pattern_class_t<- pattern_class[[2]]
  pattern_count_t <- table(pattern_class_t$class)
  pattern_count_df_t <- as.data.frame(pattern_count_t)
  colnames(pattern_count_df_t) <- c("Pattern", "Count")
  p2 <- ggplot(pattern_count_df_t, aes(x = Pattern, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Pattern Gene Count", x = "Pattern", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(list(p1,p2))
  
})

output$Pattern_specific_plot_m <- renderPlot({
  Pattern_specific_plotsave<-Pattern_specific_plotsave()
  req(Pattern_specific_plotsave)
  Pattern_specific_plotsave[[1]]
  
})

output$Pattern_specific_plot_t <- renderPlot({
  Pattern_specific_plotsave<-Pattern_specific_plotsave()
  req(Pattern_specific_plotsave)
  Pattern_specific_plotsave[[2]]
  
})



patternspecific_distribution_plotdata_m <-reactive({
  req(input$Pattern_specific_m_select)
  ion <- tryCatch(
    as.character(input$Pattern_specific_m_select),
    error = function(e) {
      showNotification("Invalid metabolite selection!", type = "error")
      return(NULL)
    }
  )
  req(ion, nzchar(ion)) 
  

  data_rds <- data_rds()
  data_m=data_rds[[1]]
  plotdata_m<-data.frame(x=data_m@meta.data$x,y=data_m@meta.data$y)
  plotdata_m$intensity <- data_m@assays$Spatial$counts[ion,]  # data_m@assays$Spatial$counts[ion,] # 某marker
  #norm
  plotdata_m$norm_intensity <-100*(plotdata_m$intensity)/max(plotdata_m$intensity)
  req(nrow(plotdata_m) > 0) 
  return(plotdata_m)
}) 
patternspecific_distribution_plot_save_m <-  reactive({
  data <- patternspecific_distribution_plotdata_m()
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))

  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = norm_intensity), size = 1) +
    scale_color_gradientn(colours = heatmap_Palette(100)) +
    xlim(min(data$x),max(data$x))+
    ylim(min(data$y),max(data$y))+
    coord_equal() +
    labs(title="Metabolite")+
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
  return(p)
})

output$patternspecific_distribution_plot <- renderPlot({
  req(patternspecific_distribution_plot_save_m())
  patternspecific_distribution_plot_save_m()
})


patternspecific_distribution_plotdata_t <-reactive({
  req(input$Pattern_specific_t_select)
  ion <- tryCatch(
    as.character(input$Pattern_specific_t_select),
    error = function(e) {
      showNotification("Invalid metabolite selection!", type = "error")
      return(NULL)
    })
  req(ion, nzchar(ion)) 

  data_rds <- data_rds()
  data_t=data_rds[[2]]
  plotdata_t<-data.frame(x=data_t@meta.data$x,y=data_t@meta.data$y)
  plotdata_t$intensity <- data_t@assays$Spatial$counts[ion,]  
  #norm
  plotdata_t$norm_intensity <-100*(plotdata_t$intensity)/max(plotdata_t$intensity)
  req(nrow(plotdata_t) > 0)  
  return(plotdata_t)
  
})
patternspecific_distribution_plot_save_t <- reactive({ 
  data <- patternspecific_distribution_plotdata_t()
  req(data)
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))

  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = norm_intensity), size = 1) +
    scale_color_gradientn(colours = heatmap_Palette(100)) +
    xlim(min(data$x),max(data$x))+
    ylim(min(data$y),max(data$y))+
    coord_equal()+
    labs(title="Gene")+
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
  return(p)
  
})
output$patternspecific_distribution_plot_t <- renderPlot({
  patternspecific_distribution_plot_save_t()
})


observe({
  pa <- input$pattern_select_m
  spatial_pattern_top <- Pattern_data()[[1]]

  req(
    spatial_pattern_top,
    pa,
    pa %in% spatial_pattern_top$class
  )

  spatial_pattern_mlist <- as.character(
    spatial_pattern_top$feature[spatial_pattern_top$class == pa]
  )
  pa_t<-input$pattern_select_t
  spatial_pattern_t <- Pattern_data()[[2]]

  req(
    spatial_pattern_t,
    pa_t,
    pa_t %in% spatial_pattern_t$class
  )

  spatial_pattern_tlist <- as.character(
    spatial_pattern_t$feature[spatial_pattern_t$class == pa_t]
  )
  
  if (length(spatial_pattern_mlist) > 300) {
    spatial_pattern_mlist <- spatial_pattern_mlist[1:300]
  } else if (length(spatial_pattern_mlist) < 3) {
    showNotification("The number of pattern-specific metabolite is less than three, please modify the threshold.")
  }
  
  if (length(spatial_pattern_tlist) > 300) {
    spatial_pattern_tlist <- spatial_pattern_tlist[1:300]
  } else if (length(spatial_pattern_tlist) < 3) {
    showNotification("The number of pattern-specific gene is less than three. please modify the threshold.")
  }
  if(length(spatial_pattern_tlist) >= 3 && length(spatial_pattern_mlist) >= 3) {

    spatial_pattern_ionlist(list(spatial_pattern_mlist, spatial_pattern_tlist))
  }else{
    spatial_pattern_ionlist(NULL)
  }
})

observe({
  pattern_class <- Pattern_data()[[1]]
  req(pattern_class)
  available_patterns <- unique(pattern_class$class)
  available_patterns <- available_patterns[available_patterns != "unclassified"]  
  available_patterns<-mixedsort(available_patterns)
  updateSelectInput(
    session = getDefaultReactiveDomain(),
    inputId = "pattern_select_m",
    choices = available_patterns,
    selected = if (length(available_patterns) > 0) available_patterns[1] else NULL
  )
  
})
observe({
  pattern_class <- Pattern_data()[[2]]
  req(pattern_class)
  available_patterns <- unique(pattern_class$class)
  available_patterns <- available_patterns[available_patterns != "unclassified"]  
  available_patterns<-mixedsort(available_patterns)
  updateSelectInput(
    session = getDefaultReactiveDomain(),
    inputId = "pattern_select_t",
    choices = available_patterns,
    selected = if (length(available_patterns) > 0) available_patterns[1] else NULL 
  )
})

observe({
  pa <- input$pattern_select_m
  spatial_pattern_top <- Pattern_data()[[1]]
  
  if (is.null(spatial_pattern_top) || is.null(pa) || !(pa %in% spatial_pattern_top$class)) {
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      "Pattern_specific_m_select",
      choices = character(0),
      selected = NULL
    )
    return()  
  }

  req(
    spatial_pattern_top,
    pa,
    pa %in% spatial_pattern_top$class
  )
  # 提取代谢物列表
  ions <- as.character(
    spatial_pattern_top$feature[spatial_pattern_top$class == pa]
  )
  ions<-sort(as.character(ions))
  updateSelectizeInput(
    session = getDefaultReactiveDomain(),
    inputId = "Pattern_specific_m_select",
    choices =ions,
    server = TRUE,
    selected = if (length(ions) > 0) ions[1] else NULL,
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

observe({
  pa <- input$pattern_select_t
  spatial_pattern_top <- Pattern_data()[[2]]
  
  if (is.null(spatial_pattern_top) || is.null(pa) || !(pa %in% spatial_pattern_top$class)) {
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      "Pattern_specific_t_select",
      choices = character(0),
      selected = NULL
    )
    return()  
  }

  req(
    spatial_pattern_top,
    pa,
    pa %in% spatial_pattern_top$class
  )

  ions <- as.character(
    spatial_pattern_top$feature[spatial_pattern_top$class == pa]
  )
  ions<-sort(as.character(ions))
  updateSelectizeInput(
    session = getDefaultReactiveDomain(),
    inputId = "Pattern_specific_t_select",
    choices =ions,
    server = TRUE,
    selected = if (length(ions) > 0) ions[1] else NULL,
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
