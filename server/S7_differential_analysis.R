# ==============================================================================
# S7_differential_analysis.R
# Server logic for Step 4: Differential Analysis
# Handles differential testing, thresholding, and visualization (Volcano, Barplot, UMAP)
# ==============================================================================

# ------------------------------------------------------------------------------
# Download Handlers
# ------------------------------------------------------------------------------

#' @title Download Differential Analysis Results
#' @description Exports the full table of differential genes and metabolites as a ZIP file.
#' Includes p-values, logFC, and significance status.
output$download_diff_omics <- downloadHandler(
  filename = function() {

    paste0("diff_data.zip")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
    diff_omics<-diff_omics()
    # Rename internal column names for user-friendly export
    diff_m<-diff_omics[[1]] %>%
      dplyr::rename(log_p_val_adj=log_test) %>%
      dplyr::select(-avg_log2FC)
    diff_t<-diff_omics[[2]] %>%
      dplyr::rename(log_p_val_adj=log_test) %>%
      dplyr::select(-avg_log2FC)
    
    diff_omics<-list(diff_m,diff_t)
    tempdir <- setwd(tempdir())
    on.exit(setwd(tempdir))
    fi=c("Metabolite_diff_data.csv","Gene_diff_data.csv")
    for (i in 1:length(diff_omics)) {
      data<-diff_omics[[i]]
      write.csv(data, fi[i], row.names = FALSE)
    }
    zip(file,fi)
    })
  })

#' @title Download Differential Bar Plot Data
#' @description Exports the source data for the differential count bar plot.
output$download_diff_count_all <- downloadHandler(
  filename = function() {

    paste0("diff_barplot.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    diff_count_all<- diff_count_all()
    diff_m_plotdata<-diff_count_all[[1]]
    diff_t_plotdata<-diff_count_all[[2]]
    diff_count<-diff_m_plotdata[[2]]
    diff_count_t<-diff_t_plotdata[[2]]
    DemoData<-rbind(diff_count,diff_count_t)
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })

#' @title Download Differential Bar Plot Image
output$download_diffbar_plot <- downloadHandler(
  filename = function() {
    "diffbar_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- diffbar_plot_save()

    ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)
  })
})

# Volcano Plot Downloads
output$download_volcano_data <- downloadHandler(
  filename = function() {

    paste0("Metabolite_volcano_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
    volcano_data<-volcano_data()
    DemoData <-volcano_data[[1]][[2]]
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_volcano_plot <- downloadHandler(
  filename = function() {
    "Metabolite_volcano_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      volcano_plot_save<-volcano_plot_save()
      p<-volcano_plot_save[[1]]

    ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)
    })
  })
output$download_volcano_datat <- downloadHandler(
  filename = function() {

    paste0("Gene_volcano_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
      volcano_data<-volcano_data()
      DemoData <-volcano_data[[2]][[2]]
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_volcano_plott <- downloadHandler(
  filename = function() {
    "Gene_volcano_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      volcano_plot_save<-volcano_plot_save()
      p<-volcano_plot_save[[2]]
  
      ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)
    })
  })

# Spatial Distribution Downloads
output$download_diffion_distribution_plot <- downloadHandler(
  filename = function() {
    "Metabolite_diff_distribution_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- diffion_distribution_plot_save()

    ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)
    })
  })
output$download_diffion_distribution_data <- downloadHandler(
  filename = function() {
    paste0("Metabolite_diff_distribution_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData <- diffion_distribution_plotdata_m() |>
      dplyr::select(x,y,intensity)
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_diffion_distribution_plot_t <- downloadHandler(
  filename = function() {
    "Gene_diff_distribution_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- diffion_distribution_plot_save_t()

    ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)
  })
})
output$download_diffion_distribution_data_t <- downloadHandler(
  filename = function() {
    paste0("Gene_diff_distribution_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    DemoData<-diffion_distribution_plotdata_t() |>
      dplyr::select(x,y,intensity)
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })


# UMAP Downloads
output$download_umap_plot <- downloadHandler(
  filename = function() {
    "Metabolite_umap_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      umap_plot_data<-umap_plot_data()
      p<-umap_plot_data[[1]][[1]]

    ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)

  })
})
output$download_umap_data <- downloadHandler(
  filename = function() {
    paste0("Metabolite_umap_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
    umap_plot_data<-umap_plot_data()
    umap_data<-umap_plot_data[[2]]
    DemoData<-umap_data[[1]] 
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_umap_plott <- downloadHandler(
  filename = function() {
    "Gene_umap_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      umap_plot_data<-umap_plot_data()
      p<-umap_plot_data[[1]][[2]]
  
      ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)
      
    })
  })
output$download_umap_datat <- downloadHandler(
  filename = function() {
    paste0("Gene_umap_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
      umap_plot_data<-umap_plot_data()
      umap_data<-umap_plot_data[[2]]
      DemoData<-umap_data[[2]] 
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
      
    })
  })

#' @title Run Initial Differential Test
#' @description Executes Seurat::FindMarkers (Wilcoxon test) on selected groups.
#' This step is computationally intensive and runs before thresholding.
#' @return Raw differential test results (p-values, logFC) for both omics.
pre_diff_omics<- eventReactive(c(input$start_diff_analysis), {
  req(data_rds_group()) # Requires data with defined groups (Treatment vs Control)
  withProgress(message = "Calculating Differential Statistics...",value=0.8,{
    data_rds_group_data <- data_rds_group()

    data_mrds_group <- data_rds_group_data[[1]]
    data_trds_group <- data_rds_group_data[[2]]
    group <- c("treatment", "control") 
    
    source("./source/OverallAnalysisFunction/Findmarker/find_marker.R")
    
    # Calculate differential features for Metabolites and Genes
    pdiff_m <- find_marker(data_mrds_group, group)#mz
    pdiff_t <- find_marker(data_trds_group,  group)
    
    pre_diff_omics<-list(pdiff_m,pdiff_t)
    gc() # Garbage collection to free memory
    return(pre_diff_omics)
    
  })
})

#' @title Classify Differential Status
#' @description Applies fold change and p-value thresholds to classify features as Up/Down/NS.
#' Uses user inputs: Fold_change, p_value_adjust.
diff_omics_allstatus<-eventReactive(c(input$start_diff_analysis), {
  req(pre_diff_omics())
  withProgress(message = "Processing data...",value=0.8,{
    pre_diff_omics <- pre_diff_omics()
    pdiff_m <- pre_diff_omics[[1]]
    pdiff_t <- pre_diff_omics[[2]]
    group <- c("treatment", "control") 
    
    # Get thresholds from UI
    FC_Threshold <- 2^input$Fold_change
    FC_Threshold_t <- 2^input$Fold_change_t
    pvalue <- input$p_value_adjust
    pvalue_t <- input$p_value_adjust_t
    
    source("./source/OverallAnalysisFunction/Findmarker/find_marker.R")
    
    # Classify features
    diff_m <- find_marker_status(pdiff_m, type = "metabolite", group, FC_Threshold, pvalue)#mz
    diff_t <- find_marker_status(pdiff_t, type = "gene", group, FC_Threshold_t, pvalue_t)
    diff_omics_allstatus<-list(diff_m,diff_t)
    return(diff_omics_allstatus)
    
  })
})

# Observer to validate differential results count
 observeEvent(input$start_diff_analysis, {
  req(diff_omics_allstatus())
  withProgress(message = "Processing data...", value = 0.8, {
    diff_omics_allstatus<-diff_omics_allstatus()
    diff_m <- diff_omics_allstatus[[1]]
    diff_t <- diff_omics_allstatus[[2]]

    diff_num_t<-diff_t %>%
      dplyr::filter(State!="Non-significant")
    if (nrow(diff_num_t) < 3) {
      showNotification("The number of differential gene is less than three. please modify the threshold.")
    }
    diff_num_m<-diff_m %>%
      dplyr::filter(State!="Non-significant")
    if (nrow(diff_num_m) < 3) {
      showNotification("The number of differential metabolites is less than three. Please adjust the threshold or select other preprocessing methods (e.g., TIC/RMS + LogNormalize).")
    }

    if(nrow(diff_num_t) >= 3 && nrow(diff_num_m) >= 3) {
    diff_omics(list(diff_m = diff_m, diff_t = diff_t))

    }else{
      diff_omics(NULL)
    }
  })
})

#' @title Generate Summary Counts for Bar Plot
#' @description Summarizes the number of significant features for visualization.
diff_count_all<-eventReactive(c(input$start_diff_analysis), {
  req(diff_omics_allstatus())
  withProgress(message = "Summarizing Results...",value=0.8,{
  diff_omics<-diff_omics_allstatus()
  diff_m<-diff_omics[[1]]
  diff_t<-diff_omics[[2]]
  group<-c("treatment","control")
  
  source("./source/OverallAnalysisFunction/Diffbarplot/diffbar_data.R")
  
  # Calculate counts for stacked bar plot
  diff_m_plotdata<-diffbar_data(diff_m,group,"Metabolite")
  diff_t_plotdata<-diffbar_data(diff_t,group,"Gene")
  
  diff_count_all<-list(diff_m_plotdata,diff_t_plotdata)
  return(diff_count_all)
  })
})


#' @title Aggregate Differential Feature Counts
#' @description Summarizes the number of significant differential features (Up/Down) for both metabolites and genes.
#' Reshapes the data into a wide format suitable for tabular display or further analysis.
#' @return A data frame with counts of differential features per sample type.
DifferNumber<-eventReactive(c(input$start_diff_analysis), {
  req(diff_count_all())
  withProgress(message = "Processing data...",value=0.8,{
    diff_count_all<-diff_count_all()
    data_count_long<-diff_count_all[[1]][[1]]
    data_count_long_t<-diff_count_all[[2]][[1]]
    data_count_long_all<-rbind(data_count_long,data_count_long_t)
    
    data_count_long_all1<-data_count_long_all %>% group_by(Sample,type) %>% 
      dplyr::summarise(total_MIDCount = sum(value)) %>%
      pivot_wider(
        names_from = type,
        values_from = total_MIDCount,
        values_fill = 0
      ) %>%
      as.data.frame()
    return(data_count_long_all1)
  })
})

#' @title Generate Differential Count Bar Plot
#' @description Creates a bar plot visualizing the number of Up and Down regulated features 
#' for both Transcriptomics and Metabolomics.
#' @return A ggplot object showing differential counts.
diffbar_plot_save<-eventReactive(c(input$start_diff_analysis), {
  req(diff_count_all())
  withProgress(message = "Processing data...",value=0.8,{
    diff_count_all<-diff_count_all()
    data_count_long<-diff_count_all[[1]][[1]]
    data_count_long_t<-diff_count_all[[2]][[1]]
    
    data_count_long_all<-rbind(data_count_long,data_count_long_t)
    
    # Check for empty data
    if(nrow(data_count_long_all) == 0) return(NULL)
    
    data_count_long_all %<>% filter(value != 0)
    if(nrow(data_count_long_all) == 0) return(NULL) # Check again after filter
    
    UD_Class <- unique(data_count_long_all$State)
    UD_len <- length(UD_Class)
    
    scale_fill_values <- "grey" # Default
    scale_fill_labels <- "None" # Default
    
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
    data_count_long_all %<>% dplyr::filter(value != 0)
    data_count_long_all$State <- factor(data_count_long_all$State, levels = c('Up', 'Down'))
    data_count_long_all$Sample_type<-data_count_long_all$type
    
    p <- ggplot(data_count_long_all,mapping = aes(x=Sample_type, y=value, fill=State)) +
      theme_bw()+ theme(panel.grid=element_blank())+
      geom_bar(stat="identity",width = 0.5, position = position_dodge(0.7)) +
      geom_text(aes(label=value),hjust = 0.3,vjust = 1,position = position_dodge(0.7))+
      scale_fill_manual(values = scale_fill_values, labels = scale_fill_labels)+
      labs(x="",y="")+
      theme(axis.text.x= element_text(angle = 70 ,hjust =1),
            text = element_text(size=15))
    return(p)
  })
})
#' @title Render Differential Count Bar Plot
#' @description Displays the bar plot generated by `diffbar_plot_save` in the UI.
observeEvent(input$start_diff_analysis,{
  req(diffbar_plot_save())
  diffbar_plot_save <- diffbar_plot_save()
  output$diffbar_plot <- renderPlot({
    withProgress(message = "Plotting...",value=0.8,{
      diffbar_plot_save()
    })
  })
})

#' @title Prepare Volcano Plot Data
#' @description Processes differential analysis results for metabolites and genes to fit the volcano plot format.
#' Adds necessary statistical columns and labels.
volcano_data <- eventReactive(c(input$start_diff_analysis), {
  req(diff_omics_allstatus())
  withProgress(message = "Processing data...",value=0.8,{
    diff_omics<-diff_omics_allstatus()
    diff_m<-diff_omics[[1]]
    diff_t<-diff_omics[[2]]
    source("./source/OverallAnalysisFunction/Volcano/volcano_data_processing.R")
    volcano_datam<-volcano_data_processing(diff_m)
    volcano_datat<-volcano_data_processing(diff_t)
    volcano_data<-list(volcano_datam,volcano_datat)
    
    return(volcano_data)
  })
})

#' @title Generate Volcano Plots
#' @description Creates volcano plots for both Metabolomics and Transcriptomics.
#' Highlights significant features based on user-defined Fold Change and P-value thresholds.
#' @return A list containing two ggplot objects (Metabolite and Gene).
volcano_plot_save <- eventReactive(c(input$start_diff_analysis), {
  req(volcano_data())
  withProgress(message = "Processing data...",value=0.8,{
    volcano_data<-volcano_data()
    volcano_data_m<-volcano_data[[1]]
    volcano_data_t<-volcano_data[[2]]
    
    group<-c("treatment","control")
    FC_Threshold=2^input$Fold_change
    FC_Threshold_t=2^input$Fold_change_t
    pvalue<-input$p_value_adjust
    pvalue_t<-input$p_value_adjust_t
    source("./source/OverallAnalysisFunction/Volcano/volcano_data_processing.R")
    volcano_datam<-volcano_plot_processing(volcano_data_m[[1]],"Metabolite",group,pvalue,FC_Threshold)
    volcano_datat<-volcano_plot_processing(volcano_data_t[[1]],"Gene",group,pvalue,FC_Threshold)
    plot<-list(volcano_datam,volcano_datat)
    
    return(plot)
  })
})
#' @title Render Metabolomics Volcano Plot
#' @description Displays the volcano plot for metabolomics data.
observeEvent(input$start_diff_analysis,{
  output$volcano_plot <- renderPlot({ 
    if(all(!is.null(volcano_plot_save()))){
      withProgress(message = "Plotting...",value=0.8,{
        volcano_plot_save<-volcano_plot_save()
        p<-volcano_plot_save[[1]]
        p
      })
    }else{
      NULL
    }
  })
})
#' @title Render Transcriptomics Volcano Plot
#' @description Displays the volcano plot for transcriptomics data.
observeEvent(input$start_diff_analysis,{
  output$volcano_plott <- renderPlot({  
  if(all(!is.null(volcano_plot_save()))){
    withProgress(message = "Plotting...",value=0.8,{
      volcano_plot_save<-volcano_plot_save()
      p<-volcano_plot_save[[2]]
      p
    })
  }else{
    NULL
  }
})
})

#' @title Update Differential Metabolite Selectize
#' @description Populates the dropdown with significant differential metabolites.
observeEvent(input$start_diff_analysis, {
  req(diff_omics_allstatus())
  diff_omics<-diff_omics_allstatus()
  diff_m<-diff_omics[[1]] %>%
    filter(State!="Non-significant")
  if(all(!is.na(diff_m))){
    ions <- sort(as.character(diff_m$metabolite))
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      inputId = "diff_ion_select",
      choices = ions,
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
    
  }else{
    
    updateSelectizeInput(session = getDefaultReactiveDomain(), "diff_ion_select", choices = NA,server = TRUE,selected = NULL)
  }
})
observeEvent(input$start_diff_analysis, {
  req(diff_omics_allstatus())
  diff_omics<-diff_omics_allstatus()
  diff_t<-diff_omics[[2]] %>%
    filter(State!="Non-significant")
  if(all(!is.na(diff_t))){
    ions_t <-  sort(as.character(diff_t$gene))
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      inputId = "diff_ion_select_t",
      choices = ions_t,
      server = TRUE,
      selected = if (length(ions_t) > 0) ions_t[1] else NULL,
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
    
  }else{
    
    updateSelectizeInput(session = getDefaultReactiveDomain(), "diff_ion_select_t", choices = NA,server = TRUE,selected = NULL)
  }
})
#' @title Prepare Gene Spatial Data
#' @description Extracts spatial expression data for the selected differential gene.
#' Normalizes intensity for visualization.
diffion_distribution_plotdata_t <- eventReactive(c(input$diff_ion_select,input$diff_ion_select_t,input$start_diff_analysis), {
  req(input$diff_ion_select_t)
  withProgress(message = "Processing data...",value=0.8,{
    print("diffion_distribution_plotdata run")
    data_rds_group <- data_rds_group()
    plotdatalist<-list()
    #t
    data_t=data_rds_group[[2]]
    plotdata_t <- data_t@meta.data
    if(!is.null(input$diff_ion_select_t) && input$diff_ion_select_t!=""){
      ion_t <- as.character(input$diff_ion_select_t)
    }else{
      ion_t <- NA
    }
    if (!is.na(ion_t)) {
      # Use FetchData for safe retrieval
      intensity_df <- tryCatch({
        Seurat::FetchData(data_t, vars = ion_t, layer = "counts")
      }, error = function(e) return(NULL))
      
      if (!is.null(intensity_df) && ncol(intensity_df) > 0) {
        # Ensure row order matches if needed, though FetchData returns consistent rows
        plotdata_t$intensity <- intensity_df[[1]]
        plotdata_t$norm_intensity <- 100 * (plotdata_t$intensity) / max(plotdata_t$intensity, na.rm=TRUE)
      } else {
        plotdata_t <- NA
      }
    } else {
      plotdata_t <- NA
    }
    
    print("diffion_distribution_plotdata_t done")
    return(plotdata_t)
  }) 
}) 
#' @title Prepare Metabolite Spatial Data
#' @description Extracts spatial intensity data for the selected differential metabolite.
#' Normalizes intensity for visualization.
diffion_distribution_plotdata_m <- eventReactive(c(input$diff_ion_select,input$diff_ion_select_t,input$start_diff_analysis), {
  req(input$diff_ion_select)
  withProgress(message = "Processing data...",value=0.8,{
    print("diffion_distribution_plotdata run")
    data_rds_group <- data_rds_group()
    plotdatalist<-list()
    #m
    data_m=data_rds_group[[1]]
    plotdata_m <- data_m@meta.data
    if(!is.null(input$diff_ion_select) && input$diff_ion_select!=""){
      ion <- as.character(input$diff_ion_select)
    }else{
      ion <- NA
    }
    if (!is.na(ion)) {
      # Use FetchData for safe retrieval
      intensity_df <- tryCatch({
        Seurat::FetchData(data_m, vars = ion, layer = "counts")
      }, error = function(e) return(NULL))
      
      if (!is.null(intensity_df) && ncol(intensity_df) > 0) {
        plotdata_m$intensity <- intensity_df[[1]]
        plotdata_m$norm_intensity <- 100 * (plotdata_m$intensity) / max(plotdata_m$intensity, na.rm=TRUE)
      } else {
        plotdata_m <- NA
      }
    } else {
      plotdata_m <- NA
    }
    
    print("diffion_distribution_plotdata_m done")
    return(plotdata_m)
  }) 
})

#' @title Generate Metabolite Spatial Plot
#' @description Creates a spatial heatmap for the selected differential metabolite.
diffion_distribution_plot_save <- eventReactive(c(input$diff_ion_select,input$start_diff_analysis), {
  req(diffion_distribution_plotdata_m())
  withProgress(message = "Processing data...",value=0.8,{
    data <- diffion_distribution_plotdata_m()
    if(any(!is.na(data))){
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
    }else{
      p<-NULL
    }
    return(p)
  })
})
#' @title Render Metabolite Spatial Plot
#' @description Displays the spatial heatmap for the selected metabolite.
observeEvent(input$start_diff_analysis,{
  output$diffion_distribution_plot <- renderPlot({
    req(diffion_distribution_plot_save())
    withProgress(message = "Plotting...",value=0.8,{
      diffion_distribution_plot_save()
    })
  })
})

#' @title Generate Gene Spatial Plot
#' @description Creates a spatial heatmap for the selected differential gene.
diffion_distribution_plot_save_t <- eventReactive(c(input$diff_ion_select_t,input$start_diff_analysis), {
  req(diffion_distribution_plotdata_t())
  withProgress(message = "Processing data...",value=0.8,{
    data <- diffion_distribution_plotdata_t()
    if(any(!is.na(data))){
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
    }else{
      p<-NULL
    }
    return(p)
  })
})
observeEvent(input$start_diff_analysis,{
  output$diffion_distribution_plot_t <- renderPlot({
    req( diffion_distribution_plot_save_t())
    withProgress(message = "Plotting...",value=0.8,{
      diffion_distribution_plot_save_t()
    })
  })
})


#' @title Subset Differential Seurat Objects
#' @description Creates subsets of the Metabolomics and Transcriptomics Seurat objects containing 
#' only the samples from the defined Treatment and Control groups. 
#' Recalculates UMAP for visualization of these groups.
#' @return A list of subsetted Seurat objects.
diff_rds<-eventReactive(c(input$start_diff_analysis), {
  
  req(data_rds_group())
  withProgress(message = "Processing data...",value=0.8,{
    data_rds_group<-data_rds_group()
    data_mrds_group<-data_rds_group[[1]]
    data_trds_group<-data_rds_group[[2]]
    
    group<-c("treatment","control")
    diff_rds_m<-subset(data_mrds_group,groups %in% group)
    diff_rds_t<-subset(data_trds_group,groups %in% group)
    
    if(is.null(diff_rds_m@reductions$umap)){
      diff_rds_m<-RunUMAP(diff_rds_m, dims = 1:30)
    }
    #umap
    if(is.null(diff_rds_t@reductions$umap)){
      diff_rds_t<-RunUMAP(diff_rds_t, dims = 1:30)
    }
    
    source("./source/main_program/slim_seurat.R")
    diff_rds_m <- slim_seurat(seurat_obj=diff_rds_m ,keep_reductions = TRUE)
    diff_rds_t <- slim_seurat(seurat_obj=diff_rds_t, keep_reductions = TRUE)
    
    diff_rds<-list(diff_rds_m,diff_rds_t)
    
    return(diff_rds)
  })
})


#' @title Generate Differential UMAP Plots
#' @description Creates UMAP visualizations specifically for the Treatment vs Control comparison.
#' Helps visualize the separation between the groups in feature space.
#' @return A list containing UMAP plots and the underlying embedding data.
umap_plot_data <-eventReactive(c(input$start_diff_analysis), {
  req(diff_rds())
  withProgress(message = "Processing data...",value=0.8,{
    diff_rds<-diff_rds()
    diff_rds_m<-diff_rds[[1]]
    diff_rds_t<-diff_rds[[2]]
    source("./source/OverallAnalysisFunction/Clustering/clusterplot.R")
    umap = diff_rds_m@reductions$umap@cell.embeddings %>%
      as.data.frame()
    umapdata_m<-data.frame(x=diff_rds_m$x,y=diff_rds_m$y,umap,group=diff_rds_m$groups)
    umap = diff_rds_t@reductions$umap@cell.embeddings %>%
      as.data.frame()
    umapdata_t<-data.frame(x=diff_rds_t$x,y=diff_rds_t$y,umap,group=diff_rds_t$groups)
    
    umapdata<-list(umapdata_m,umapdata_t)
    
    plotm <- umap_plot(diff_rds_m,pt.size=1) + ggtitle("Metabolite")
    plott <- umap_plot(diff_rds_t,pt.size=1) + ggtitle("Gene")
    
    plot<-list(plotm,plott)
    umap_plot_data<-list(plot,umapdata)
    return(umap_plot_data)
  })
})

observeEvent(input$start_diff_analysis,{
  req(umap_plot_data())
  umap_plot_data<-umap_plot_data()
  umap_plot<-umap_plot_data[[1]]
  output$umap_plot <- renderPlot({
    withProgress(message = "Plotting...",value=0.8,{
      p<-umap_plot[[1]]
      p
    })
  })
})
#' @title Render Gene UMAP
#' @description Displays the UMAP plot for the transcriptomics comparison.
observeEvent(input$start_diff_analysis,{
  req(umap_plot_data())
  umap_plot_data<-umap_plot_data()
  umap_plot<-umap_plot_data[[1]]
  output$umap_plott <- renderPlot({
    withProgress(message = "Plotting...",value=0.8,{
      p<-umap_plot[[2]]
      p
    })
  })
})

