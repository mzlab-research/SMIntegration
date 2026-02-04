
# ==============================================================================
# S4_clustering.R
# Server logic for Step 3: Clustering Analysis and Cell Annotation
# Implements:
# - Step 3.1: Preprocessing & Data Assessment (Normalization, Integration, PCA)
# - Step 3.2: Clustering Algorithm Selection & Visualization (Louvain, SLM, K-means, Sankey)
# ==============================================================================

#' @title Download Cluster Plot
#' @description Handles the downloading of the combined clustering plot as a PNG file.
#' @return A PNG file named "cluster_plot.png".
output$download_cluster_plot  <- downloadHandler(
  filename = function() {
    "cluster_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- get("cluster_plot", envir = .GlobalEnv)

    ggsave(file, plot = p, device = "png", width = 24, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })

#' @title Download Cluster Data
#' @description Exports the clustering results for Metabolite, Gene, and Merged datasets.
#' @return A ZIP file containing three CSV files: "Metabolite_cluster_data.csv", "Gene_cluster_data.csv", and "Merge_cluster_data.csv".
output$download_cluster_data <- downloadHandler(
  filename = function() {
    "cluster_data.zip"
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
    cluster_plotlist<-cluster_plotlist()
    tempdir <- setwd(tempdir())
    on.exit(setwd(tempdir))
    fi=c("Metabolite_cluster_data.csv","Gene_cluster_data.csv","Merge_cluster_data.csv")
    for (i in 1:length(fi)) {
      data<-cluster_plotlist[[i]]
      write.csv(data[[2]], fi[i], row.names = FALSE)
    }

zip(file,fi)
    })
  })

#' @title Download Sankey Plot
#' @description Saves the interactive Sankey diagram as an HTML file.
#' @return An HTML file named "cluster_sanky_plot.html".
output$download_cluster_sanky_plot <- downloadHandler(
  filename = function() {
    "cluster_sanky_plot.html"
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
      p <- cluster_sanky_plot_save()
      saveNetwork(p,file)

      
    })
  })

#' @title Download Sankey Data
#' @description Exports the data used to generate the Sankey diagram.
#' @return A tab-separated text file named "cluster_sanky_data.txt".
output$download_cluster_sanky_plot_data <- downloadHandler(
  filename = function() {
    paste0("cluster_sanky_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData <- cluster_sanky()
      write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })



#' @title Preprocessing and Integration Reactive
#' @description This reactive block handles the preprocessing of metabolomics and transcriptomics data
#' based on user-selected parameters (normalization, transformation) and integrates them.
#' It returns a list containing processed Seurat objects for Metab, Trans, and Combined data.
rds_norm<- reactive({
  req(data_rds())#data_crds()
  withProgress(message = "Processing data...",value=0.8,{
    data_rds<-data_rds()

    # Apply preprocessing options to Metabolomics
    # Inputs from U2_upload.R: metab_norm_method, metab_transform_method, metab_scale
    norm_method <- if(!is.null(input$metab_norm_method)) input$metab_norm_method else "None"
    trans_method <- if(!is.null(input$metab_transform_method)) input$metab_transform_method else "LogNormalize"
    # Scaling is forced to TRUE for robust PCA/Clustering
    do_scale <- TRUE
    source("./source/preprocessing/run_prerds.R")
    # Process Metabolomics with custom options
    data_mrds <- RUNSCT(data_rds[[1]], norm_method = norm_method, transform_method = trans_method, scale_data = do_scale)
    
    # Process Transcriptomics with standard LogNormalize
    data_trds <- RUNSCT(data_rds[[2]], norm_method = "None", transform_method = "LogNormalize", scale_data = TRUE)

    # Perform Data Integration (Merge)
    decon_mtrx = data_mrds@assays$Spatial$counts
    decon_ttrx = data_trds@assays$Spatial$counts
    source("./source/preprocessing/pre_merge.R")
    data_crds<-run_merge_rds(decon_mtrx,decon_ttrx)
    
    # Process Combined Object

    data_crds<-run_prerds(data=data_crds)
    data_crds=RUNSCT(data_crds, norm_method = "None", transform_method = "LogNormalize", scale_data = TRUE)
    
    # Combine normalized data and scale data slots for integrated analysis
    data_crds@assays$SCT$data=rbind(
      data_mrds@assays$SCT$data,
      data_trds@assays$SCT$data
    )
    data_crds@assays$SCT$scale.data=rbind(
      data_mrds@assays$SCT$scale.data,
      data_trds@assays$SCT$scale.data
    )
    # Combine Variable Features
    hv_metab <- VariableFeatures(data_mrds)
    hv_trans <- VariableFeatures(data_trds)
    hv_combined <- c(hv_metab, hv_trans)
    VariableFeatures(data_crds) <- hv_combined


    rds_norm=list(data_mrds=data_mrds,data_trds=data_trds,combine_rds=data_crds)

    return(rds_norm)
  })
})

#' @title Cluster Resolution UI
#' @description Dynamically renders a slider input for cluster resolution or number of clusters
#' based on the selected clustering method (e.g., k-means vs. Louvain/SLM).
#' @return A UI element: sliderInput.
output$culster_resolution_button_container <- renderUI({
  if (input$cluster_select=="UMAP-kmeans" || input$cluster_select=="PCA-kmeans") {
    sliderInput("cluster_resolution","Number of Clusters",value = 7,min=2,max=30,step = 1)
  } else {
    sliderInput("cluster_resolution","Resolution of clusters",value =2,min=0.2,max=2,step = 0.1)
  }
})
#' @title Clustering Execution Reactive
#' @description Runs the selected clustering algorithm (Louvain, SLM, K-means, etc.) on the processed data.
#' Triggered by 'Start clustering computation' button.
cluster_plotlist <-eventReactive(c(input$start_cluster), {
  req(rds_norm())
  withProgress(message = "Processing data...",value=0.8,{
    if(!is.null(input$cluster_resolution)){
      rds_norm<-rds_norm()
      source("./source/OverallAnalysisFunction/Clustering/clusterplot.R")
      data_mrds=rds_norm[["data_mrds"]]
      data_trds=rds_norm[["data_trds"]]
      combine_rds=rds_norm[["combine_rds"]]
      
      # Determine clustering type based on selection
      if(input$cluster_select=="UMAP-kmeans"){
        clustertype=0
        k=7 # Note: This k seems unused if resolution slider provides K
      }else if(input$cluster_select=="PCA-kmeans"){
        clustertype=4
        k=7
      }else if(input$cluster_select=="LV"){
        clustertype=1
      }else if(input$cluster_select=="LM"){
        clustertype=2
      }else{
        # SLM
        clustertype=3
      }
      
      # Run clustering visualization/calculation
      resolution=input$cluster_resolution
      mplot<-run_clusterplot(data_mrds,type="Metabolite",pointSize=1,breakseq=50,resolution=resolution,clustertype)
      tplot<-run_clusterplot(data_trds,type="Gene",pointSize=1,breakseq=50,resolution=resolution,clustertype)
      cplot<-run_clusterplot(combine_rds,type="Merge",pointSize=1,breakseq=50,resolution=resolution,clustertype)
      plotlist<-list(mplot,tplot,cplot)
      return(plotlist)
      gc()
    }else{
      return(NULL)
    }
    
  })
})

#' @title Render Cluster Plots
#' @description Triggers the rendering of clustering plots for Metabolite, Gene, and Merged datasets
#' when 'Start clustering computation' is clicked. It arranges them in a grid.
#' @return Updates 'output$cluster_plot' with the combined plot.
observeEvent(c(input$start_cluster),{
  req(cluster_plotlist())
  cluster_plotlist<-cluster_plotlist()
  mplot<-cluster_plotlist[[1]][[1]]
  tplot<-cluster_plotlist[[2]][[1]]
  cplot<-cluster_plotlist[[3]][[1]]
  
  output$cluster_plot <- renderPlot({
    withProgress(message = "Plotting...",value=0.8,{
      p<-gridExtra::grid.arrange(mplot,tplot,cplot,ncol=3)
      p
      assign("cluster_plot", p, envir = .GlobalEnv)
    })
  })
})
#' @title Sankey Data Reactive
#' @description Prepares the data for the Sankey diagram by calculating overlaps and transitions
#' between Metabolite, Gene, and Merged clusters.
#' @return A dataframe suitable for generating a Sankey diagram, containing source, target, value, and color columns.
cluster_sanky <-eventReactive(c(input$start_cluster), {
  req(cluster_plotlist())
  withProgress(message = "Processing data...",value=0.8,{
    cluster_plotlist<-cluster_plotlist()
    mcluster<-cluster_plotlist[[1]][[3]]
    tcluster<-cluster_plotlist[[2]][[3]]
    ccluster<-cluster_plotlist[[3]][[3]]
    source("./source/main_program/slim_seurat.R")
    mcluster <- slim_seurat(seurat_obj=mcluster)
    tcluster <- slim_seurat(seurat_obj=tcluster)
    ccluster <- slim_seurat(seurat_obj=ccluster)
    gc()
    mcluster_data<-data.frame(x_y=mcluster$x_y,clusters_m=paste0("Metabolite_",mcluster$seurat_clusters))
    tcluster_data<-data.frame(x_y=tcluster$x_y,clusters_t=paste0("Gene_",tcluster$seurat_clusters))
    ccluster_data<-data.frame(x_y=ccluster$x_y,clusters_c=paste0("Merge_",ccluster$seurat_clusters))
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
    return(cluster_sanky)
  })
})
#' @title Sankey Plot Generation Reactive
#' @description Generates the Sankey diagram object using the prepared data.
#' @return A networkD3 object representing the Sankey diagram.
cluster_sanky_plot_save <- eventReactive(c(input$start_cluster), {
  req(cluster_sanky())
  withProgress(message = "Processing data...",value=0.8,{
    source("./source/OverallAnalysisFunction/Sanky/sankyplot.R")
    links<-cluster_sanky()
    p=sankyplot_ID(links)
    
    return(p)
  })
})
#' @title Render Sankey Plot
#' @description Renders the Sankey diagram in the UI.
#' @return Updates 'output$cluster_sanky_plot' with the interactive plot.
observeEvent(input$start_cluster, {
  output$cluster_sanky_plot <- renderUI({
    req( cluster_sanky_plot_save())
    withProgress(message = "Plotting...", value = 0.8, {
      cluster_sanky_plot_save()
    })
  })
})

#' @title PCA Assessment Reactive
#' @description Computes Principal Component Analysis (PCA) for data assessment.
#' Generates Elbow Plots (variance explained) and Spatial PC Plots.
#' Triggered by 'Run Data Assessment' button.
pca_analysis_reactive <- eventReactive(input$run_assessment, {
  req(rds_norm())
  withProgress(message = "Calculating PCA...", value = 0.5, {
    rds_list <- rds_norm()
    
    # Run PCA on each dataset (Metab, Trans, Merge)
    m_obj <- RunPCA(rds_list$data_mrds, verbose = FALSE, features = VariableFeatures(rds_list$data_mrds), npcs = 30)
    t_obj <- RunPCA(rds_list$data_trds, verbose = FALSE, features = VariableFeatures(rds_list$data_trds), npcs = 30)
    c_obj <- RunPCA(rds_list$combine_rds, verbose = FALSE, features = VariableFeatures(rds_list$combine_rds), npcs = 30)
    
    # Function to generate Elbow Plot (Standard Deviation vs PC)
    get_elbow_plot <- function(obj, title) {
      stdev <- obj@reductions$pca@stdev
      
      # Plot all computed PCs (30)
      ndims_to_plot <- length(stdev) 
      data <- data.frame(PC = 1:ndims_to_plot, Stdev = stdev)
      
      p <- ggplot(data, aes(x = PC, y = Stdev)) +
        geom_line(color = "black") +
        geom_point(color = "black") +
        theme(panel.grid=element_blank())+
        theme_classic()+
        theme(axis.line.x = element_line(color="black", size = 0.5),
              axis.line.y = element_line(color="black", size = 0.5),
              plot.title = element_text(hjust = 0.5))+
        labs(title = title, y = "Standard Deviation", x = "PC") +
        theme(plot.title = element_text(hjust = 0.5))
      
      return(p)
    }

    elbow_m <- get_elbow_plot(m_obj, "Metabolite Variance")
    elbow_t <- get_elbow_plot(t_obj, "Gene Variance")
    elbow_c <- get_elbow_plot(c_obj, "Merge Variance")
    
    # Function to generate Spatial Feature Plot for top PCs
    get_spatial_pc_plot <- function(obj, title_prefix) {
      # Extract top 3 PCs
      embeddings <- Embeddings(obj, reduction = "pca")[, 1:3]
      colnames(embeddings) <- c("PC_1", "PC_2", "PC_3")
      plot_data <- cbind(obj@meta.data, embeddings)
      
      plot_list <- list()
      for(i in 1:3) {
        pc_name <- paste0("PC_", i)
        p <- ggplot(plot_data, aes(x = x, y = y, color = .data[[pc_name]])) +
          geom_point(shape = 19, size = 0.5) +
          scale_color_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)) +
          coord_equal() +
          theme_void() +
          ggtitle(paste0(title_prefix, " ", pc_name)) +
          theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "right") 
        plot_list[[i]] <- p
      }
      return(plot_list)
    }

    spatial_m <- get_spatial_pc_plot(m_obj, "Metabolite")
    spatial_t <- get_spatial_pc_plot(t_obj, "Gene")
    spatial_c <- get_spatial_pc_plot(c_obj, "Merge")
    
    list(
        elbow = list(m=elbow_m, t=elbow_t, c=elbow_c),
        spatial = list(m=spatial_m, t=spatial_t, c=spatial_c)
    )
  })
})

#' @title Render PCA Elbow Plot
#' @description Renders the Elbow plots (Standard Deviation vs PC) for all three datasets.
#' @return Updates 'output$pca_elbow_plot'.
output$pca_elbow_plot <- renderPlot({
  req(pca_analysis_reactive())
  plots <- pca_analysis_reactive()$elbow
  p <- gridExtra::grid.arrange(plots$m, plots$t, plots$c, ncol = 3)
  assign("pca_elbow_plot_obj", p, envir = .GlobalEnv)
  p
})

#' @title Render Spatial PCA Plot
#' @description Renders the spatial distribution of the top 3 Principal Components for all three datasets.
#' @return Updates 'output$pca_spatial_plot'.
output$pca_spatial_plot <- renderPlot({
  req(pca_analysis_reactive())
  withProgress(message = "Rendering Spatial PCA Plots...", value = 0.5, {
      plots <- pca_analysis_reactive()$spatial
      # Arrange 3x3: Rows = Dataset (Metab, Gene, Merge), Cols = PC1, PC2, PC3
      p <- gridExtra::grid.arrange(
        plots$m[[1]], plots$m[[2]], plots$m[[3]],
        plots$t[[1]], plots$t[[2]], plots$t[[3]],
        plots$c[[1]], plots$c[[2]], plots$c[[3]],
        ncol = 3, nrow = 3
      )
      assign("pca_spatial_plot_obj", p, envir = .GlobalEnv)
      p
  })
})

#' @title Download PCA Elbow Plot
#' @description Saves the PCA Elbow plot as a PNG file.
#' @return A PNG file named "pca_elbow_plot.png".
output$download_pca_elbow_plot <- downloadHandler(
  filename = function() {
    "pca_elbow_plot.png"
  },
  content = function(file) {
    req(exists("pca_elbow_plot_obj", envir = .GlobalEnv))
    p <- get("pca_elbow_plot_obj", envir = .GlobalEnv)
    ggsave(file, plot = p, device = "png", width = 24, height = 8, bg = "#FFFFFF", dpi = 300)
  }
)

#' @title Download Spatial PCA Plot
#' @description Saves the Spatial PCA plot as a PNG file.
#' @return A PNG file named "pca_spatial_plot.png".
output$download_pca_spatial_plot <- downloadHandler(
  filename = function() {
    "pca_spatial_plot.png"
  },
  content = function(file) {
    req(exists("pca_spatial_plot_obj", envir = .GlobalEnv))
    p <- get("pca_spatial_plot_obj", envir = .GlobalEnv)
    ggsave(file, plot = p, device = "png", width = 20, height = 20, bg = "#FFFFFF", dpi = 300)
  }
)
