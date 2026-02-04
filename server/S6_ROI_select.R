# ==============================================================================
# S6_ROI_select.R
# Server logic for Step 4: ROI Selection
# Handles definition of Regions of Interest (ROI) via clustering, cell types, or interactive selection.
# ==============================================================================

# ------------------------------------------------------------------------------
# UI Renderers
# ------------------------------------------------------------------------------

output$roi_prerequisite_warning <- renderUI({
  # Check if upstream analysis is complete
  # cell_annotation_rds is the output of the previous step
  # We use tryCatch or simply check for NULL/error state
  is_ready <- tryCatch({
    !is.null(cell_annotation_rds())
  }, error = function(e) FALSE)
  
  if (!is_ready) {
    div(style = "color: red; font-weight: bold; margin-bottom: 10px;", 
        "⚠️ Prerequisite: Please complete 'Clustering Analysis' and 'Cell Annotation' before starting ROI Selection.")
  } else {
    NULL
  }
})

# ------------------------------------------------------------------------------
# Download Handlers
# ------------------------------------------------------------------------------

#' @title Download ROI Assignments
output$download_samplelist <- downloadHandler(
  filename = function() {
    
    paste0("samplelist.txt")
  },
  content = function(file) {
    DemoData<-samplelist()
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
  }
)
#' @title Download Segmentation Plot
output$download_Segmentation_plot <- downloadHandler(
  filename = function() {
    "Segmentation_plot.png"
  },
  content = function(file) {
    
    p <- get("Segmentation_plot", envir = .GlobalEnv)
    
    ggsave(file, plot = p, device = "png", width = 8, height = 8,bg = "#FFFFFF", dpi = 300)
    
  })
#' @title Download Group Spectrum Plot
output$download_group_spectrum_reactive <- downloadHandler(
  filename = function() {
    "Segmentation_group_plot.png"
  },
  content = function(file) {
    p <- group_spectrum_reactive()
    ggsave(file, plot = p, device = "png", width = 12, height = 8,bg = "#FFFFFF", dpi = 300)
    
  })
# ------------------------------------------------------------------------------
# Core ROI Logic
# ------------------------------------------------------------------------------

#' @title Prepare Data for ROI Selection
#' @description Extracts relevant Seurat object (Metabolite, Gene, or Merge) based on user selection.
#' Ensures metadata columns (clusters, celltype) are character type.
peak <- eventReactive(c(input$start_cell_annotation,input$clusterdata_select,input$start_cluster), {
  cell_annotation_rds<-cell_annotation_rds()
  req(cell_annotation_rds())
  if(!is.null(input$clusterdata_select) && input$clusterdata_select=="Metabolite"){
    data <-cell_annotation_rds[[1]]
  }else if(!is.null(input$clusterdata_select) && input$clusterdata_select=="Gene"){
    data <-cell_annotation_rds[[2]]
  }else{
    data <-cell_annotation_rds[[3]]
  }
  
  withProgress(message = "Processing data...",value=0.8,{
    
    data@meta.data$seurat_clusters<-as.character(data@meta.data$seurat_clusters)
    data@meta.data$celltype<-as.character(data@meta.data$celltype)

    
    return(data)
    
  })
})

#' @title Update Feature Selection Dropdown
#' @description Populates the "Select a feature" dropdown with available ions/features 
#' when "Use selection tool" mode is active.
observeEvent(c(input$plot_select,peak()), {
  req(peak())
  if (input$plot_select == "Use selection tool") {
    data <- peak()
    ions <- rownames(data@assays$Spatial$counts)
    ions<-sort(as.character(ions))
    ions <-c("total counts",ions)
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      inputId = "ion_select",
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
  }
})

#' @title Update Group Selection Choices
#' @description Updates the available choices for Treatment/Control groups based on
#' the selected annotation mode (Cluster or Cell Type).
observeEvent(c(peak(),input$plot_select,input$clusterdata_select,input$clear_selection), {
  req(peak())
  data <- peak()
  if(input$plot_select !="Use selection tool"){
    
    if(input$plot_select=="Select clustering groups"){
      clusters <- unique(data@meta.data$seurat_clusters)
      clusters<-sort(as.character(clusters))
    }else if(input$plot_select=="Select cell groups"){
      clusters <- unique(data$celltype)
      clusters<-sort(as.character(clusters))
    }
    
    
    if (length(clusters) > 0) {
      updateSelectizeInput(session = getDefaultReactiveDomain(), "auto_add_to_treatment", 
                           choices = clusters, selected = NULL)#clusters[1]
      
      if (length(clusters) > 1) {
        updateSelectizeInput(session = getDefaultReactiveDomain(), "auto_add_to_control",
                             choices = clusters, selected = NULL)
        
      }
    }
  }
  
})
#' @title Filter Control Group Choices
#' @description Dynamically removes selected Treatment groups from the Control group options
#' to prevent overlapping selections.
observeEvent(input$auto_add_to_treatment, {
  req(peak())
  data <- peak()
  if(input$plot_select!="Use selection tool"){
    
    if(input$plot_select=="Select clustering groups"){
      clusters <- unique(data@meta.data$seurat_clusters)
      clusters<-sort(as.character(clusters))
    }else if(input$plot_select=="Select cell groups"){
      clusters <- unique(data@meta.data$celltype)
      clusters<-sort(as.character(clusters))
    }
    
    updateSelectInput(session=getDefaultReactiveDomain(), "auto_add_to_control",
                      choices = setdiff(clusters, input$auto_add_to_treatment),
                      selected = input$auto_add_to_control)
  }
  
})


#' @title Update Feature Selection Dropdown
#' @description Populates the "Select a feature" dropdown with available ions/features 
#' when "Use selection tool" mode is active.
observeEvent(c(input$plot_select,peak()), {
  req(peak())
  if (input$plot_select == "Use selection tool") {
    data <- peak()
    ions <- rownames(data@assays$Spatial$counts)
    ions<-sort(as.character(ions))
    ions <-c("total counts",ions)
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      inputId = "ion_select",
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
  }
})

#' @title Update Group Selection Choices
#' @description Updates the available choices for Treatment/Control groups based on
#' the selected annotation mode (Cluster or Cell Type).
observeEvent(c(peak(),input$plot_select,input$clusterdata_select,input$clear_selection), {
  req(peak())
  data <- peak()
  if(input$plot_select !="Use selection tool"){
    
    if(input$plot_select=="Select clustering groups"){
      clusters <- unique(data@meta.data$seurat_clusters)
      clusters<-sort(as.character(clusters))
    }else if(input$plot_select=="Select cell groups"){
      clusters <- unique(data$celltype)
      clusters<-sort(as.character(clusters))
    }
    
    
    if (length(clusters) > 0) {
      updateSelectizeInput(session = getDefaultReactiveDomain(), "auto_add_to_treatment", 
                           choices = clusters, selected = NULL)#clusters[1]
      
      if (length(clusters) > 1) {
        updateSelectizeInput(session = getDefaultReactiveDomain(), "auto_add_to_control",
                             choices = clusters, selected = NULL)
        
      }
    }
  }
  
})
#' @title Filter Control Group Choices
#' @description Dynamically removes selected Treatment groups from the Control group options
#' to prevent overlapping selections.
observeEvent(input$auto_add_to_treatment, {
  req(peak())
  data <- peak()
  if(input$plot_select!="Use selection tool"){
    
    if(input$plot_select=="Select clustering groups"){
      clusters <- unique(data@meta.data$seurat_clusters)
      clusters<-sort(as.character(clusters))
    }else if(input$plot_select=="Select cell groups"){
      clusters <- unique(data@meta.data$celltype)
      clusters<-sort(as.character(clusters))
    }
    
    updateSelectInput(session=getDefaultReactiveDomain(), "auto_add_to_control",
                      choices = setdiff(clusters, input$auto_add_to_treatment),
                      selected = input$auto_add_to_control)
  }
  
})


#' @title Filter Treatment Group Choices
#' @description Dynamically removes selected Control groups from the Treatment group options
#' to prevent overlapping selections.
observeEvent(input$auto_add_to_control, {
  data <- peak()
  if(input$plot_select !="Use selection tool"){
    if(input$plot_select=="Select clustering groups"){
      clusters <- unique(data@meta.data$seurat_clusters)
      clusters<-sort(as.character(clusters))
    }else if(input$plot_select=="Select cell groups"){
      clusters <- unique(data@meta.data$celltype)
      clusters<-sort(as.character(clusters))
    }
    updateSelectInput(session=getDefaultReactiveDomain(), "auto_add_to_treatment",
                      choices = setdiff(clusters, input$auto_add_to_control),
                      selected = input$auto_add_to_treatment)
  }
  
})

#' @title Prepare Plot Data for Interactive Selection
#' @description Extracts coordinate and intensity data for visualization.
#' Supports visualizing Total Counts or specific features (Genes/Metabolites).
plotdata <-eventReactive(c(input$ion_select,input$start_cell_annotation,
                           input$clusterdata_select), {
                             req(peak())
                             
                             withProgress(message = "Processing plotdata",value=0.8,{
                               
                               
                               print("plotdata start")
                               data <- peak()
                               plotdata <- data@meta.data
                               
                               
                               
                               
                               if(!is.null(input$ion_select) && input$ion_select!="" && input$ion_select!="total counts"){
                                 ion <- as.character(input$ion_select)
                                 
                                 plotdata$intensity<- data@assays$Spatial$counts[ion,]
                                 plotdata$norm_intensity <-100*(plotdata$intensity)/max(plotdata$intensity)
                               }else{
                                 plotdata$intensity<- plotdata$nCount_Spatial
                                 plotdata$norm_intensity <-100*(plotdata$intensity)/max(plotdata$intensity)
                               }
                               
                               return(plotdata)
                               print("plotdata done")
                               
                             })
                             
                           })  
#' @title Generate Color Mapping
#' @description Creates a consistent color palette for clusters or cell types.
color_mapping <- eventReactive(c(input$start_cell_annotation,input$plot_select,
                                 input$clusterdata_select), {
                                   req(plotdata())
                                   
                                   withProgress(message = "Processing data...",value=0.8,{
                                     
                                     data <- plotdata()
                                     qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
                                     cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
                                     if(input$plot_select=="Select cell groups"){
                                       k=length(unique(data$celltype))
                                       color_mapping <- setNames(cluster_Palette[1:k], unique(data$celltype))
                                       
                                     }else{
                                       k=length(unique(data$seurat_clusters))
                                       
                                       color_mapping <- setNames(cluster_Palette[1:k], unique(data$seurat_clusters))
                                     }
                                     
                                     return(color_mapping)
                                     
                                   })
                                   
                                 })
#' @title Render Interactive Scatter Plot
#' @description Generates an interactive Plotly scatter plot for manual ROI selection.
#' Allows users to select regions using the lasso tool. Points are colored by intensity.
output$scatter_plot <- renderPlotly({
  validate(
    need(try(peak(), silent = TRUE), "Please complete 'Clustering Analysis' and 'Cell Annotation' first. Data is missing.")
  )
  req(plotdata())
  data <- plotdata()
  
  plot_ly(data, x = ~x, y = ~y, color = ~norm_intensity, colors =colorRampPalette(c("blue", "red"))(100),
          type = 'scatter', mode = 'markers', source = "scatter",marker = list(size = 5)) %>%
    layout(
      xaxis = list(
        range = c(min(data$x),max(data$x)),
        constrain = "domain"
      ),
      yaxis = list(
        range = c(min(data$y),max(data$y)),
        scaleanchor = "x",
        scaleratio = 1
      ),
      dragmode = "lasso"
    )
  
  
})

#' @title Generate Static Cluster/Cell Plot
#' @description Creates a static ggplot of the spatial data, colored by Cluster ID or Cell Type.
#' Used when manual lasso selection is NOT active.
scatter_plot_cluster_save<- eventReactive(c(input$start_cell_annotation,input$plot_select,
                                            input$clusterdata_select), {
                                              validate(
                                                need(try(peak(), silent = TRUE), "Please complete 'Clustering Analysis' and 'Cell Annotation' first.")
                                              )
                                              req(c(plotdata(),color_mapping()))
                                              data <- plotdata()
                                              color_mapping <- color_mapping()
                                              
                                              if(input$plot_select=="Select clustering groups"){
                                                p<-ggplot(data, aes(x = x, y = y)) +
                                                  geom_point(aes(color = seurat_clusters), size = 1) +
                                                  scale_color_manual(values = color_mapping) +
                                                  guides(colour = guide_legend(title = "Group",override.aes = list(size=3), nrow = 10))+
                                                  theme_minimal() +
                                                  xlim(min(data$x), max(data$x)) +
                                                  ylim(min(data$y), max(data$y)) +
                                                  coord_equal() 
                                              }else if(input$plot_select=="Select cell groups"){
                                                if(length(unique(data$celltype))>30){
                                                  pointSize=0.5
                                                }else{
                                                  pointSize=1
                                                }
                                                p<-ggplot(data, aes(x = x, y = y)) +
                                                  geom_point(aes(color = celltype), size = pointSize) +
                                                  scale_color_manual(values = color_mapping) +
                                                  guides(colour = guide_legend(title = "Group",override.aes = list(size=3), nrow = 10))+
                                                  theme_minimal() +
                                                  xlim(min(data$x), max(data$x)) +
                                                  ylim(min(data$y), max(data$y)) +
                                                  coord_equal() 
                                              }
                                              
                                              
                                              return(p)
                                              
                                            })
#' @title Render Static Cluster Plot
#' @description Displays the static cluster/cell type visualization.
output$scatter_plot_cluster <- renderPlot({
  scatter_plot_cluster_save()
})
#' @title Define Control Group (Cluster/Cell Type)
#' @description Filters data based on user-selected clusters/cell types for the Control group.
control_cluster<-eventReactive(c(input$auto_add_to_control), {
  req(peak())
  data <- peak()
  if (input$plot_select != "Use selection tool") {
    if(input$plot_select=="Select clustering groups"){
      control_cluster<-data@meta.data %>% filter(seurat_clusters %in% input$auto_add_to_control)
    }else if(input$plot_select=="Select cell groups"){
      control_cluster<-data@meta.data %>% filter(celltype %in% input$auto_add_to_control)
    }
    return(control_cluster)
  }
})
#' @title Define Treatment Group (Cluster/Cell Type)
#' @description Filters data based on user-selected clusters/cell types for the Treatment group.
treatment_cluster<-eventReactive(c(input$auto_add_to_treatment), {
  req(peak())
  data <- peak()
  if (input$plot_select != "Use selection tool") {
    if(input$plot_select=="Select clustering groups"){
      treatment_cluster<-data@meta.data %>% filter(seurat_clusters %in% input$auto_add_to_treatment)
    }else if(input$plot_select=="Select cell groups"){
      treatment_cluster<-data@meta.data %>% filter(celltype %in% input$auto_add_to_treatment)
    }
    return(treatment_cluster)
  }
})
#' @title Cluster Selection Table Data
#' @description Reactive value to store metadata of selected groups in Cluster/Cell Type mode.
selection_data_c <- reactiveVal(data.frame(
  Group = character(),
  Region = character(),
  X_Min = numeric(),
  X_Max = numeric(),
  Y_Min = numeric(),
  Y_Max = numeric(),
  Point_Count = integer(),
  stringsAsFactors = FALSE
))
#' @title Update Cluster Selection Table
#' @description Adds or updates a row in the selection table when a cluster/cell type is added to a group.
update_selection_data_cluster <- function(group_name, selected_df,cluster_annotation_select) {
  if(cluster_annotation_select=="Select cell groups"){
    region_name <- paste0("Cluster_", paste(unique(selected_df$celltype), collapse = "_"))
    
  }else{
    
    region_name <- paste0("Cluster_", paste(unique(selected_df$seurat_clusters), collapse = "_"))
  }
  
  x_min <- min(selected_df$x)
  x_max <- max(selected_df$x)
  y_min <- min(selected_df$y)
  y_max <- max(selected_df$y)
  
  point_count <- nrow(selected_df)
  
  current_selection_data <- selection_data_c()
  new_row <- data.frame(
    Group = group_name,
    Region = region_name,
    X_Min = x_min,
    X_Max = x_max,
    Y_Min = y_min,
    Y_Max = y_max,
    Point_Count = point_count,
    stringsAsFactors = FALSE
  )
  
  existing_row_index <- which(current_selection_data$Group == new_row$Group[1])
  
  if (length(existing_row_index) > 0) {
    current_selection_data[existing_row_index, ] <- new_row
    selection_data_c(current_selection_data)
  } else{
    selection_data_c(rbind(current_selection_data, new_row))
  }
  
  output$selection_table <- renderTable({
    selection_data_c()
  })
}

#' @title Manual Selection Storage
#' @description Reactive lists to store point coordinates selected via the lasso tool.
control <- reactiveVal(list())  
treatment <- reactiveVal(list())  

#' @title Manual Selection Table Data
#' @description Reactive value to store metadata of manual lasso selections.
selection_data <- reactiveVal(data.frame(
  Group = character(),
  Region = character(),
  X_Min = numeric(),
  X_Max = numeric(),
  Y_Min = numeric(),
  Y_Max = numeric(),
  Point_Count = integer(),
  stringsAsFactors = FALSE
))

#' @title Current Lasso Selection
#' @description Stores the points currently highlighted by the user's lasso tool in Plotly.
selected_points <- reactiveVal()

#' @title Capture Plotly Selection Event
#' @description Listens for 'plotly_selected' events and extracts the corresponding data points.
observe({
  data <- plotdata()
  selected_data <- event_data("plotly_selected", source = "scatter")
  
  if (is.null(selected_data)) {
    selected_points(NULL)
  } else {
    selected_points(data[selected_data$pointNumber + 1, ])
  }
})

region_count <- reactiveVal(1)

#' @title Update Manual Selection Table
#' @description Adds a new row to the selection table representing a manually defined region.
update_selection_data <- function(group_name, selected_df) {
  region_name <- paste0("Region ", region_count())
  region_count(region_count() + 1)
  
  x_min <- min(selected_df$x)
  x_max <- max(selected_df$x)
  y_min <- min(selected_df$y)
  y_max <- max(selected_df$y)
  
  point_count <- nrow(selected_df)
  
  current_selection_data <- selection_data()
  new_row <- data.frame(
    Group = group_name,
    Region = region_name,
    X_Min = x_min,
    X_Max = x_max,
    Y_Min = y_min,
    Y_Max = y_max,
    Point_Count = point_count,
    stringsAsFactors = FALSE
  )
  
  selection_data(rbind(current_selection_data, new_row))
  output$selection_table <- renderTable({
    selection_data()
  })
}

#' @title Add Selection to Control Group
#' @description Moves the currently selected points (from lasso) into the Control group list.
observeEvent(input$add_to_control, {
  if (input$plot_select == "Use selection tool") {
    if (!is.null(selected_points())) {
      
      new_control <- control()
      new_control <- append(new_control, list(selected_points()))
      control(new_control) 
      update_selection_data("control", selected_points())
      showNotification("The selection has been added to the control group.")
    } else {
      showNotification("No regions have been selected.")
    }
  }
  
})

#' @title Add Selection to Treatment Group
#' @description Moves the currently selected points (from lasso) into the Treatment group list.
observeEvent(input$add_to_treatment, {
  if (input$plot_select == "Use selection tool") {
    if (!is.null(selected_points())) {
      new_treatment <- treatment()
      new_treatment <- append(new_treatment, list(selected_points()))
      treatment(new_treatment)  
      update_selection_data("treatment", selected_points())
      showNotification("The selection has been added to the treatment group.")
    } else {
      showNotification("No regions have been selected.")
    }
  }
})

#' @title Validate Final Selection
#' @description Checks if both Treatment and Control groups have been populated before proceeding.
observeEvent(input$finish_selection, {
  treatment <- treatment()
  control <- control()
  treatment_cluster <- treatment_cluster()
  control_cluster <- control_cluster()
  if (input$plot_select == "Use selection tool"){
    validate(
      need(length(treatment)>0,"Please select an treatment area and a control area.")
    )
    validate(
      need(length(control)>0,"Please select an treatment area and a control area.")
    )
  }else{
    validate(
      need(length(treatment_cluster)>0,"Please select an treatment area and a control area.")
    )
    validate(
      need(length(control_cluster)>0,"Please select an treatment area and a control area.")
    )
  }
})

#' @title Process Manual Groups
#' @description Combines and labels the manually selected points for Treatment and Control.
#' Handles duplicates and returns a split list.
split_groupunique <- eventReactive(input$finish_selection, {
  
  treatment <- treatment()
  control <- control()
  
  control_all <- bind_rows(control) %>%
    mutate(groups = "control") %>%
    distinct(cell, .keep_all = TRUE)
  
  treatment_all <- bind_rows(treatment) %>%
    mutate(groups = "treatment") %>%
    distinct(cell, .keep_all = TRUE)
  
  groupunique <- rbind(treatment_all, control_all) %>%
    distinct(cell, .keep_all = TRUE) 
  
  split_groupunique <- split(groupunique, groupunique$groups)
  return(split_groupunique)
  
})
#' @title Process Cluster Groups
#' @description Combines and labels the selected Clusters/Cell Types for Treatment and Control.
#' Returns a split list.
split_groupunique_cluster <- eventReactive(input$finish_selection, {
  
  treatment_cluster <- treatment_cluster()
  control_cluster <- control_cluster()
  str(control_cluster)
  control_all <- control_cluster %>%
    mutate(groups = "control") %>%
    distinct(cell, .keep_all = TRUE)
  
  treatment_all <- treatment_cluster %>%
    mutate(groups = "treatment") %>%
    distinct(cell, .keep_all = TRUE)
  
  groupunique <- rbind(treatment_all, control_all) %>%
    distinct(cell, .keep_all = TRUE)
  
  split_groupunique_cluster <- split(groupunique, groupunique$groups)
  
  return(split_groupunique_cluster)
  
  
})

#' @title Annotate Metadata with Groups
#' @description Adds a "groups" column (Treatment/Control/0) to the main metadata frame.
plotdata_group <- eventReactive(input$finish_selection, {
  req(plotdata())
  withProgress(message = "Processing data...",value=0.8,{
    data <- plotdata()
    if (input$plot_select == "Use selection tool"){
      split_groupunique<-split_groupunique()
    }else{
      split_groupunique<-split_groupunique_cluster()
    }
    
    if(!"groups" %in% colnames(data)) {
      data$groups <- "0"
    }
    for(i in names(split_groupunique)){
      k <- which(row.names(data) %in% row.names(split_groupunique[[i]]))
      data$groups[k] <- i
    }
    
    return(data)
  })
})

#' @title Update Seurat Objects with Groups
#' @description Adds the "groups" metadata column to both Metabolomics and Transcriptomics Seurat objects.
#' @return A list containing the updated Seurat objects.
data_rds_group <- eventReactive(input$finish_selection, {
  cluster_rds <- cluster_plotlist()
  req(cluster_rds)
  data_mrds<-cluster_rds[[1]][[3]]
  data_trds <- cluster_rds[[2]][[3]]
  req(data_mrds,data_trds)
  withProgress(message = "Processing data...",value=0.8,{
    
    if (input$plot_select == "Use selection tool"){
      split_groupunique<-split_groupunique()
    }else{
      split_groupunique<-split_groupunique_cluster()
    }
    change_rds_group<-function(data,split_groupunique){
      if(!"groups" %in% colnames(data@meta.data)) {
        data@meta.data$groups <- "0"
      }
      for(i in names(split_groupunique)){
        k <- which(row.names(data@meta.data) %in% row.names(split_groupunique[[i]]))
        data@meta.data$groups[k] <- i
      }
      return(data)
    }
    data_mrds_group=change_rds_group(data_mrds,split_groupunique)
    data_trds_group=change_rds_group(data_trds,split_groupunique)
    
    data_rds_group<-list(data_mrds=data_mrds_group,data_trds=data_trds_group)
    
    return(data_rds_group)
  })
})

#' @title Create Final Sample List
#' @description Generates a mapping dataframe of spots and their assigned group (Treatment/Control).
#' Used for downstream differential analysis.
samplelist <- eventReactive(input$finish_selection, {
  req(plotdata_group())
  withProgress(message = "Processing data...",value=0.8,{
    
    plotdata_group <- plotdata_group()
    print(unique(plotdata_group$groups))
    samplelist<- plotdata_group %>%
      mutate(gene = paste0("x",plotdata_group$x,"_y", plotdata_group$y)) 
    
    samplelist<-samplelist %>%
      mutate(metabolomics = samplelist$gene) %>%
      mutate(group = samplelist$groups) %>%
      dplyr::select(gene,metabolomics,group) 
    
    return(samplelist)
  })
})

rv <- reactiveValues(clear = FALSE)

#' @title Generate Group Visualization
#' @description Creates a spatial plot showing the final assigned groups (Red=Treatment, Blue=Control).
group_spectrum_reactive <- eventReactive(c(input$finish_selection,input$clear_selection), {
  req(plotdata_group())
  if (rv$clear) {
    rv$clear <- FALSE 
    return(NULL)
  }
  withProgress(message = "Processing data...",value=0.8,{
    data <- plotdata_group()
    
    data <- data %>%
      group_by(groups) 
    group_color<-c("treatment"="red","control"="blue","0"="grey")
    p <-  ggplot(data, aes(x = x, y = y)) +
      geom_point(aes(color = groups), size = 1) +
      guides(colour = guide_legend(title = "Group",override.aes = list(size=3), nrow = 10))+
      scale_color_manual(values = group_color) +
      theme_minimal()+
      xlim(min(data$x),max(data$x))+
      ylim(min(data$y),max(data$y))+
      coord_equal()
    
    return(p)
  })
})    

#' @title Render Group Plot
#' @description Displays the final group assignment plot.
output$group_spectrum <- renderPlot({
  req(group_spectrum_reactive())
  p<-group_spectrum_reactive()
  p
  
})






#' @title Add Cluster to Control Group
#' @description Adds the selected clusters/cell types to the Control group.
#' Updates the selection table with the new group information.
observeEvent(input$auto_add_to_control, {
  if (input$plot_select != "Use selection tool") {
    if (!is.null(control_cluster())) {
      update_selection_data_cluster("control", control_cluster()
                                    ,cluster_annotation_select=input$plot_select)
      showNotification("The selection has been added to the control group.")
    } else {
      showNotification("No clusters have been selected.")
    }
  }
})
#' @title Add Cluster to Treatment Group
#' @description Adds the selected clusters/cell types to the Treatment group.
#' Updates the selection table with the new group information.
observeEvent(input$auto_add_to_treatment, {
  if (input$plot_select != "Use selection tool") {
    if (!is.null(treatment_cluster())) {
      update_selection_data_cluster("treatment", treatment_cluster(),
                                    cluster_annotation_select=input$plot_select)
      
      showNotification("The selection has been added to the treatment group.")
    } else {
      showNotification("No clusters have been selected.")
    }
  }
})


#' @title Trigger Clear Selection Flag
#' @description Sets a reactive flag to indicate that the selection needs to be cleared.
observeEvent(input$clear_selection, {
  rv$clear <- TRUE 
})
#' @title Reset Selection State
#' @description Clears all selection data (treatment/control lists, tables) and resets the UI.
observeEvent(input$clear_selection, {
  
  treatment(list())  
  control(list())  
  
  selection_data(data.frame(
    Group = character(),
    Region = character(),
    X_Min = numeric(),
    X_Max = numeric(),
    Y_Min = numeric(),
    Y_Max = numeric(),
    Point_Count = integer(),
    stringsAsFactors = FALSE
  ))
  
  region_count(1)
  
  output$selection_table <- renderTable({
    selection_data()
  })
  
  showNotification("All selections have been cleared.")
})



#' @title Toggle Cluster Selection UI
#' @description Hides or shows the cluster data selection input based on the chosen plot mode.
observeEvent(input$plot_select, {
  if(input$plot_select == "Select cell groups"){
    shinyjs::hide(id = "clusterdata_select")
  }else{
    shinyjs::show(id = "clusterdata_select")
  }
})

#' @title Dynamic Feature Selection UI
#' @description Renders the feature selection dropdown only when "Use selection tool" is active.
output$ion_select_button_container <- renderUI({
  if (input$plot_select=="Use selection tool") {
    selectizeInput("ion_select", "Select a feature to plot spatial distributions:", choices = NULL,options = list(server = TRUE),selected = NULL)
  } else{
    NULL
  }
})

#' @title Dynamic Action Buttons UI
#' @description Renders the appropriate action buttons (Add, Clear, Finish) based on the current selection mode.
output$auto_add_button_container <- renderUI({
  if (input$plot_select=="Use selection tool"){
    tagList(
      actionButton("add_to_treatment", "Add as treatment group"),
      actionButton("add_to_control", "Add as control group"),
      actionButton("clear_selection", "Clear selection"),  
      actionButton("finish_selection", "Finish selection")
    )
  }else{
    tagList(
      selectizeInput("auto_add_to_treatment", "Select clusters as treatment group:",
                     choices =NULL,
                     selected = NULL,
                     multiple = TRUE),
      selectizeInput("auto_add_to_control", "Select clusters as control group:",
                     choices =NULL,
                     selected = NULL,
                     multiple = TRUE),
      actionButton("clear_selection", "Clear selection"), 
      actionButton("finish_selection", "Finish selection")
      
    )
  }
  
})

#' @title Dynamic Plot Container UI
#' @description Switches between the interactive Plotly widget and the static plot output based on the selection mode.
output$scatter_plot_button_container <- renderUI({
  
  if (input$plot_select=="Use selection tool" ){
    plotlyOutput("scatter_plot",width = "100%", height = "600px")
  }else{
    plotOutput("scatter_plot_cluster",width = "100%", height = "600px")
  }
  
})
