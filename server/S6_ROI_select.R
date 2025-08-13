
output$download_samplelist <- downloadHandler(
  filename = function() {

    paste0("samplelist.txt")
  },
  content = function(file) {
    DemoData<-samplelist()
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
  }
)
output$download_Segmentation_plot <- downloadHandler(
  filename = function() {
    "Segmentation_plot.png"
  },
  content = function(file) {

    p <- get("Segmentation_plot", envir = .GlobalEnv)

    ggsave(file, plot = p, device = "png", width = 8, height = 8,bg = "#FFFFFF", dpi = 300)

  })
output$download_group_spectrum_reactive <- downloadHandler(
  filename = function() {
    "Segmentation_group_plot.png"
  },
  content = function(file) {
    p <- group_spectrum_reactive()
    ggsave(file, plot = p, device = "png", width = 12, height = 8,bg = "#FFFFFF", dpi = 300)
    
  })

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


observeEvent(input$auto_add_to_control, {
  req(peak())
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


output$scatter_plot <- renderPlotly({
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

scatter_plot_cluster_save<- eventReactive(c(input$start_cell_annotation,input$plot_select,
                                            input$clusterdata_select), {
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
output$scatter_plot_cluster <- renderPlot({
  scatter_plot_cluster_save()
})

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

control <- reactiveVal(list())  
treatment <- reactiveVal(list())  

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

selected_points <- reactiveVal()

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

plotdata_group <- eventReactive(input$finish_selection, {
  req(plotdata())
  withProgress(message = "Processing data...",value=0.8,{
  data <- plotdata()
  if (input$plot_select == "Use selection tool"){
  split_groupunique<-split_groupunique()
  }else{
    split_groupunique<-split_groupunique_cluster()
  }
  for(i in names(split_groupunique)){
   k <- which(row.names(data) %in% row.names(split_groupunique[[i]]))
   data$groups[k] <- i
  }
  
  return(data)
  })
})

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

output$group_spectrum <- renderPlot({
  req(group_spectrum_reactive())
  p<-group_spectrum_reactive()
  p

})

observeEvent(input$clear_selection, {
  rv$clear <- TRUE
})

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
  


observeEvent(input$plot_select, {
  if(input$plot_select == "Select cell groups"){
    shinyjs::hide(id = "clusterdata_select")
  }else{
    shinyjs::show(id = "clusterdata_select")
  }
})

output$ion_select_button_container <- renderUI({
  if (input$plot_select=="Use selection tool") {
    selectizeInput("ion_select", "Select a feature to plot spatial distributions:", choices = NULL,options = list(server = TRUE),selected = NULL)
  } else{
    NULL
  }
})

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
  
  output$scatter_plot_button_container <- renderUI({
    
    if (input$plot_select=="Use selection tool" ){
      plotlyOutput("scatter_plot",width = "100%", height = "600px")
    }else{
      plotOutput("scatter_plot_cluster",width = "100%", height = "600px")
    }
    
  })
