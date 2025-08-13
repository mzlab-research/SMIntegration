########download

output$download_pseudocolor_data <- downloadHandler(
  filename = function() {

    paste0("pseudocolor_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    Co_visualisationlist<-cor_gene_metabolite_Co_visualisation_plot_save()
    datalist<-Co_visualisationlist[[2]]
    DemoData=datalist[[1]]
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')

  })
})
output$download_pseudocolor_plot <- downloadHandler(
  filename = function() {
    "pseudocolor_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- get("pseudocolor_plot", envir = .GlobalEnv)

    ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
    
    })
  })

output$download_Co_visualisation_data <- downloadHandler(
  filename = function() {

    paste0("co-expression_pattern_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      Co_visualisationlist<-cor_gene_metabolite_Co_visualisation_plot_save()
    datalist<-Co_visualisationlist[[2]]
    DemoData=datalist[[2]]
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    
  })
})
output$download_Co_visualisation_plot <- downloadHandler(
  filename = function() {
    "co-expression_pattern_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- get("Co_visualisation_plot", envir = .GlobalEnv)

    ggsave(file, plot = p, device = "png", width = 10, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })



##############test
# data_rds<- reactive({
#   data_rds<-readRDS("data_rds.rds")
#   return(data_rds)
# })
#########

spatial_coordlist<-eventReactive(c(input$start_Co_visualisation_analysis), {
  req(data_rds())

  withProgress(message = "Processing data...",value=0.8,{
    print("spatial_coordlist start")
    data_rds<-data_rds()
    data_mrds<-data_rds[[1]]
    data_trds<-data_rds[[2]]
    decon_mtrx = data_mrds@assays$Spatial$counts
    decon_ttrx = data_trds@assays$Spatial$counts
    stopifnot(identical(colnames(decon_mtrx), colnames(decon_ttrx)))
    combined_matrix <- rbind2(decon_mtrx, decon_ttrx)
    meta.data<-data_mrds@meta.data

  source("./source/CorrelationAnalysisFunction/Co_visualisation/Co_visualisation.R")

  combine_spatial_coord<-spatial_coord_processing(combined_matrix,meta.data)
  spatial_coordlist<-combine_spatial_coord
  print("spatial_coordlist done")
return(spatial_coordlist)
 })
})
spatial_coord_features<-reactive({
  req(data_rds())

  withProgress(message = "Processing data...",value=0.8,{
    print("spatial_coordlist start")
    data_m<-data_rds()[[1]]
    data_t<-data_rds()[[2]]
    gene<-unique(as.character(rownames(data_t@assays$Spatial$counts)))
    gene <- sort(gene)
    req(data_m)
    mz<-unique(as.character(rownames(data_m@assays$Spatial$counts)))
    spatial_coord_features<-c(mz,gene)
    return(spatial_coord_features)
  })
})

cor_gene_metabolite_Co_visualisation_plot_save<-eventReactive(c(input$start_Co_visualisation_analysis,input$cor_F1_select,input$cor_F2_select,input$cor_F3_select), {

  req(spatial_coordlist())
  withProgress(message = "Processing data...",value=0.8,{
  print("cor_gene_metabolite_Co_visualisation_plot_save start")
   spatial_coord<-spatial_coordlist()

   source("./source/CorrelationAnalysisFunction/Co_visualisation/Co_visualisation.R")
   
   if (!is.null(input$cor_F1_select) && input$cor_F1_select!=""){
     F1=input$cor_F1_select
   }else{
     F1=NA
   }
   if (!is.null(input$cor_F2_select) && input$cor_F2_select!=""){
     F2=input$cor_F2_select
   }else{
     F2=NA
   }
   if (!is.null(input$cor_F3_select) && input$cor_F3_select!=""){
     F3=input$cor_F3_select
   }else{
     F3=NA
   }
     pair<-c(F1,F2,F3)
     print(pair)
     if(sum(is.na(pair))==3){
       showNotification("Please select at least one feature!")
      return(NULL)
     }

     p<-Co_visualisation_plot(spatial_coord,pair)
     p0<-multiple_ions_plot(spatial_coord,pair)
     plot<-list(p0[[1]],p[[1]])
     data<-list(p0[[2]],p[[2]])
     all<-list(plot,data)

     return(all)
     print("cor_gene_metabolite_Co_visualisation_plot_save done")

  })
})

observeEvent(c(input$start_Co_visualisation_analysis,input$cor_F1_select,input$cor_F2_select,input$cor_F3_select),{
  req(cor_gene_metabolite_Co_visualisation_plot_save())
  plotlist<-cor_gene_metabolite_Co_visualisation_plot_save()[[1]]
  mplot<-plotlist[[1]]
  tplot<-plotlist[[2]]

  output$pseudocolor_plot <- renderPlot({

    withProgress(message = "Plotting...",value=0.8,{
      
      p<-gridExtra::grid.arrange(mplot,ncol=1)
      p
      assign("pseudocolor_plot", p, envir = .GlobalEnv)
    })
  })
  output$cor_gene_metabolite_Co_visualisation_plot <- renderPlot({

    withProgress(message = "Plotting...",value=0.8,{
      p<-gridExtra::grid.arrange(tplot,ncol=1)
      p
      assign("Co_visualisation_plot", p, envir = .GlobalEnv)
    })
  })

})

observe({
  req(spatial_coord_features())
  f=spatial_coord_features()

  updateSelectizeInput(
    session = getDefaultReactiveDomain(),
    inputId = "cor_F1_select",
    choices = f,
    server = TRUE,
    selected = if (length(f) > 0) f[1] else NULL,
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
  updateSelectizeInput(
    session = getDefaultReactiveDomain(),
    inputId = "cor_F2_select",
    choices = f,
    server = TRUE,
    selected = if (length(f) > 1) f[2] else NULL,
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
  updateSelectizeInput(
    session = getDefaultReactiveDomain(),
    inputId = "cor_F3_select",
    choices = f,
    server = TRUE,
    selected = if (length(f) > 2) f[3] else NULL,
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
