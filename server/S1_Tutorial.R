
demometab <- reactive({

  data <- fread("./example_data/metab_bin_filter.txt")
  
  return(data)
})
demotrans <- reactive({

  data <- fread("./example_data/trans_bin_filter.txt")
  return(data)
})
demometab_rds <- reactive({

  data <- readRDS("./example_data/pre_metab1.rds")

  return(data)
})
demotrans_rds <- reactive({

  data <- readRDS("./example_data/pre_trans1.rds")

  return(data)
})

#downloadData--------------------------------
output$download_demometab <- downloadHandler(
  filename = function() {
    paste0("spatial_metabolomics_demo.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    DemoData<- demometab()
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_demotrans <- downloadHandler(
  filename = function() {
    paste0("spatial_transcriptomics_demo.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    DemoData<-demotrans()
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
  })
})

output$download_demometab_rds <- downloadHandler(
  filename = function() {
    paste0("spatial_metabolomics_demo.rds")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    DemoData<- demometab_rds()
    saveRDS(DemoData, file)
    
  })
})
output$download_demotrans_rds <- downloadHandler(
  filename = function() {
    paste0("spatial_transcriptomics_demo.rds")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    DemoData<-demotrans_rds()
    saveRDS(DemoData, file)
    
  })
})


