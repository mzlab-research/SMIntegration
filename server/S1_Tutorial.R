# ==============================================================================
# S1_Tutorial.R
# Server logic for the Tutorial/Introduction tab.
# Handles loading of built-in demo datasets (TXT and RDS formats) and download handlers for tutorial materials.
# ==============================================================================

# ------------------------------------------------------------------------------
# Demo Data Loaders
# ------------------------------------------------------------------------------

#' @title Load Demo Metabolomics Data (TXT)
#' @description Reads the pre-processed spatial metabolomics demo data in text format.
demometab <- reactive({
  data <- fread("./example_data/metab_bin_filter.txt")
  return(data)
})

#' @title Load Demo Transcriptomics Data (TXT)
demotrans <- reactive({
  data <- fread("./example_data/trans_bin_filter.txt")
  return(data)
})


#' @title Load Demo Metabolomics Data (RDS/Seurat)
#' @description Reads the pre-processed spatial metabolomics demo data as a Seurat object.
demometab_rds <- reactive({
  data <- readRDS("./example_data/pre_metab1.rds")
  return(data)
})

#' @title Load Demo Transcriptomics Data (RDS/Seurat)
demotrans_rds <- reactive({
  data <- readRDS("./example_data/pre_trans1.rds")
  return(data)
})

#' @title Load Unregistered Metabolomics Data (RDS/Seurat)
noreg_demometab_rds <- reactive({
  data <- readRDS("./example_data/mbefore_metabolite.rds")
  return(data)
})

# ------------------------------------------------------------------------------
# Download Handlers
# ------------------------------------------------------------------------------

#' @title Download Demo Metabolomics (TXT)
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

#' @title Download Demo Transcriptomics (TXT)
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

#' @title Download Demo Metabolomics (RDS)
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

#' @title Download Demo Transcriptomics (RDS)
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

#' @title Download Registration Tutorial PDF
output$download_registration_tutorial <- downloadHandler(
  filename = "registration_tutorial.pdf",
  content = function(file) {
    file.copy("www/registration_tutorial.pdf", file)
  }
)