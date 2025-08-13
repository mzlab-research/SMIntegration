## app.R ##
suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinydashboard))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(plotly))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(rjson))
suppressMessages(library(RColorBrewer))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(readr))
suppressMessages(library(data.table))
suppressMessages(library(gridExtra))
suppressMessages(library(zip))
suppressMessages(library(DT))
suppressMessages(library(viridis))
suppressMessages(library(networkD3))
suppressMessages(library(webshot))
suppressMessages(library(svglite))
suppressMessages(library(ggrepel))
suppressMessages(library(ade4))
suppressMessages(library(pheatmap))
suppressMessages(library(scales))
suppressMessages(library(tibble))
suppressMessages(library(mixOmics))
suppressMessages(library(caret))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(igraph))
suppressMessages(library(ggnewscale))
suppressMessages(library(grid))
suppressMessages(library(KEGGREST))
suppressMessages(library(purrr))
suppressMessages(library(pbapply))
suppressMessages(library(crayon))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(SpaGene))
suppressMessages(library(reshape2))
suppressMessages(library(parallel))
suppressMessages(library(Hmisc))
suppressMessages(library(VennDiagram))
suppressMessages(library(magick))
suppressMessages(library(SingleR))
suppressMessages(library(future))
suppressMessages(library(MetaboCoreUtils)) 
suppressMessages(library(enviPat)) 
suppressMessages(library(msentropy)) 
suppressMessages(library(Spectra)) 
suppressMessages(library(gtools))
suppressMessages(library(ggraph))
suppressMessages(library(igraph))
suppressMessages(library(spdep))
suppressMessages(library(sf))


options(ggrepel.max.overlaps = Inf)
options(encoding='UTF-8')
ui <- dashboardPage(
  dashboardHeader(
    title = span(icon("chart-area"), "Spatial Multi-omics Integration Platform (SMIntegration)"),
    titleWidth = 550,
    tags$li(class = "dropdown",
            tags$style(".main-header .logo {height: 50px; line-height: 50px;}")
    )
  ),
  dashboardSidebar(width = 350,
                   sidebarMenu(
                     id = "tabs",
                     menuItem("Tutorial", tabName = "Tutorial", icon = icon("dashboard")),
                     menuItem("Overall Distribution Analysis", tabName = "File_upload", icon = icon("chart-bar")),
                     menuItem("Spatial Pattern Analysis", tabName = "Spatial_Pattern", icon = icon("map-marked-alt")),
                     menuItem("Clustering Analysis and Cell Annotation", tabName = "Clustering", icon = icon("project-diagram"), startExpanded = TRUE,
                              menuSubItem("Clustering Analysis", tabName = "overall_cluster"),
                              menuSubItem("Cell Annotation", tabName = "cell_anno")
                     ),
                     menuItem("Differential Analysis", tabName = "Differential", icon = icon("balance-scale"), startExpanded = TRUE,
                              menuSubItem("Comparison Group Selection", tabName = "diff_select"),
                              menuSubItem("Differential Screening and Visualization", tabName = "diff_analysis"),
                              menuSubItem("Group-Specific Network", tabName = "con_net")
                     ),
                     # menuItem("Correlation Analysis", tabName = "diff_cor", icon = icon("th")),
                     
                     menuItem("Functional Association Analysis", tabName = "Pathway_analysis", icon = icon("sitemap")),
                     menuItem("Data Visualization", tabName = "datavi", icon = icon("image"), startExpanded = TRUE,
                              menuItem("Single feature visualization", tabName = "visual_sin"),
                              menuItem("Co-visualisation Analysis", tabName = "co_visual")
                     )
                   ),
                   tags$style(HTML("
      .sidebar-menu li a { 
        font-size: 15px;
        padding: 12px 15px;
      }
      .sidebar-menu .treeview-menu li a {
        padding: 10px 15px 10px 35px;
      }
    "))
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$title("SMIntegration"),
      tags$link(rel = "stylesheet", type = "text/css", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.3/css/all.min.css"),
      tags$style(HTML("
        .content-wrapper {
          background-color: #f9fbfd;
        }
        .box {
          border-radius: 8px;
          box-shadow: 0 2px 6px rgba(0,0,0,0.08);
          border-top: 3px solid #3c8dbc;
        }
        .box.box-primary {
          border-top-color: #3c8dbc;
        }
        .box.box-info {
          border-top-color: #00c0ef;
        }
        .box.box-success {
          border-top-color: #00a65a;
        }
        .box.box-warning {
          border-top-color: #f39c12;
        }
        .box.box-danger {
          border-top-color: #dd4b39;
        }
        h2, h3, h4 {
          color: #2c3e50;
          font-weight: 600;
        }
        .btn {
          border-radius: 4px;
          font-weight: 500;
        }
      "))
    ),
    tabItems(
      source(file.path("ui","U1_Tutorial.R"),  local = TRUE)$value,
      
      source(file.path("ui","U2_upload.R"),  local = TRUE)$value,
      source(file.path("ui","U3_spatial_pattern.R"),  local = TRUE)$value,
      source(file.path("ui","U4_clustering.R"),  local = TRUE)$value,
      source(file.path("ui","U5_cell.R"),  local = TRUE)$value,
      source(file.path("ui","U6_ROI_select.R"),  local = TRUE)$value,
      source(file.path("ui","U7_differential_analysis.R"),  local = TRUE)$value,
      source(file.path("ui","U8_network.R"),  local = TRUE)$value,
      # source(file.path("ui","U5.R"),  local = TRUE)$value,
      source(file.path("ui","U9_pathway.R"),  local = TRUE)$value,
      source(file.path("ui","U10_visualization.R"),  local = TRUE)$value,
      source(file.path("ui","U11_co_visualization.R"),  local = TRUE)$value
      
      
      
    )#tabItems
  )#dashboardBody
)#ui

server <- function(input, output) {
  
  options(shiny.maxRequestSize=10 * 1024^3)
  options(future.globals.maxSize = 30 * 1024^3)#30GB
  source(file.path("server","S1_Tutorial.R"),  local = TRUE)$value
  
  source(file.path("server","S2_upload.R"),  local = TRUE)$value
  source(file.path("server","S3_spatial_pattern.R"),  local = TRUE)$value
  source(file.path("server","S4_clustering.R"),  local = TRUE)$value
  source(file.path("server","S5_cell.R"),  local = TRUE)$value
  source(file.path("server","S6_ROI_select.R"),  local = TRUE)$value
  source(file.path("server","S7_differential_analysis.R"),  local = TRUE)$value
  source(file.path("server","S8_network.R"),  local = TRUE)$value
  source(file.path("server","S9_pathway.R"),  local = TRUE)$value
  source(file.path("server","S10_visualization.R"),  local = TRUE)$value
  source(file.path("server","S11_co_visualization.R"),  local = TRUE)$value
  
}
shinyApp(ui, server)