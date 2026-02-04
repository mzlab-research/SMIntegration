# ==============================================================================
# U4_clustering.R
# UI definition for the "Clustering Analysis and Cell Annotation" tab (Step 3).
#
# Purpose:
#   Provides the interface for identifying spatial domains (clusters) using multi-modal data.
#
# Key Features:
#   - Preprocessing Configuration: Options for Normalization (TIC/RMS) and Transformation (LogNormalize).
#   - Data Assessment: Tools (Elbow Plot, Spatial PC Plots) to guide parameter selection.
#   - Algorithm Selection: Choice of 5 clustering algorithms (LV, LM, SLM, Kmeans variants).
#   - Visualization:
#     - Spatial Cluster Mapping: Side-by-side view of clusters from Metabolomics, Transcriptomics, and Integrated data.
#     - Cross-Modal Correspondence: Sankey diagram showing relationships between the different clustering results.
#
# Structure:
#   - Detailed description of the pipeline.
#   - Control boxes for Preprocessing (Step 3.1) and Clustering (Step 3.2).
#   - Collapsible box for Data Assessment results.
#   - Output boxes for Spatial Maps and Sankey Diagrams.
# ==============================================================================

#' @title Clustering Analysis UI
#' @description Defines the layout for spatial domain identification.
#' Includes preprocessing controls, PCA assessment, algorithm selection, and cluster visualization.
tabItem(tabName = "overall_cluster",
        fluidRow(
          h2("Step3: Clustering Analysis and Cell Annotation"),
          h3("Spatial Domain Identification"),
          p("This module identifies biologically distinct tissue regions through multi-modal clustering. The pipeline includes:"),
          tags$ol(
            tags$li(strong("Preprocessing:"), 
                    tags$ul(
                      tags$li("Metabolomics: (1) Optional TIC or RMS normalization to correct pixel-wise intensity; (2) Log-transformation to stabilize variance; and (3) Scaling with variable feature selection."),
                      tags$li("Transcriptomics: Standard LogNormalize, top 2000 HVG selection, and Scaling.")
                    )
            ),
            tags$li(strong("Data Integration:"), "Integration of transcriptomic/metabolomic matrices"),
            tags$li(strong("Dimensionality Reduction:"), "PCA/UMAP on top 30 principal components. PCA is utilized to reduce data complexity. Visualizations include Elbow Plots to aid in optimal dimension selection by observing the 'knee' point, and spatial plots of top PCs (PC1-3) to identify global spatial patterns."),
            tags$li(strong("Clustering:"), "Five algorithms: Seurat-LV (original Louvain algorithm), Seurat-LM (Louvain algorithm with multilevel refinement), 
          Seurat-SLM (Smart Local Moving algorithm), UMAP-Kmeans, PCA-Kmeans"),
            tags$li(strong("Visualizing:"), "Visualize the relationships among clustering results of three datasets (spatial transcriptomics, spatial metabolomics, and multi-omics integration) through Sankey diagrams")
            
          ),
          p("Select a method and click 'Start clustering computation' to identify spatial domains across three data views. Adjust resolution to control cluster granularity (higher values = more clusters).")
        ),
        fluidRow(
          box(
            title = "Step 3.1: Preprocessing & Data Assessment",
            width = 6,
            status = "primary",
            solidHeader = TRUE,
            p("Configure the preprocessing pipeline for spatial metabolomics data."),
            selectInput("metab_norm_method", "Normalization Method:",
                        choices = c("None" = "None",
                                    "Total Ion Current (TIC)" = "TIC",
                                    "Root Mean Squared (RMS)" = "RMS"),
                        selected = "TIC"),
            selectInput("metab_transform_method", "Transformation:",
                        choices = c("LogNormalize" = "LogNormalize",
                                    "None" = "None"),
                        selected = "LogNormalize"),
            helpText("These settings will be applied before clustering."),
            hr(),
            p("Run this first to determine optimal parameters:"),
            tags$ul(
              tags$li("Use Elbow Plot to visualize data variance captured by top 30 PCs and estimate optimal cluster number.")
            ),
            actionButton("run_assessment", "Run Data Assessment (PCA)", class = "btn-warning", width = "100%")
          ),
          box(
            title = "Step 3.2: Clustering Algorithm Selection",
            width = 6,
            status = "primary",
            solidHeader = TRUE,
            selectizeInput("cluster_select", "Select clustering algorithm", choices = c("LV","LM","SLM","UMAP-kmeans","PCA-kmeans"), selected = "LV"),
            uiOutput("culster_resolution_button_container"),
            p("Note: Computational time varies with data size. Please avoid duplicate submissions."),
            actionButton("start_cluster", "Start clustering computation"),
          ),
          column(12, br()),
          box(
            title = "Data Assessment Results (PCA)",
            width = 12,
            collapsible = TRUE,
            collapsed = FALSE,
            p("Visualizations to guide parameter selection and assess data quality."),
            tabsetPanel(
              tabPanel("Variance (Elbow Plot)",
                br(),
                p("The elbow plot shows the standard deviation of each principal component. Observe the “knee” in the plot, where the standard deviation to level off, to help estimate the optimal cluster number."),
                plotOutput("pca_elbow_plot", width = "1200px", height = "400px"),
                downloadButton("download_pca_elbow_plot", "Download Elbow Plot")
              ),
              tabPanel("Spatial Patterns (Top 3 PCs)",
                br(),
                plotOutput("pca_spatial_plot", width = "1200px", height = "1000px"),
                downloadButton("download_pca_spatial_plot", "Download Spatial Plot")
              )
            )
          ),
          column(12, br()) ,
          box(
            title = "Spatial Cluster Mapping",
            p("Visualization of clustering results mapped to original tissue architecture:"),
            tags$ul(
              tags$li(strong("Panels:"), "Metabolomics-only (left), Transcriptomics-only (center), Integrated multi-omics (right)")
            ),
            style = "height: 580px; overflow-y: auto;",
            width = 12,
            plotOutput("cluster_plot", width = "1200px", height = "400px"),
            column(3,downloadButton("download_cluster_plot", "Download image")),
            column(3,downloadButton("download_cluster_data", "Export cluster assignments"))
          ),
          column(12, br()) ,
          box(
            title = "Cross-Modal Cluster Correspondence",
            p("Sankey diagram revealing relationships between clustering solutions:"),
            tags$ul(
              tags$li(strong("Left:"), "Metabolomics-derived clusters"),
              tags$li(strong("Center:"), "Integrated multi-omics clusters"),
              tags$li(strong("Right:"), "Transcriptomics-derived clusters"),
              tags$li(strong("Connections:"), "Proportional flow between cluster assignments")
            ),
            width = 12,
            style = "height: 700px; overflow-y: auto;",
            htmlOutput("cluster_sanky_plot", width = "800px", height = "400px"),
            column(3,downloadButton("download_cluster_sanky_plot", "Download image")),
            column(3,downloadButton("download_cluster_sanky_plot_data", "Export data"))
          )
        )
        
)  # tabitem