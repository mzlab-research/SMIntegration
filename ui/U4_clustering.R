tabItem(tabName = "overall_cluster",
        fluidRow(
          h2("Step3: Clustering Analysis and Cell Annotation"),
          h3("Spatial Domain Identification"),
          p("This module identifies biologically distinct tissue regions through multi-modal clustering. The pipeline includes:"),
          tags$ol(
            tags$li(strong("Preprocessing:"), "LogNormalize scaling, top 2000 HVG/HVM selection, and feature scaling via ScaleData"),
            tags$li(strong("Data Integration:"), "Integration of transcriptomic/metabolomic matrices"),
            tags$li(strong("Dimensionality Reduction:"), "PCA/UMAP on top 30 principal components"),
            tags$li(strong("Clustering:"), "Four algorithms: Seurat-LV (original Louvain algorithm), Seurat-LM (Louvain algorithm with multilevel refinement), 
          Seurat-SLM (Smart Local Moving algorithm), UMAP-Kmeans"),
            tags$li(strong("Visualizing:"), "Visualize the relationships among clustering results of three datasets (spatial transcriptomics, spatial metabolomics, and multi-omics integration) through Sankey diagrams")
            
          ),
          p("Select a method and click 'Start clustering computation' to identify spatial domains across three data views. Adjust resolution to control cluster granularity (higher values = more clusters).")
        ),
        fluidRow(
          box(
            width = 6,
            selectizeInput("cluster_select", "Select clustering algorithm", choices = c("LV","LM","SLM","UMAP-kmeans"), selected = "LV"),
            uiOutput("culster_resolution_button_container"),
            p("Note: Computational time varies with data size. Please avoid duplicate submissions."),
            actionButton("start_cluster", "Start clustering computation"),
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