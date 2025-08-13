tabItem(tabName = "diff_analysis",
        fluidRow(
          h2("Step4: Differential Analysis"),
          h3("Differential Feature Identification"),
          p("Identify differentially expressed genes (DEGs) and differentially abundant metabolites (DAMs) using Seurat's FindMarkers function:"),
          tags$ul(
            tags$li("Wilcoxon rank-sum test with Bonferroni correction"),
            tags$li("Thresholds: |log2FC| > 0.26 & adjusted p-value < 0.05"),
            tags$li("Pixel-level analysis treating each spatial spot as independent sample")
          ),
          p("Set thresholds below and click 'Start differential screening'.")
        ),
        fluidRow(
          box(
            width = 6,
            style = "height: 590px; overflow-y: auto;",
            numericInput("Fold_change","Metabolomics |log2FC| threshold:",min=0.1,max=2,value=0.26,step=0.01),
            numericInput("Fold_change_t","Transcriptomics |log2FC| threshold:",min=0.1,max=2,value=0.26,step=0.01),
            numericInput("p_value_adjust","Metabolomics FDR threshold:",min=0.001,max=1,value=0.05,step=0.01),
            numericInput("p_value_adjust_t","Transcriptomics FDR threshold:",min=0.001,max=1,value=0.05,step=0.01),
            p("Note: Computation scales with feature count. Please avoid duplicate submissions."),
            actionButton("start_diff_analysis", "Start differential screening"),
          ),
          box(
            title = "Differential Feature Summary",
            p("Barplot showing significantly upregulated (red) and downregulated (blue) features:"),
            tags$ul(
              tags$li("Height indicates feature counts"),
              tags$li("Annotations show total differential features per modality")
            ),
            width = 6,
            style = "height: 550px; overflow-y: auto;",
            plotOutput("diffbar_plot", width = "400px", height = "400px"),
            column(3,downloadButton("download_diffbar_plot", "Download image")),
            column(3,downloadButton("download_diff_omics", "Export differential features"))
          )
          ),
          fluidRow(
            h3("Differential Feature Characterization"),
            p("Multi-dimensional visualization of differential features:"),
            column(12, br()),
            box(
              title = "UMAP: Metabolite",
              tags$ul(
                tags$li("Points represent spatial spots"),
                tags$li("Color indicates treatment (red) vs control (blue) groups"),
                tags$li("Separation distance reflects group dissimilarity")
              ),
              width = 6,
              style = "height: 600px; overflow-y: auto;",
              plotOutput("umap_plot", width = "400px", height = "400px"),
              column(3,downloadButton("download_umap_plot", "Download image")),
              column(3,downloadButton("download_umap_data", "Export coordinates"))
            ),

            box(
              title = "UMAP: Gene",
              tags$ul(
                tags$li("Points represent spatial spots"),
                tags$li("Color indicates treatment (red) vs control (blue) groups"),
                tags$li("Separation distance reflects group dissimilarity")
              ),
              width = 6,
              style = "height: 600px; overflow-y: auto;",
              plotOutput("umap_plott", width = "400px", height = "400px"),
              column(3,downloadButton("download_umap_plott", "Download image")),
              column(3,downloadButton("download_umap_datat", "Export coordinates"))
            ),
            column(12, br()),
            box(
              title = "Volcano Plot: Metabolite",
              tags$ul(
                tags$li("X-axis: log2 fold change (effect size)"),
                tags$li("Y-axis: -log10(adjusted p-value) (significance)"),
                tags$li("Red: significantly upregulated features"),
                tags$li("Blue: significantly downregulated features")
              ),
              width = 6,
              style = "height: 600px; overflow-y: auto;",
              plotOutput("volcano_plot", width = "400px", height = "400px"),
              column(3,downloadButton("download_volcano_plot", "Download image")),
              column(3,downloadButton("download_volcano_data", "Export data"))
            ),
            box(
              title = "Volcano Plot: Gene",
              tags$ul(
                tags$li("X-axis: log2 fold change (effect size)"),
                tags$li("Y-axis: -log10(adjusted p-value) (significance)"),
                tags$li("Red: significantly upregulated features"),
                tags$li("Blue: significantly downregulated features")
              ),
              width = 6,
              style = "height: 600px; overflow-y: auto;",
              plotOutput("volcano_plott", width = "400px", height = "400px"),
              column(3,downloadButton("download_volcano_plott", "Download image")),
              column(3,downloadButton("download_volcano_datat", "Export data"))
            )
            ),
            fluidRow(
              h3("Spatial Distribution of Differential Features"),
              p("Examine spatial expression patterns of significant features:"),
              box(
                title = "Differential Metabolite Distribution",
                width = 6,
                style = "height: 550px; overflow-y: auto;",
                selectizeInput("diff_ion_select", "Select DAM to visualize:", choices = NULL,selected = NULL,options = list(server = TRUE)),
                plotOutput("diffion_distribution_plot", width = "400px", height = "400px"),
                column(3,downloadButton("download_diffion_distribution_plot", "Download image")),
                column(3,downloadButton("download_diffion_distribution_data", "Export intensity data"))
              ),
              box(
                title = "Differential Gene Distribution",
                width = 6,
                style = "height: 550px; overflow-y: auto;",
                selectizeInput("diff_ion_select_t", "Select DEG to visualize:", choices = NULL,selected = NULL,options = list(server = TRUE)),
                plotOutput("diffion_distribution_plot_t", width = "400px", height = "400px"),
                column(3,downloadButton("download_diffion_distribution_plot_t", "Download image")),
                column(3,downloadButton("download_diffion_distribution_data_t", "Export expression data"))
              )
            )
)#tabitem