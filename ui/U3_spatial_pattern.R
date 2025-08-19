tabItem(tabName = "Spatial_Pattern",
        fluidRow(
          h2("Step2: Spatial Pattern Analysis"),
          p("This module employs the SpaGene algorithm to identify spatially co-varying molecular modules and discover cross-modal pattern correlations."),
          tags$ul(
            tags$li("Step 1: Construct k-nearest neighbor spatial networks for each omics layer"),
            tags$li("Step 2: Identify spatially variable features by comparing actual Earth mover's distance (EMDg) against a null distribution from 500 random permutations"),
            tags$li("Step 3: Cluster significant features into spatial modules"),
            tags$li("Step 4: Quantify cross-omics pattern correlations using Moran's I (spdep)"),
            tags$li("Step 5: Screen pattern-specific features")
            
          ),
          p("Click 'Start spatial pattern analysis' to initiate the workflow.")
        ),
        fluidRow(
          box(
            width = 6,
            p("Note: Computational time varies with feature count. Please avoid duplicate submissions."),
            actionButton("start_Spatial_Pattern", "Start spatial pattern analysis")
          )
        ),
        column(12, br()),
        fluidRow(
          h3("Spatial Pattern Modules"),
          p("Identified spatial pattern modules for metabolites (left) and genes (right). "),
          box(
            title = "Metabolite Spatial Pattern Modules",
            width = 6,
            style = "height: 520px; overflow-y: auto;",
            plotOutput("spatial_pattern_plot_m", width = "600px", height = "450px"),
            column(3, downloadButton("download_spatial_spatial_pattern_polt_m", "Download image"))
          ),
          box(
            title = "Gene Spatial Pattern Modules",
            width = 6,
            style = "height: 520px; overflow-y: auto;",
            plotOutput("spatial_pattern_plot_t", width = "600px", height = "450px"),
            column(3, downloadButton("download_spatial_spatial_pattern_polt_t", "Download image"))
          )
        ),
        fluidRow(
          h3("Cross-Omics Pattern Correlation"),
          p("Moran's I correlation matrix quantifying spatial co-variation between metabolite and gene pattern modules. Significance levels: *p<0.05, **p<0.01, ***p<0.001."),
          box(
            width = 6,
            style = "height: 600px; overflow-y: auto;",
            div(
              style = "display: flex; justify-content: center; align-items: center; height: 500px;",
              plotOutput("spatial_heatmap_plot", width = "100%", height = "100%")
            ),
            column(3, downloadButton("download_spatial_heatmap_plot", "Download image")),
            column(3, downloadButton("download_spatial_heatmap_data", "Export correlation matrix"))
          )
        ),
        fluidRow(
          h3("Pattern-Specific Feature Classification"),
          p("Features are assigned to spatial modules based on spatial correlation:"),
          tags$ul(
            tags$li("Correlation ≥0.5: Assigned to best-matching module"),
            tags$li("Correlation <0.5: Categorized as 'unclassified'")
          )

        ),
        fluidRow(
          box(
            title = "Metabolite Distribution Across Spatial Modules",
            width = 6,
            style = "height: 700px; overflow-y: auto;",
            p("Barplot showing metabolite counts per spatial pattern module after correlation thresholding."),
            plotOutput("Pattern_specific_plot_m", width = "600px", height = "600px"),
            column(3, downloadButton("download_Pattern_specific_polt_m", "Download image")),
            column(3, downloadButton("download_Pattern_specific_data_m", "Export metabolite assignments"))
          ),
          box(
            title = "Gene Distribution Across Spatial Modules",
            width = 6,
            style = "height: 700px; overflow-y: auto;",
            p("Barplot showing gene counts per spatial pattern module after correlation thresholding."),
            plotOutput("Pattern_specific_plot_t", width = "600px", height = "600px"),
            column(3, downloadButton("download_Pattern_specific_polt_t", "Download image")),
            column(3, downloadButton("download_Pattern_specific_data_t", "Export gene assignments"))
          )
        ),
        column(12, br()),
        fluidRow(
          h3("Spatial Distribution of Pattern-Specific Features"),
          p("Visualize spatial expression patterns of top module-associated features. Select modules and features below:"),
          p("The top 300 features per module (ranked by correlation) are designated as 'pattern-specific features' for functional analysis."),
          box(
            title = "Metabolite Spatial Distribution",
            width = 6,
            style = "height: 700px; overflow-y: auto;",
            selectInput("pattern_select_m", "Select metabolite spatial module:", choices = NULL, selected = NULL),
            selectizeInput("Pattern_specific_m_select", "Select metabolite:", choices = NULL, selected = NULL,options = list(server = TRUE)),
            plotOutput("patternspecific_distribution_plot", width = "400px", height = "400px"),
            column(3, downloadButton("download_Pattern_specific_m_plot", "Download image")),
            column(3, downloadButton("download_Pattern_specific_m_data", "Export intensity data"))
          ),
          box(
            title = "Gene Spatial Distribution",
            width = 6,
            style = "height: 700px; overflow-y: auto;",
            selectInput("pattern_select_t", "Select gene spatial module:", choices = NULL, selected = NULL),
            selectizeInput("Pattern_specific_t_select", "Select gene:", choices = NULL, selected = NULL,options = list(server = TRUE)),
            plotOutput("patternspecific_distribution_plot_t", width = "400px", height = "400px"),
            column(3, downloadButton("download_Pattern_specific_t_plot", "Download image")),
            column(3, downloadButton("download_Pattern_specific_t_data", "Export expression data"))
          )
        )
)#tabitem
