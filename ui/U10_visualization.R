# ==============================================================================
# U10_visualization.R
# UI definition for the "Single Molecule Spatial Imaging" tab (Step 6 Part 1).
#
# Purpose:
#   Provides the interface for exploring the spatial distribution of individual molecules (Genes or Metabolites).
#   Automatically identifies and visualizes top correlated features (co-localization).
#
# Key Features:
#   - Feature Selection: Choose modality (Gene/Metabolite) and specific feature name.
#   - Primary Visualization: Spatial heatmap of the selected feature.
#   - Correlation Analysis:
#     - Displays top 6 positively/negatively correlated Metabolites.
#     - Displays top 6 positively/negatively correlated Genes.
#   - Export: Download buttons for all plots and underlying data.
#
# Structure:
#   - Control box for feature selection.
#   - Output box for the main spatial plot.
#   - Grid of 4 boxes showing correlated features (Pos/Neg for Metabolites/Genes).
# ==============================================================================

#' @title Single Feature Visualization UI
#' @description Defines the layout for single-molecule spatial analysis.
#' Includes a selector for the target feature and panels displaying its distribution 
#' along with the most correlated neighbors (spatial co-expression/co-abundance).
tabItem(tabName = "visual_sin",
        fluidRow(
          h2("Single Molecule Spatial Imaging")
        ),
        fluidRow(
          box(
            title = "Molecular Feature Selection",
            p("Step 1: Select molecular feature type and specific gene or metabolite for spatial visualization."),
            p("Step 2: Visualize spatial distribution of selected feature."),
            p("Step 3: Identify the top 6 positively and negatively correlated metabolites and genes, ranked by the strength of their correlation with the selected feature."),
            width = 6,
            style = "height: 400px; overflow-y: auto;",
            selectInput("single_ions_type", "Feature type:", choices = c("Gene","Metabolite"),selected = "Gene"),
            selectizeInput("select_single_ions", "Select feature:", choices = NULL, selected = NULL,options = list(server = TRUE)),
            actionButton("start_visualisation_analysis", "Generate visualization")
            
          ),
          box(
            title = "Spatial Distribution Visualization",
            tags$ul(
              tags$li("Color gradient: Molecular abundance (red: high, blue: low)")
            ),
            width = 6,
            style = "height: 400px; overflow-y: auto;",
            plotOutput("single_ion_visualizaiton", width = "250px", height = "250px"),
            column(3,downloadButton("download_single_feature_plot", "Download image")),
            column(3,downloadButton("download_single_feature_data", "Export data"))
          )
        ),
        
        fluidRow(
          box(
            title = "Top Spatially Co-localized Metabolites",
            width = 6,
            p("Identify spatially co-varying metabolites:"),
            tags$ul(
              tags$li("Positive correlation may suggests functional association"),
              tags$li("Top 6 metabolites by correlation strength")
            ),
            style = "height: 650px; overflow-y: auto;",
            plotOutput("positively_correlated_ions_m_plot", width = "700px", height = "500px"),
            column(3,downloadButton("download_positively_correlated_ions_m_plot", "Download image")),
            column(3,downloadButton("download_positively_correlated_ions_m_data", "Export data"))
          ),
          box(
            title = "Top Spatially Anti-correlated Metabolites",
            width = 6,
            p("Identify spatially anti-correlated metabolites:"),
            tags$ul(
              tags$li("Negative correlation may indicate inhibitory relationships"),
              tags$li("Top 6 metabolites by correlation strength")
            ),
            style = "height: 650px; overflow-y: auto;",
            plotOutput("negatively_correlated_ions_m_plot", width = "700px", height = "500px"),
            column(3,downloadButton("download_negatively_correlated_ions_m_plot", "Download image")),
            column(3,downloadButton("download_negatively_correlated_ions_m_data", "Export data"))
          ),
          box(
            title = "Top Spatially Co-expressed Genes",
            width = 6,
            style = "height: 650px; overflow-y: auto;",
            p("Identify spatially co-expressed genes:"),
            tags$ul(
              tags$li("Positive correlation may suggests functional association"),
              tags$li("Top 6 genes by correlation strength")
            ),
            plotOutput("positively_correlated_ions_t_plot", width = "700px", height = "500px"),
            column(3,downloadButton("download_positively_correlated_ions_t_plot", "Download image")),
            column(3,downloadButton("download_positively_correlated_ions_t_data", "Export data"))
          ),
          box(
            title = "Top Spatially Anti-expressed Genes",
            width = 6,
            style = "height: 650px; overflow-y: auto;",
            p("Identify spatially anti-expressed genes:"),
            tags$ul(
              tags$li("Negative correlation may indicate inhibitory relationships"),
              tags$li("Top 6 genes by correlation strength")
            ),
            plotOutput("negatively_correlated_ions_t_plot", width = "700px", height = "500px"),
            column(3,downloadButton("download_negatively_correlated_ions_t_plot", "Download image")),
            column(3,downloadButton("download_negatively_correlated_ions_t_data", "Export data"))
          )
        )
)#tabitem