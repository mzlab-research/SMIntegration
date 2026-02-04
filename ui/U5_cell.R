

# ==============================================================================
# U5_cell.R
# UI definition for the "Cell Annotation" sub-tab (Step 3.3).
#
# Purpose:
#   Provides the interface for assigning cell types to spatial spots/cells.
#
# Key Features:
#   - Annotation Methods:
#     - SingleR: Automated annotation using reference databases (Mouse/Human).
#     - Manual Upload: Users can upload a custom annotation table.
#   - Visualization: Spatial plot showing the distribution of annotated cell types.
#   - Data Propagation: Notes that annotations are transferred to metabolomics data.
#
# Structure:
#   - Parameter selection box (Method, File Upload).
#   - Visualization box for the spatial cell type map.
# ==============================================================================

#' @title Cell Annotation UI
#' @description Defines the layout for cell type annotation.
#' Includes controls for selecting the annotation method (SingleR/Upload) and visualizing the results.
tabItem(tabName = "cell_anno",
        fluidRow(
          h2("Step3: Clustering Analysis and Cell Annotation")
        ),
        fluidRow(
          h3("Cell Type Deconvolution"),
          p("Spatial transcriptomics data is annotated using:"),
          tags$ol(
            tags$li(strong("SingleR automated annotation:"), "Leverages reference atlases (MouseRNAseqData/HumanPrimaryCellAtlasData) for cell type prediction"),
            tags$li(strong("Custom annotation upload:"), "User-provided annotations mapping spatial coordinates to cell types (download template below)")
          ),
          p("Annotation results are propagated to metabolomics pixels for downstream cross-modal analysis.")
        ),	
        fluidRow(
          box(
            width = 6,
            style = "height: 520px; overflow-y: auto;",
            title = "Step 3.3: Annotation Parameters",
            status = "primary",
            solidHeader = TRUE,
            selectizeInput("cell_annotation_select", "Select annotation method", choices = c("SingleR","Upload_annotation_table"), selected = "Upload_annotation_table"),
            uiOutput("cell_annotation_button_container"),
            p("Note: If you are using the demo data, you can proceed directly by clicking 'Annotate Cell Types' without uploading a file."),
            p("Note: Computational time varies with data size. Please avoid duplicate submissions."),
            actionButton("start_cell_annotation", "Annotate Cell Types"),
          ),
          box(
            width = 6,
            style = "height: 520px; overflow-y: auto;",
            title = "Spatial Cell Type Mapping",
            tags$ul(
              tags$li("Color-coded cell type assignments")
            ),
            plotOutput("cell_annotation_plot", width = "100%", height = "400px"),
            column(3,downloadButton("download_cell_annotation_plot", "Download image")),
            column(3,downloadButton("download_cell_annotation_data", "Export annotations"))
          )
        )
        
)  # tabitem