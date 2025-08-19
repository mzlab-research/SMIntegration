tabItem(tabName = "Tutorial",
        fluidRow(

          h2("SMIntegration: Spatial Multi-omics Integration Platform"),
          p(strong("SMIntegration"), "is an innovative open-source platform for integrated analysis of spatial transcriptomics and metabolomics data. It addresses critical challenges in spatial multi-omics by providing:"),
          tags$ul(
            tags$li(strong("Comprehensive integrated spatial analysis:"), "Joint analysis pipeline from initial visualization to biological interpretation"),
            tags$li(strong("Flexible differential analysis:"), "Comparative analysis supporting multiple region definition strategies"),
            tags$li(strong("Multimodal pattern discovery:"), "Multiple algorithms for identifying spatially co-varying regions and gene-metabolite modules"),
            tags$li(strong("Flexible differential analysis:"), "Comparative analysis supporting multiple region definition strategies"),
            tags$li(strong("Inter-group network comparison:"), "Construction and comparison of condition-specific correlation networks for differential molecules (DEGs/DAMs)"),
            tags$li(strong("Interactive visualization:"), "Dynamic exploration of spatial distributions and co-localization patterns")
          ),
          

          h2("Core Analytical Modules"),
          tags$div(
            style = "display: grid; grid-template-columns: repeat(2, 1fr); gap: 15px;",
            tags$div(
              h4(icon("chart-bar"), " Overall Distribution Analysis"),
              p("Initial data visualization of spatial distributions")
            ),
            tags$div(
              h4(icon("map-marked-alt"), " Spatial Pattern Analysis"),
              p("Identification of spatially co-varying gene-metabolite modules")
            ),
            tags$div(
              h4(icon("project-diagram"), " Clustering Analysis and Cell Annotation"),
              p("Spatial domain detection with multiple algorithms and automated cell type mapping")
            ),
            tags$div(
              h4(icon("balance-scale"), " Differential Analysis"),
              p("DEGs/DAMs discovery and region-specific network construction")
            ),
            tags$div(
              h4(icon("sitemap"), " Functional Association Analysis"),
              p("Pathway enrichment and cross-omics functional integration")
            ),
            tags$div(
              h4(icon("image"), " Data Visualization"),
              p("Advanced spatial imaging including co-localization and multi-feature views")
            )
          )
          ),
          fluidRow(
            column(6, 
                   h3("Analysis Workflow")
            ),
            column(6,
                   h3("Typical Analysis Workflows")
            )
          ),

         fluidRow(
          box(width = 6,
              tags$img(src = "workflow.png", width = "100%", style = "border: 1px solid #eee;"),
              tags$div(
                class = "callout callout-info",
                p("Figure: End-to-end analysis pipeline reflecting the main menu structure")
              )
              ),


          box(
            width = 6,
            tags$div(
              class = "panel-group",
              tags$div(
                class = "panel panel-default",
                tags$div(
                  class = "panel-heading",
                  h4("Workflow 1: Comprehensive Spatial Profiling", class = "panel-title")
                ),
                tags$div(
                  class = "panel-body",
                  tags$ol(
                    tags$li("Upload data in ", strong("Overall Distribution Analysis")),
                    tags$li("Identify patterns in ", strong("Spatial Pattern Analysis")),
                    tags$li("Perform ", strong("Clustering and Cell Annotation")),
                    tags$li("Compare regions in ", strong("Differential Analysis")),
                    tags$li("Interpret results in ", strong("Functional Association Analysis")),
                    tags$li("Visualize findings in ", strong("Data Visualization"))
                  )
                )
              ),
              tags$div(
                class = "panel panel-default",
                tags$div(
                  class = "panel-heading",
                  h4("Workflow 2: Targeted Region Comparison", class = "panel-title")
                ),
                tags$div(
                  class = "panel-body",
                  tags$ol(
                    tags$li("Define ROIs in ", strong("Differential Analysis > Comparison Group Selection")),
                    tags$li("Identify DEGs/DAMs in ", strong("Differential Analysis > Differential Screening")),
                    tags$li("Build networks in ", strong("Differential Analysis > Group-Specific Network")),
                    tags$li("Interpret results in ", strong("Functional Association Analysis"))
                  )
                )
              )
            )
          ),

          column(6, h3("Demo Datasets")),
          box(
            width = 6,
            p("Preprocessed mouse brain data demonstrating key features:"),
            tags$ul(
              tags$li(strong("Tissue:"), "Coronal section"),
              tags$li(strong("Transcriptomics:"), "500 HVGs from Stereo-seq (0.5μm → 50μm binned)"),
              tags$li(strong("Metabolomics:"), "500 HVMs from AFAI-MSI (positive mode)"),
              tags$li(strong("Preprocessing:"), "Full spatial registration and identification")
            ) 
          ),
          ),
          fluidRow(

          h2("Data Preparation Requirements"),
          tags$div(
            class = "alert alert-warning",
            h4("Critical Preprocessing Steps:"),
            tags$ul( 
              tags$li(strong("Spatial Registration:"), "Please align coordinates using ", 
                      tags$a(href="https://github.com/mzlab-research/SMIntegration/SpatialData", "SpatialData"), 
                      " (scripts provided on GitHub)"),
              tags$li(strong("Resolution Harmonization:"), "Aggregate higher-resolution data (e.g., bin100 for 500nm→50μm conversion)"),
              tags$li(strong("Coordinate System:"), "Identical coordinate units and origin point for both modalities"),
              tags$li(strong("Metabolite Identification:"), "Requires annotated metabolite names (not m/z values)")
            )
          )),
        fluidRow(

          h2("Input File Specifications"),
          h3("Option 1: Seurat Object Format (Recommended)"),
          p("You need to upload two separate Seurat objects (.rds files):"),
          box(width = 8,
              style = "height: 320px; overflow-y: auto;",

                tags$div(
                  class = "col-md-6",
                  h4("Spatial Metabolomics Object"),
                  tags$table(
                    class = "table table-bordered",
                    tags$thead(tags$tr(
                      tags$th("Component"), tags$th("Content"), tags$th("Slot Name")
                    )),
                    tags$tbody(
                      tags$tr(tags$td("Assay"), tags$td("Metabolite intensity matrix"), tags$td("Spatial")),
                      tags$tr(tags$td("Meta Data"), tags$td("Spatial coordinates"), tags$td("x, y in meta.data")),
                      tags$tr(tags$td("Active Assay"), tags$td("Must be named 'Spatial'"), tags$td("active.assay"))
                    )
                  )
                ),
                tags$div(
                  class = "col-md-6",
                  h4("Spatial Transcriptomics Object"),
                  tags$table(
                    class = "table table-bordered",
                    tags$thead(tags$tr(
                      tags$th("Component"), tags$th("Content"), tags$th("Slot Name")
                    )),
                    tags$tbody(
                      tags$tr(tags$td("Assay"), tags$td("Gene expression matrix"), tags$td("Spatial")),
                      tags$tr(tags$td("Meta Data"), tags$td("Spatial coordinates"), tags$td("x, y in meta.data")),
                      tags$tr(tags$td("Active Assay"), tags$td("Must be named 'Spatial'"), tags$td("active.assay"))
                    )
                  )
                ),
              downloadButton("download_demometab_rds", "Download metab Seurat demo", class = "btn-sm"),
              downloadButton("download_demotrans_rds", "Download trans Seurat demo", class = "btn-sm")
          ),
            box(
              width = 4,
              style = "height: 320px; overflow-y: auto;",
              tags$img(src = "feature_matrix_rds.png", width = "100%", style = "max-width: 400px;"),

            )
          ),

        fluidRow(
          
          tags$div(
            class = "alert alert-info",
            h4("Important Note:"),
            p("For both formats, the assay name must be 'Spatial' in each Seurat object. 
    The platform expects this naming convention to correctly identify the data modality.")
          ),
          hr(),
          h3("Option 2: Text Matrix Format"),
          p("You need to upload two separate files:"),
          box(width = 8,
              style = "height: 320px; overflow-y: auto;",

            tags$div(
              class = "col-md-6",
              h4("Spatial Metabolomics File"),
              tags$table(
                class = "table table-bordered",
                tags$thead(tags$tr(
                  tags$th("Column"), tags$th("Description"), tags$th("Example")
                )),
                tags$tbody(
                  tags$tr(tags$td("1"), tags$td("Metabolite Name"), tags$td("metabolite 1")),
                  tags$tr(tags$td("2"), tags$td("X Coordinate"), tags$td("123")),
                  tags$tr(tags$td("3"), tags$td("Y Coordinate"), tags$td("134")),
                  tags$tr(tags$td("4"), tags$td("Intensity"), tags$td("2.0"))
                )
              )
            ),
            tags$div(
              class = "col-md-6",
              h4("Spatial Transcriptomics File"),
              tags$table(
                class = "table table-bordered",
                tags$thead(tags$tr(
                  tags$th("Column"), tags$th("Description"), tags$th("Example")
                )),
                tags$tbody(
                  tags$tr(tags$td("1"), tags$td("Gene ID"), tags$td("gene 1")),
                  tags$tr(tags$td("2"), tags$td("X Coordinate"), tags$td("123")),
                  tags$tr(tags$td("3"), tags$td("Y Coordinate"), tags$td("134")),
                  tags$tr(tags$td("4"), tags$td("MIDCount"), tags$td("2"))
                )
              )
              
            ),
            downloadButton("download_demometab", "Download demo spatial metabolomics data", class = "btn-sm"),
            downloadButton("download_demotrans", "Download demo spatial transcriptomics data", class = "btn-sm"),
            ),
          box(width = 4,
              style = "height: 320px; overflow-y: auto;",
            tags$img(src = "feature_matrix_txt.png", width = "100%", style = "max-width: 400px;"))),
          fluidRow(

          h2("Additional Resources"),
          tags$div(
            class = "list-group",
            a(href = "https://github.com/mzlab-research/SMIntegration.git", 
              class = "list-group-item",
              icon("github"), " Source Code & Local Installation Guide")
          ),

          h2("Technical Support"),
          p("For immediate assistance:"),
          tags$ul(
            tags$li(strong("Documentation:"), "Comprehensive guides available in Help section"),
            tags$li(strong("Email:"), tags$a(href = "mailto:denghaoke@genomics.cn", "denghaoke@genomics.cn"))
          ),
          hr(),
          p(em("Current Version: 1.0.0 | Last Updated: 2025-07-12"), style = "text-align: center;")
        )
)


