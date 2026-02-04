# ==============================================================================
# U2_upload.R
# UI definition for the "Overall Distribution Analysis" and "File Upload" tab.
#
# Purpose:
#   Provides the interface for:
#   - Uploading spatial multi-omics data (Transcriptomics + Metabolomics).
#   - Supporting both raw text/CSV and Seurat (.rds) formats.
#   - Configuring spatial registration parameters (flip, rotate, translation).
#   - Displaying initial QC plots:
#     - Raw data preview.
#     - Registration comparison (before/after).
#     - Overall distribution heatmaps (Total Counts).
#     - Basic dataset statistics (dimensions, features).
#
# Interaction:
#   - Users select data source (Demo vs Upload).
#   - Toggles for "Perform Registration" enabling complex alignment UI.
#   - Submit button triggers data loading and preprocessing in S2_upload.R.
# ==============================================================================

#' @title File Upload & Registration UI
#' @description Defines the layout for data ingestion.
#' Includes conditional panels for standard upload vs. registration workflows.
#' Displays QC plots (Total Counts) and basic dataset statistics.
tabItem(tabName = "File_upload",
        
        fluidRow(
          h2("Step1: Overall Distribution Analysis"),
          h3("File upload"),
          conditionalPanel(
            condition = "input.do_registration == false",
            p("Upload spatially registered multi-omics datasets or use demo data for initial data check and visualization. Select parameters and click 'Submit' to begin.")
          ),
          conditionalPanel(
            condition = "input.do_registration == true",
            p("Upload raw multi-omics datasets (unregistered). Use the Registration Configuration panel to align them. Once aligned, downstream analysis begins automatically.")
          ),
          
          # Row 1: Upload and Registration Settings
          fluidRow(
            # Left: Upload Box
            box(
              title = "Upload Files or Try Demo Data",
              width = 6,
              p(strong("Data Requirement:"), "For standard analysis, input pre-registered data. For registration, input raw data."),
              selectizeInput("demo_select", "Select input data type", choices = c("Use demo data","Upload rds data","Upload txt data"), selected = "Use demo data"),
              
              # Registration Option
              checkboxInput("do_registration", "Perform Registration", value = FALSE),
              
              uiOutput("file_button_container"), 
              conditionalPanel(
                condition = "input.demo_select == 'Upload rds data' || input.demo_select == 'Upload txt data'",
                selectizeInput("speciesname_select", "Please select species name", 
                               choices = c("Homo sapiens" = "hsa", "Mus musculus" = "mmu"), 
                               selected = "mmu"),
                selectizeInput("metab_mode", "Select metabolic mode", 
                               choices = c("pos", "neg"), 
                               selected = "pos")
              ),
              p("Note: Processing time varies with data size. Please wait patiently."),
              
              # Submit Button (Only if NOT registering)
              conditionalPanel(
                condition = "input.do_registration == false",
                actionButton("Submit", "Submit")
              ),
              
              # Registration Instructions (If registering)
              conditionalPanel(
                condition = "input.do_registration == true",
                helpText("Step 1: Use demo data or upload data above."),
                helpText("Step 2: Click 'Preprocess Preview' on the right to align the two omics to approximate angles."),
                helpText("Step 3: Click 'Run Registration' on the right to align and start analysis.")
              )
            ),
            
            # Right: Registration Settings (Conditional)
            conditionalPanel(
              condition = "input.do_registration == true",
              box(
                title = "Registration Configuration",
                width = 6,
                status = "warning",
                solidHeader = TRUE,
                
                # Preprocessing Settings
                h4("Preprocessing"),
                checkboxInput("do_translate", "Apply Translation", value = TRUE),
                selectInput("flip_type", "Flip Type:", 
                            choices = c("No Flip" = "none", 
                                        "Vertical Flip" = "vertical", 
                                        "Horizontal Flip" = "horizontal")),
                selectInput("rotate_type", "Rotation Type:", 
                            choices = c("No Rotation" = "none",
                                        "90 Degree CCW" = "90cw",
                                        "180 Degree CCW" = "180cw", 
                                        "270 Degree CCW" = "270cw")),
                actionButton("preprocess_preview", "Preprocess Preview", class = "btn-info", width = "100%"),
                
                hr(),
                
                # Registration Settings
                h4("Registration"),
                selectInput("metab_display", "Metabolomics Display:",
                            choices = c("TotalCounts" = "TIC", "Select m/z Channel" = "specific")),
                uiOutput("mz_selector_reg"),
                checkboxInput("constrain_transform", "Constrain Transformation", value = FALSE),
                conditionalPanel(
                  condition = "input.constrain_transform",
                  numericInput("max_rotation", "Max Rotation Angle (degrees)", value = 10, min = 0, max = 45)
                ),
                actionButton("run_registration", "Run Registration", class = "btn-primary", width = "100%")
              )
            )
          )
        ),
        
        # Row 2: Registration Visualizations (Conditional)
        conditionalPanel(
          condition = "input.do_registration == true",
          
          # Row A: Raw Preview (8) + Preprocessing (4)
          fluidRow(
            box(
              title = "Raw Data Preview", status = "info", solidHeader = TRUE, width = 6,
              collapsible = TRUE, collapsed = FALSE,
              style = "height: 500px; overflow-y: auto;",
                column(6, h4("Original Transcriptomics"), plotOutput("pre_transcript_plot", height = "400px")),
                column(6, h4("Original Metabolomics"), plotOutput("pre_metabolic_plot", height = "400px"))
              
            ),
            box(
              title = "Preprocessing Results", status = "info", solidHeader = TRUE, width = 6,
              collapsible = TRUE, collapsed = FALSE,
              style = "height: 500px; overflow-y: auto;",
                h4("Preprocessed Metabolomics"),
                plotOutput("preprocessed_metabolic_plot", height = "400px")
              
            )
          ),
          fluidRow(
            box(
              title = "Registered Metabolomics", status = "info", solidHeader = TRUE, width = 6,
              collapsible = TRUE, collapsed = FALSE,
              plotOutput("aligned_metabolic", height = "400px")
            ),
            
            box(
              title = "RGB Overlay", status = "info", solidHeader = TRUE, width = 6,
              collapsible = TRUE, collapsed = FALSE,
              plotOutput("rgb_overlay", height = "400px")
            )
          ),
          fluidRow(
            box(
              title = "Registration Comparison Matrix", status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = FALSE,
              style = "height: 500px; overflow-y: auto;",
              column(6, h4("Comparison Matrix"),plotOutput("aligned_compare", height = "400px")),
              column(6, h4("Quality Metrics"),verbatimTextOutput("data_info_reg"),
                     downloadButton("download_reg_plots_zip", "Download All Registration Plots (ZIP)", class = "btn-success", style = "margin-top: 10px;"))
              
            )
          )
        ),
        
        fluidRow(
          h3(" Overall Distribution Analysis"),
          p("Spatial heatmaps below provide initial visualization for data quality assessment:")
        ),
        
        # Row 3: Final Data Display (Side by Side)
        fluidRow(

          box(
            title = "Overall distribution plot",
            p("This visualization performs initial data inspection through spatial intensity mapping:"),
            tags$ul(
              tags$li("Color gradient indicates molecular abundance (red: high, blue: low)"),
              tags$li("Left panel: total ion intensity in metabolomics"),
              tags$li("Right panel: total gene expression in transcriptomics")
            ),
            width = 6,
            style = "height: 600px; overflow-y: auto;",
            plotOutput("totalcounts_plot", width = "100%", height = "400px"),
              column(3, downloadButton("download_totalcounts_plot", "Download image")),
              column(3, downloadButton("download_totalcounts_data", "Download data"))
            
          ),
          box(
            title = "Basic information of data",
            style = "height: 600px; overflow-y: auto;",
            width = 6,
            tableOutput("basic_info")
          )
        )
)  # tabitem