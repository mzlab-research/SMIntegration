tabItem(tabName = "File_upload",

        
        fluidRow(
          h2("Step1:Overall Distribution Analysis"),
          h3("File upload"),
          p("Upload spatially registered multi-omics datasets or use demo data for initial data check and visualization. Select parameters and click 'Submit' to begin."),
          # file upload -----------------------
          box(
            title = "Upload Files or Try Demo Data",
            width = 6,
            p(strong("Data Requirement:"), "Input datasets must be pre-registered with consistent resolution. Refer to Tutorial for alignment protocols using SpatialData/ImageJ."),
            selectizeInput("demo_select", "Select input data type", choices = c("Use demo data","Upload rds data","Upload txt data"), selected = "Use demo data"),
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
            p("Note: Processing time varies with data size. Please wait patiently after submission."),
            actionButton("Submit", "Submit")
            
          ),
          box(
            title = "Basic information of data",
            style = "height: 210px; overflow-y: auto;",
            width = 6,
            tableOutput("basic_info")
          )  # box2

        ), 
        
        fluidRow(
          h3(" Overall Distribution Analysis"),
          p("Spatial heatmaps below provide initial visualization for data quality assessment:")
        ),
        fluidRow(
          box(
            title = "Overall distribution plot",
            p("This visualization performs initial data inspection through spatial intensity mapping:"),
            tags$ul(
              tags$li("Color gradient indicates molecular abundance (red: high, blue: low)"),
              tags$li("Left panel: total ion intensity in metabolomics"),
              tags$li("Right panel: total gene expression in transcriptomics")
            ),
            width = 8,
            style = "height: 600px; overflow-y: auto;",
            plotOutput("totalcounts_plot", width = "800px", height = "400px"),
            column(3,downloadButton("download_totalcounts_plot", "Download image")),
            column(3,downloadButton("download_totalcounts_data", "Download data"))
          )  
        )
)  # tabitem