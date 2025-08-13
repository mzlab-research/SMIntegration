tabItem(tabName = "diff_select",
        fluidRow(
          h2("Step4: Differential Analysis"),
          h3("Region of Interest (ROI) Definition"),
          p("Define spatial regions for comparative analysis using three flexible approaches:"),
          tags$ol(
            tags$li(strong("Cluster-based:"), "Select domains from spatial clustering results"),
            tags$li(strong("Cell type-based:"), "Leverage annotated cell types for region definition"),
            tags$li(strong("Interactive selection:"), "Manually delineate regions using lasso or rectangular tools")
          ),
          p("After defining ROIs, click 'Add as treatment group' or 'Add as control group'. Use 'Clear selection' to reset. Finalize with 'Finish Selection' for downstream analysis.")
        ),        
        
        fluidRow(
          useShinyjs(),
          box(
            title = "ROI Selection Parameters",width = 6,
            style = "height: 400px; overflow-y: auto;",
            selectizeInput("plot_select", "ROI definition method:", choices = c( "Cluster-based"="Select clustering groups","Cell type-based"="Select cell groups","Interactive selection"="Use selection tool"), selected = "Select clustering groups"),
            selectInput("clusterdata_select", "Base visualization modality:", choices = c("Metabolite","Gene","Merge"),selected = "Gene"),
            uiOutput("ion_select_button_container"),
            h5("Define regions and assign to treatment and control groups:"),
            uiOutput("auto_add_button_container")
          ),
          box(
            title = "ROI Assignment Summary",width = 6,
            style = "height: 400px; overflow-y: auto;",
            tableOutput("selection_table"),
          )),
        fluidRow(
          box(
            width = 6,
            style = "height: 700px; overflow-y: auto;",
            uiOutput("scatter_plot_button_container")
          ),
          box(
            title ="Spatial Region Assignment",
            width = 6,
            style = "height: 660px; overflow-y: auto;",
            p("Final ROI assignment for differential analysis:"),
            tags$ul(
              tags$li(strong("Treatment:"), "Experimental regions (e.g., disease foci)"),
              tags$li(strong("Control:"), "Reference regions (e.g., healthy tissue)")
            ),
            plotOutput("group_spectrum",width = "500px", height = "500px"),
            downloadButton("download_samplelist", "Export ROI assignments"),
            downloadButton("download_group_spectrum_reactive", "Download image")
          )
        )
)#tabitem