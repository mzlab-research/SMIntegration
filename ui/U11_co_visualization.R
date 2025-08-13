tabItem(tabName = "co_visual",
        fluidRow(
          h2("Multi-Molecule Spatial Co-visualization")
        ),
        fluidRow(
          p("Integrated visualization of molecular co-localization patterns using RGB channel mapping:"),
          box(
            title = "Molecular Feature Selection",
            p("Step 1: Select up to three molecular features for co-visualization:"),
            tags$ul(
              tags$li("Feature 1 assigned to red channel"),
              tags$li("Feature 2 assigned to green channel"),
              tags$li("Feature 3 assigned to blue channel")
            ),
            width = 6,
            selectizeInput("cor_F1_select", "Feature 1 (Red):", choices = NULL,selected = NULL,options = list(server = TRUE)),
            selectizeInput("cor_F2_select", "Feature 2 (Green):", choices = NULL,selected = NULL,options = list(server = TRUE)),
            selectizeInput("cor_F3_select", "Feature 3 (Blue):", choices = NULL,selected = NULL,options = list(server = TRUE)),
            p("Note: Processing time depends on image resolution and feature count."),
            actionButton("start_Co_visualisation_analysis", "Generate Co-visualization")
          ),
          column(12, br()),
          box(
            title="RGB Spatial Overlay",
            width = 6,
            style = "height: 600px; overflow-y: auto;",
            p("Step 2: Pseudocolor representation of molecular co-distribution:"),
            tags$ul(
              tags$li("Red: Feature 1 intensity"),
              tags$li("Green: Feature 2 intensity"),
              tags$li("Blue: Feature 3 intensity")
            ),
            plotOutput("pseudocolor_plot", width = "400px", height = "400px"),
            column(3,downloadButton("download_pseudocolor_plot", "Download image")),
            column(3,downloadButton("download_pseudocolor_data", "Export data"))
          ),
          box(
            title = "Co-expression Pattern Analysis",
            p("Step 3: Quantitative assessment of combinatorial expression states:"),
            tags$ul(
              tags$li("Expression quartiles: High (>75%), Medium, Low (<25%)"),
              tags$li("8 possible combinatorial states for 2 features"),
              tags$li("27 possible combinatorial states for 3 features"),
              tags$li("Dominant patterns reveal functional co-regulation")
            ),
            width = 6,
            style = "height: 600px; overflow-y: auto;",
            plotOutput("cor_gene_metabolite_Co_visualisation_plot", width = "400px", height = "400px"),
            column(3,downloadButton("download_Co_visualisation_plot", "Download image")),
            column(3,downloadButton("download_Co_visualisation_data", "Export data"))
          )
        )
)#tabitem