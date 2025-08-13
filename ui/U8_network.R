tabItem(tabName = "con_net",
        fluidRow(
          h2("Step4: Differential Analysis"),
          h3("Group-Specific Correlation Networks"),
          p("Reveal spatial co-expression relationships between top differential features:"),
          tags$ul(
            tags$li("Constructs separate networks for each group"),
            tags$li("Computes Spearman correlations using group-specific pixels"),
            tags$li("Multiple testing correction: Applies Bonferroni adjustment to p-values"),
            tags$li("Filters: |r| > 0.6 & adjusted p-value < 0.01")
          ),
          p("By default, top 20 differential features per modality (by p-value) are included.")
        ),
        fluidRow(
          box(
            width = 6,
            numericInput("cor_Coefficient1","Correlation threshold (|r|):",min=0,max=1,value=0.6,step=0.1),
            numericInput("cor_p_value1","adjusted p-value threshold:",min=0.001,max=1,value=0.01,step=0.01),
            numericInput("topdiffnum","Features per modality:",min=10,max=100,value=20,step=1)
          ),
          box(
            width = 6,
            selectizeInput("layout_type", "Network layout:", choices = c("Force-directed (fr)"="fr","Kamada-Kawai (kk)"="kk","Automatic (nicely)"="nicely"), selected = "kk"),
            selectizeInput("show_node_name", "Display node labels:", choices = c("yes"=TRUE,"no"=FALSE), selected = TRUE),
            numericInput("node_size","Node size:",min=1,max=10,value=6,step=0.5),
            numericInput("node_name_size","Label size:",min=1,max=8,value=3,step=0.1)
          )
        ),
        fluidRow(
          box(
            title = "Treatment Group Network",
            tags$ul(
              tags$li(strong("Edges:"), "Red: positive correlation, Green: negative correlation"),
              tags$li(strong("Nodes Shapes:"), "Triangles: metabolites, Circles: genes"),
              tags$li(strong("Node Colors:"), "Feature abundance gradient (blue-low to red-high)")
            ),
            width = 6,
            style = "height: 900px; overflow-y: auto;",
            plotOutput("cornetwork_plot1", width = "700px", height = "700px"),
            column(3,downloadButton("download_cornetwork_plot1", "Download image")),
            column(3,downloadButton("download_cornetwork_data1", "Export network"))
          ),
          box(
            title = "Control Group Network",
            tags$ul(
              tags$li("Identical node layout for direct comparison"),
              tags$li("Edge differences indicate altered regulatory relationships"),
              tags$li("Node color differences reflect abundance changes")
            ),
            width = 6,
            style = "height: 900px; overflow-y: auto;",
            plotOutput("cornetwork_plot2", width = "700px", height = "700px"),
            column(3,downloadButton("download_cornetwork_plot2", "Download image")),
            column(3,downloadButton("download_cornetwork_data2", "Export network"))
          )
        )
)#tabitem