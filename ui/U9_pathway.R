# ==============================================================================
# U9_pathway.R
# UI definition for the "Functional Association Analysis" tab (Step 5).
#
# Purpose:
#   Provides the interface for interpreting biological mechanisms via KEGG pathway analysis.
#   Integrates features from either Differential Analysis or Spatial Pattern Analysis.
#
# Key Features:
#   - Input Selection: Choose source features (Differential vs Pattern) and direction (Up/Down).
#   - Annotation Source: Use built-in KEGG database or upload custom annotation files.
#   - Pathway Overlap Analysis: Venn diagram showing pathways annotated by Genes vs Metabolites.
#   - Pathway Annotation Analysis: Dot plot ranking pathways by feature count.
#   - Pathway-Centric Integration:
#     - Detailed KEGG map visualization with nodes colored by feature status.
#     - Feature tables listing mapped genes/metabolites.
#     - Spatial distribution plots for specific pathway features.
#
# Structure:
#   - Parameter control box (Input type, Annotation source).
#   - Summary box.
#   - Row for Venn and Dot plots.
#   - Detailed pathway exploration section (Map image, Tables, Spatial plots).
# ==============================================================================

#' @title Functional Analysis UI
#' @description Defines the layout for pathway annotation analysis.
#' Includes controls for feature selection, annotation source, and multiple visualization levels (Global Overlap -> Ranked List -> Detailed Map -> Spatial Features).
tabItem(tabName = "Pathway_analysis",
        fluidRow(
          h2("Step5: Functional Association Analysis"),
          p("This module integrates differential or pattern-specific features to elucidate biological mechanisms through three complementary approaches:"),
          tags$ol(
            tags$li(strong("Pathway Overlap Analysis:"), "Count the number of pathways co-annotated by differentially expressed genes/metabolites or pattern-specific genes/metabolites, and those annotated solely by either genes or metabolites."),
            tags$li(strong("Pathway Annotation Analysis:"), "Visualize pathways with the highest number of annotated differential/pattern-specific features."),
            tags$li(strong("Pathway-Centric Feature Integration:"), "Integrate pathway topology with spatial distributions of key molecules")
          ),
          p("By analyzing pathways co-annotated by spatially significant features, we uncover critical biological mechanisms underlying tissue phenotypes. The platform supports two feature input types:"),
          tags$ul(
            tags$li(strong("Differential features:"), "Top 300 DEGs/DAMs (by adjusted p-value) from differential analysis"),
            tags$li(strong("Pattern-specific features:"), "Top 300 features per spatial module (by correlation score)")
          ),
          p("Default KEGG annotation is provided. For custom pathways, select 'Upload annotation data'.")
        ),
        fluidRow(
          box(
            title = "Analysis Parameters",
            width = 6,
            uiOutput("functional_data_button_container"),
            uiOutput("functional_diffdata_button_container"),
            selectizeInput("annotation_select", "Annotation source", choices = c("Built-in KEGG database"="use built-in database","Custom annotation upload"="upload your annotation data"), selected = "use built-in database"),
            uiOutput("annotationfile_button_container"), 
            p("Note: Pathway mapping scales with feature count. Please avoid duplicate submissions."),
            actionButton("start_annotation", "Launch Functional Analysis")
          ),
          box(
            title = "Annotation Summary",
            width = 6,
            style = "height: 250px; overflow-y: auto;",
            column(3,downloadButton("download_annotation_data_all", "Export Results")),
            tableOutput("annotation_info")
          )
        ),
        fluidRow(
          column(6, 
                 h3("Pathway Overlap Analysis"),
                 p("Venn diagram showing shared and modality-specific KEGG pathway annotations.")
          ),
          column(6,
                 h3("Pathway Annotation Analysis"),
                 p("Dot plot visualizing pathways with the most mapped features (Top 20 by Count).")
          )
        ),
        fluidRow(
          box(
            title = "Pathway Annotation Overlap",width = 6,
            style = "height: 580px; overflow-y: auto;",
            p("Intersection of pathways annotated by input features:"),
            tags$ul(
              tags$li("Left circle: Gene-annotated pathways"),
              tags$li("Right circle: Metabolite-annotated pathways"),
              tags$li("Overlap: Co-annotated pathways")
            ),
            plotOutput("venn_plot", width = "400px", height = "400px"),
            column(3,downloadButton("download_venn_plot", "Download image"))
          ),
          box(
            title = "Pathway Annotation Count",width = 6,
            style = "height: 580px; overflow-y: auto;",
            selectInput("bubble_pathway_types", "Filter by pathway type:", choices = NULL,selected = NULL),
            plotOutput("dotplot_run", width = "600px", height = "400px"),
            column(3,downloadButton("download_dotplot", "Download image")),
            column(3,downloadButton("download_dotplot_data", "Export data"))
          )
        ),
        column(12, br()),
        h3("Pathway-Centric Feature Integration"),
        p("Co-annotation network mapping of genes and metabolites within biological pathways. Select pathways below to explore:"),
        fluidRow(
          box(
            title = "Pathway Topology Mapping",width = 12,
            style = "height: 800px; overflow-y: auto;",
            column(6,selectInput("Pathway_types", "Annotation type:", choices = c("Co-annotated pathways"="gm",
                                                                         "Gene-only pathways"="onlyg" ,"Metabolite-only pathways"="onlym"),selected = "gm")),
                   column(6,selectInput("Pathway_select", "Select pathway:", choices = NULL,selected = NULL)),
            column(3,downloadButton("download_pathway_annotation_plot", "Download pathway")),
            column(3,downloadButton("download_pathway_annotation_data", "Export features")),
            imageOutput("pathway_annotation_plot", width = "100%", height = "400px")
          ),
          column(12, br()),
          box(
            title = "Pathway Feature Table: Metabolites",width = 6,
            style = "height: 350px; overflow-y: auto;",
            p("Key metabolites in selected pathway:"),
            tableOutput("annotation_namemap_info_m"),
          ),
          box(
            title = "Pathway Feature Table: Gene",width = 6,
            style = "height: 350px; overflow-y: auto;",
            p("Key genes in selected pathway:"),
            tableOutput("annotation_namemap_info_t")
          ),
          h3("Spatial Activity Visualization"),
          p("Examine spatial distributions of key pathway components:"),
          box(
            title = "Metabolite Spatial Distribution",width = 6,
            style = "height: 550px; overflow-y: auto;",
            selectizeInput("Pathway_m_select", "Select pathway metabolite:", choices = NULL,selected = NULL,options = list(server = TRUE)),
            plotOutput("pathway_annotation_m_plot", width = "400px", height = "400px"),
            column(3,downloadButton("download_pathway_annotation_m_plot", "Download image")),
            column(3,downloadButton("download_pathway_annotation_m_data", "Export data"))
          ),
          box(
            title = "Gene Spatial Distribution",width = 6,
            style = "height: 550px; overflow-y: auto;",
            selectizeInput("Pathway_t_select", "Select pathway gene:", choices = NULL,selected = NULL,options = list(server = TRUE)),
            plotOutput("pathway_annotation_t_plot", width = "400px", height = "400px"),
            column(3,downloadButton("download_pathway_annotation_t_plot", "Download image")),
            column(3,downloadButton("download_pathway_annotation_t_data", "Export data"))
          )
        )
)#tabitem