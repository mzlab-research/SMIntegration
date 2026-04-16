# ==============================================================================
# S2_upload.R
# Server logic for Data Upload, Preprocessing, and Spatial Registration.
#
# Purpose:
#   Handles the ingestion of user data, performing validity checks, 
#   preprocessing (creation of Seurat objects), and optional spatial registration.
#
# Inputs:
#   - input$transfile, input$metabfile: User uploaded files.
#   - input$demo_select: Selection of demo data.
#   - Registration parameters (rotation, flip, translation).
#
# Logic:
#   1. Data Loading: Reads TXT or RDS files using 'fread' or 'readRDS'.
#   2. Validation: Checks for required columns (x, y, gene/metabolite, intensity).
#   3. Registration (Optional):
#      - Uses 'RNiftyReg' for affine/rigid transformation.
#      - Aligns Metabolomics (source) to Transcriptomics (target) coordinates.
#      - visualizes overlay and calculates MSE/Correlation metrics.
#   4. Seurat Object Creation: Converts raw data into Seurat objects ('CreateSeuratObject').
#   5. QC Visualization: Generates 'totalcounts_plot' and basic info tables.
#
# Outputs:
#   - Reactive objects: 'transdata', 'metabdata', 'data_rds' (list of Seurat objects).
#   - Plots: 'aligned_metabolic', 'rgb_overlay', 'totalcounts_plot'.
# ==============================================================================

# Reactive values to hold state across modules
diff_omics <- reactiveVal(NULL)
spatial_pattern_ionlist <- reactiveVal(NULL)

#' @title Load Transcriptomics Data
#' @description Reads and validates transcriptomics data from user upload (txt/rds) or demo files.
#' @return A dataframe or Seurat object containing transcriptomics data (x, y, counts).
transdata<-eventReactive(c(input$Submit), {
  source("./source/preprocessing/runxy.R")
  
  withProgress(message = "Loading data ...",value=0.8,{

  if(input$demo_select == "Use demo data"){
    data<-demotrans_rds()
  }else if(input$demo_select == "Upload txt data"){
    if(!is.null(input$transfile) && input$transfile$name != ""){
      data <- fread(input$transfile$datapath)
      if ("geneID" %in% colnames(data) & "x" %in% colnames(data) & "y" %in% colnames(data) & "MIDCount" %in% colnames(data)) {
        data <- data[, c("geneID", "x", "y", "MIDCount")]
      } else {
        showNotification("Required columns not found in the data.")
      }
      shiny::validate(
        shiny::need(is.numeric(data$x), "x column must be numeric"),
        shiny::need(is.numeric(data$y), "y column must be numeric")
      )
      data<-colname_change(data)
    }
  }else if(input$demo_select == "Upload rds data"){
    if(!is.null(input$transfile) && input$transfile$name != ""){
      data <- readRDS(input$transfile$datapath)
      has_counts <- !is.null(data@assays$Spatial$counts)
      has_data <- !is.null(data@assays$Spatial$data)
      shiny::validate(
        need(data@active.assay == "Spatial" && "Spatial" %in% names(data@assays),
             "Active assay must be named 'Spatial'!"),
        need(c("x", "y") %in% colnames(data@meta.data),
             "The meta.data of Seurat object must contain two columns named 'x' and 'y'."),
        need(is.numeric(data@meta.data$x), "x column must be numeric"),
        need(is.numeric(data@meta.data$y), "y column must be numeric"),
        
        need("Spatial" %in% names(data@assays), 
             "Spatial assay not found in the Seurat object."),
        need(!is.null(data@assays$Spatial$counts), 
             "Counts matrix in Spatial assay is NULL."),
        need(is.matrix(data@assays$Spatial$counts) || inherits(data@assays$Spatial$counts, "dgCMatrix"),
             "Counts must be a dense matrix or a dgCMatrix (sparse)."),
        need(ncol(data@assays$Spatial$counts) > 0, 
             "Counts matrix has no columns (no spots)."),
        need(
          (has_counts && all(grepl("_", colnames(data@assays$Spatial$counts)))) ||
            (has_data && all(grepl("_", colnames(data@assays$Spatial$data)))),
             "Column names of the counts matrix must contain '_' to separate x and y coordinates (e.g., sample:10_20).")
      )
    }
  }
    return(data)
  }) 
}) 

#' @title Load Metabolomics Data
#' @description Reads and validates metabolomics data from user upload (txt/rds) or demo files.
#' Handles switching between raw data (if registration is needed) and processed data.
#' @return A dataframe or Seurat object containing metabolomics data.
metabdata<-eventReactive(c(input$Submit), {
  source("./source/preprocessing/runxy.R")
  withProgress(message = "Loading data ...",value=0.8,{

  if(input$demo_select == "Use demo data"){
    # Conditional loading based on do_registration
    if(input$do_registration){ # if registration is ON, load the raw demo
      data<-noreg_demometab_rds()

    } else { # if registration is OFF, load the processed demo
      data<-demometab_rds()
    }
  }else if(input$demo_select == "Upload txt data"){
    if(!is.null(input$metabfile) && input$metabfile$name != ""){
      data <- fread(input$metabfile$datapath)
      if ("metabolite" %in% colnames(data) & "x" %in% colnames(data) & "y" %in% colnames(data) & "Intensity" %in% colnames(data)) {
        data <- data[, c("metabolite", "x", "y", "Intensity")]
      } else if ("mz" %in% colnames(data) & "x" %in% colnames(data) & "y" %in% colnames(data) & "Intensity" %in% colnames(data)) {
        data <- data[, c("mz", "x", "y", "Intensity")]
      } else {
        showNotification("Required columns not found in the data.")
      }
      shiny::validate(
        shiny::need(is.numeric(data$x), "x column must be numeric"),
        shiny::need(is.numeric(data$y), "y column must be numeric")
      )
      data<-colname_change(data)
    }
    }else if(input$demo_select == "Upload rds data"){
      if(!is.null(input$metabfile) && input$metabfile$name != ""){
        data <- readRDS(input$metabfile$datapath)
        has_counts <- !is.null(data@assays$Spatial$counts)
        has_data <- !is.null(data@assays$Spatial$data)
        shiny::validate(
          need(data@active.assay == "Spatial" && "Spatial" %in% names(data@assays),
               "Active assay must be named 'Spatial'!"),
          need(c("x", "y") %in% colnames(data@meta.data),
               "The meta.data of Seurat object must contain two columns named 'x' and 'y'."),
          need(is.numeric(data@meta.data$x), "x column must be numeric"),
          need(is.numeric(data@meta.data$y), "y column must be numeric"),
          
          need("Spatial" %in% names(data@assays), 
               "Spatial assay not found in the Seurat object."),
          need(!is.null(data@assays$Spatial$counts), 
               "Counts matrix in Spatial assay is NULL."),
          need(is.matrix(data@assays$Spatial$counts) || inherits(data@assays$Spatial$counts, "dgCMatrix"),
               "Counts must be a dense matrix or a dgCMatrix (sparse)."),
          need(ncol(data@assays$Spatial$counts) > 0, 
               "Counts matrix has no columns (no spots)."),
          need(
            (has_counts && all(grepl("_", colnames(data@assays$Spatial$counts)))) ||
              (has_data && all(grepl("_", colnames(data@assays$Spatial$data)))),
               "Column names of the counts matrix must contain '_' to separate x and y coordinates (e.g., sample:10_20).")
        )
      }
  }

  return(data)
  })
})

#' @title Download Total Counts Plot
#' @description Saves the Total Counts QC plot as a PNG file.
#' @return A PNG file named "totalcounts_plot.png".
output$download_totalcounts_plot <- downloadHandler(
  filename = function() {
    "totalcounts_plot.png"
  },
  content = function(file) {

    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- get("totalcounts_plot", envir = .GlobalEnv)

    ggsave(file, plot = p, device = "png", width = 12, height = 8,bg = "#FFFFFF", dpi = 300)
    })
  })

#' @title Download Total Counts Data
#' @description Exports the Total Counts data for both Metabolomics and Transcriptomics.
#' @return A ZIP file containing two CSV files: "Metabolite_data.csv" and "Gene_data.csv".
output$download_totalcounts_data <- downloadHandler(
  filename = function() {
    "totalcounts_data.zip"
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
    data_rds<-data_rds()
    data_mrds=data_rds[[1]]
    data_trds=data_rds[[2]]

    data_m<-data.frame(x=data_mrds$x,y=data_mrds$y,nCount_Spatial=data_mrds$nCount_Spatial)
    data_t<-data.frame(x=data_trds$x,y=data_trds$y,nCount_Spatial=data_trds$nCount_Spatial)

    datalist<-list(data_m,data_t)
    tempdir <- setwd(tempdir())
    on.exit(setwd(tempdir))
    fi=c("Metabolite_data.csv","Gene_data.csv")
    for (i in 1:length(fi)) {
      data<-datalist[[i]]
      write.csv(data, fi[i], row.names = FALSE)
    }
    
    zip(file,fi)
    })
  })




#' @title Preprocess Metabolomics Data List
#' @description Generates a list of preprocessed metabolomics data.
#' If registration is enabled, it returns NULL (waiting for registration).
#' Otherwise, it processes data from TXT or RDS inputs.
#' @return A list containing the processed metabolomics Seurat object.
pre_metabolomicslist<-eventReactive(c(input$Submit), {
  # If registration is selected, we don't produce the "analysis-ready" list yet.
  # We wait for the user to run registration.
  if (input$do_registration) {
    return(NULL)
  }
  
  if(input$demo_select == "Upload txt data" ){
    if (is.null(transdata()) || is.null(metabdata())) {
      showNotification("Please ensure that both files have been uploaded.")
    } else {
      withProgress(message = "Preprocessing metabolomics data ...",value=0.8,{
        pre_metabdata <- metabdata()
        binsize<-1
        source("./source/preprocessing/runxy.R")
        metabolomicslist<-runxy(input=pre_metabdata,binsize=binsize)
        return(metabolomicslist)
      })
    }
  }else if(input$demo_select == "Upload rds data" || input$demo_select == "Use demo data"){
    if (is.null(transdata()) || is.null(metabdata())) {
      showNotification("Please ensure that three files have been uploaded.")
    } else {
      withProgress(message = "Preprocessing metabolomics data ...",value=0.8,{
        pre_metabdata <- metabdata()
        source("./source/preprocessing/run_prerds.R")
        metabolomicslist<-run_prerds(data=pre_metabdata)
        return(metabolomicslist)
      })
    }
  }
})

# Switcher for downstream analysis
#' @title Finalize Metabolomics Data for Downstream
#' @description Switches between registered data and standard preprocessed data based on user choice.
#' If registration is performed, it uses the registered output.
#' Otherwise, it aligns metabolomics and transcriptomics based on shared coordinates (intersection).
#' @return A Seurat object ready for downstream analysis.
final_metabolomics_list <- reactive({
  if (input$do_registration) {
    req(metab_reg())
    return(metab_reg())
  } else {
    req(pre_metabolomicslist())
    req(pre_genelist())
    
    pre_metabolomicslist=pre_metabolomicslist()
    pre_genelist=pre_genelist()
    m= as.data.frame(t(pre_metabolomicslist@assays$Spatial$counts)) %>%
      mutate(X=pre_metabolomicslist$x,Y=pre_metabolomicslist$y) %>%
      dplyr::select(X,Y,everything())
    g=data.frame(X=pre_genelist$x,Y=pre_genelist$y)
    source("./source/preprocessing/runxy.R")
    m$X_Y <- paste0("sample:",m$X, "_", m$Y)
    g$X_Y <- paste0("sample:",g$X, "_", g$Y)
    intersect_xy <- intersect(g$X_Y, m$X_Y)
    m <- m %>%
      dplyr::filter(X_Y %in% intersect_xy) %>%
      dplyr::arrange(match(X_Y, intersect_xy)) 
    m$X_Y <- NULL
    m[is.na(m)] <- 0
    mz_mat <- create_sparse_from_wide(m) 
    coords <- data.frame(x=m$X,y=m$Y)
    metab <- create_seurat(mz_mat, coords, binsize=1)
    metab$x_y=paste0("x",metab$x,"_y",metab$y)
    
    return(metab)
    
  }
  
})

#' @title Preprocess Transcriptomics Data List
#' @description Generates a list of preprocessed transcriptomics data.
#' If registration is disabled, it waits for metabolomics data to ensure synchronization.
#' @return A list containing the processed transcriptomics Seurat object.
pre_genelist<-eventReactive(c(input$Submit), {
  # If doing registration, we don't wait for pre_metabolomicslist (which is NULL)
  if (!input$do_registration) {
    req(pre_metabolomicslist())
  }
  
  withProgress(message = "Preprocessing transcriptomics data ...",value=0.8,{
    if(input$demo_select == "Upload txt data" ){
      pre_transdata <- transdata()
      binsize<-1
      source("./source/preprocessing/runxy.R")
      genelist<-runxy(input=pre_transdata,binsize=binsize)
      return(genelist)
    }else if(input$demo_select == "Upload rds data" || input$demo_select == "Use demo data"){
      pre_transdata <- transdata()
      source("./source/preprocessing/run_prerds.R")
      genelist<-run_prerds(data=pre_transdata)
      return(genelist)
    }
  })
})
###############################
# Reactive values to store preprocessing state
preprocess_done <- reactiveVal(FALSE)
preprocessed_data <- reactiveVal(NULL)

# Read Transcript Data (Raw for Registration)
#' @title Load Raw Transcriptomics for Registration
#' @description Reads transcriptomics data (TXT, RDS, or Demo) specifically for the registration module.
#' Validates input format and aggregates counts if necessary.
#' @return A data frame with columns: X, Y, TotalCounts.
transcript_data_reg <- reactive({
  
  # Determine source data based on selection
  if(input$demo_select == "Use demo data"){
    data <- transdata()
    req(data)
  } else {
    # Upload mode - Read directly to allow preview before Submit
    req(input$transfile)
    source("./source/preprocessing/runxy.R")
    
    if(input$demo_select == "Upload txt data"){
      data <- fread(input$transfile$datapath)
      if ("geneID" %in% colnames(data) & "x" %in% colnames(data) & "y" %in% colnames(data) & "MIDCount" %in% colnames(data)) {
        data <- data[, c("geneID", "x", "y", "MIDCount")]
      }
      shiny::validate(
        shiny::need(is.numeric(data$x), "x column must be numeric"),
        shiny::need(is.numeric(data$y), "y column must be numeric")
      )
      data <- colname_change(data)
      
    } else if(input$demo_select == "Upload rds data"){
      data <- readRDS(input$transfile$datapath)
      shiny::validate(
        shiny::need(data@active.assay == "Spatial" && "Spatial" %in% names(data@assays), "Active assay must be named 'Spatial'!"),
        shiny::need(all(c("x","y") %in% colnames(data@meta.data)), "The meta.data of Seurat object must contain 'x' and 'y'."),
        shiny::need(is.numeric(data@meta.data$x), "x column must be numeric"),
        shiny::need(is.numeric(data@meta.data$y), "y column must be numeric")
      )
    }
  }
  
  withProgress(message = "Loading Transcriptomics for Registration...", value = 0.5, {
    # Check format (Seurat vs Data.Frame)
    if (inherits(data, "Seurat")) {
      meta <- data@meta.data
      if ("nCount_Spatial" %in% colnames(meta)) {
        counts <- meta$nCount_Spatial
      } else {
        counts <- colSums(GetAssayData(data, assay = "Spatial", layer = "counts"))
      }
      df <- data.frame(X = meta$x, Y = meta$y, TotalCounts = counts)
    } else {
      # Assuming data frame with columns: geneID, x, y, MIDCount (from runxy check)
      # We need to aggregate to spots
      if("MIDCount" %in% colnames(data)) {
        df <- data %>% dplyr::group_by(x, y) %>% dplyr::summarise(TotalCounts = sum(MIDCount), .groups='drop')
        colnames(df) <- c("X", "Y", "TotalCounts")
      } else {
        # Maybe already processed? Just check X, Y
        if("x" %in% colnames(data)) colnames(data)[colnames(data)=="x"] <- "X"
        if("y" %in% colnames(data)) colnames(data)[colnames(data)=="y"] <- "Y"
        df <- data
      }
    }
    
    return(df)
  })
})

# Read Metabolic Data (Raw for Registration)
#' @title Load Raw Metabolomics for Registration
#' @description Reads metabolomics data (TXT, RDS, or Demo) specifically for the registration module.
#' Handles conversion from long to wide format if required.
#' @return A wide-format data frame with columns X, Y, and one column per metabolite/mz.
metabolic_data_reg <- reactive({
  
  # Determine source data based on selection
  if(input$demo_select == "Use demo data"){
    data <- metabdata()
    req(data)
  } else {
    # Upload mode - Read directly to allow preview before Submit
    req(input$metabfile)
    source("./source/preprocessing/runxy.R")
    
    if(input$demo_select == "Upload txt data"){
      data <- fread(input$metabfile$datapath)
      if ("metabolite" %in% colnames(data) & "x" %in% colnames(data) & "y" %in% colnames(data) & "Intensity" %in% colnames(data)) {
        data <- data[, c("metabolite", "x", "y", "Intensity")]
      } else if ("mz" %in% colnames(data) & "x" %in% colnames(data) & "y" %in% colnames(data) & "Intensity" %in% colnames(data)) {
        data <- data[, c("mz", "x", "y", "Intensity")]
      }
      shiny::validate(
        shiny::need(is.numeric(data$x), "x column must be numeric"),
        shiny::need(is.numeric(data$y), "y column must be numeric")
      )
      data <- colname_change(data)
      
    } else if(input$demo_select == "Upload rds data"){
      data <- readRDS(input$metabfile$datapath)
      shiny::validate(
        shiny::need(data@active.assay == "Spatial" && "Spatial" %in% names(data@assays), "Active assay must be named 'Spatial'!"),
        shiny::need(c("x","y") %in% colnames(data@meta.data), "The meta.data of Seurat object must contain 'x' and 'y'."),
        shiny::need(is.numeric(data@meta.data$x), "x column must be numeric"),
        shiny::need(is.numeric(data@meta.data$y), "y column must be numeric")
      )
    }
  }
  
  withProgress(message = "Loading Metabolomics for Registration...", value = 0.5, {
    if (inherits(data, "Seurat")) {
      m <- as.data.frame(t(as.matrix(GetAssayData(data, assay = "Spatial", layer = "counts"))))
      m$X <- data@meta.data$x
      m$Y <- data@meta.data$y
      df <- m %>% dplyr::select(X, Y, everything())
    } else {
      # Data frame
      # Check format: metabolite, x, y, Intensity (Long) OR mz, x, y, Intensity (Long)
      if (("metabolite" %in% colnames(data) || "mz" %in% colnames(data)) && "Intensity" %in% colnames(data)) {
        # Convert Long to Wide
        val_col <- if("metabolite" %in% colnames(data)) "metabolite" else "mz"
        df <- reshape2::dcast(data, x + y ~ get(val_col), value.var="Intensity", fill=0)
        colnames(df)[1:2] <- c("X", "Y")
      } else {
        # Assume Wide
        if("x" %in% colnames(data)) colnames(data)[colnames(data)=="x"] <- "X"
        if("y" %in% colnames(data)) colnames(data)[colnames(data)=="y"] <- "Y"
        df <- data
      }
    }
    
    return(df)
  })
})

# Get m/z columns
#' @title Extract m/z Column Names
#' @description Identifies metabolite intensity columns from the loaded metabolomics data.
#' Excludes spatial coordinate columns (X, Y).
#' @return A character vector of column names.
mz_columns_reg <-  eventReactive(c(input$Submit), {
  data <- metabolic_data_reg()
  req(data)
  setdiff(colnames(data), c("X", "Y"))
})

# Dynamic UI for m/z selection
#' @title Render m/z Selector UI
#' @description Dynamically creates a dropdown menu for selecting a specific metabolite/m/z 
#' to visualize if "Specific m/z" display mode is chosen.
output$mz_selector_reg <- renderUI({
  if(input$metab_display == "specific") {
    selectInput("selected_mz_reg", "Select m/z Channel:", choices = mz_columns_reg())
  }
})

# Preprocessing Function
#' @title Preprocess Metabolic Coordinates
#' @description Applies geometric transformations (translation, flip, rotation) to metabolic data
#' to provide an initial coarse alignment with transcriptomics data.
#' @param metabolic_data Data frame of metabolomics data.
#' @param transcript_data Data frame of transcriptomics data (reference).
#' @param do_translate Boolean, whether to align minimum coordinates.
#' @param flip_type String, "vertical", "horizontal", or "none".
#' @param rotate_type String, "90cw", "180cw", "270cw", or "none".
#' @return Transformed metabolomics data frame.
preprocess_metabolic_data <- function(metabolic_data, transcript_data, do_translate, flip_type, rotate_type) {
  metabolic_plot_data <- metabolic_data
  # Check if transcript_data has TotalCounts, otherwise use first feature or dummy
  if("TotalCounts" %in% colnames(transcript_data)){
    val <- transcript_data$TotalCounts
  } else {
    val <- rowSums(transcript_data[, 3:ncol(transcript_data)])
  }
  
  transcript_plot_data <- data.frame(
    X = transcript_data$X,
    Y = transcript_data$Y,
    expression_value = val
  )
  
  # Flip
  if(flip_type == "vertical") {
    max_y <- max(metabolic_plot_data$Y)
    metabolic_plot_data$Y <- max_y - metabolic_plot_data$Y
  } else if(flip_type == "horizontal") {
    max_x <- max(metabolic_plot_data$X)
    metabolic_plot_data$X <- max_x - metabolic_plot_data$X
  }
  
  # Rotate
  if(rotate_type != "none") {
    if(rotate_type == "90cw") {
      metabolic_plot_data <- metabolic_plot_data |>
        mutate(X1 = max(Y) - Y, Y1 = X) |>
        mutate(X = X1, Y = Y1) |>
        dplyr::select(-X1, -Y1)
    } else if(rotate_type == "180cw") {
      metabolic_plot_data <- metabolic_plot_data |>
        mutate(X1 = max(X) - X, Y1 = max(Y) - Y) |>
        mutate(X = X1, Y = Y1) |>
        dplyr::select(-X1, -Y1)
    } else if(rotate_type == "270cw") {
      metabolic_plot_data <- metabolic_plot_data |>
        mutate(X1 = Y, Y1 = max(X) - X) |>
        mutate(X = X1, Y = Y1) |>
        dplyr::select(-X1, -Y1)
    }
  }
  
  # Translate
  if(do_translate) {
    original_metabolic_coords <- metabolic_plot_data[, c("X", "Y")]
    original_transcript_coords <- transcript_plot_data[, c("X", "Y")]
    xchange <- min(original_metabolic_coords$X) - min(original_transcript_coords$X)
    ychange <- min(original_metabolic_coords$Y) - min(original_transcript_coords$Y)
    
    metabolic_plot_data <- metabolic_plot_data |>
      dplyr::mutate(X = X - xchange, Y = Y - ychange)
  }
  return(metabolic_plot_data)
}

# Preprocess Event
#' @title Trigger Preprocessing Preview
#' @description Reacts to the "Preview" button. Runs the `preprocess_metabolic_data` function
#' with current parameters (flip, rotate, translate) and stores the result for visualization.
observeEvent(input$preprocess_preview, {
  req(metabolic_data_reg(), transcript_data_reg())
  
  withProgress(message = "Preprocessing Preview...", value = 0.5, {
    processed_data <- preprocess_metabolic_data(
      metabolic_data = metabolic_data_reg(),
      transcript_data = transcript_data_reg(),
      do_translate = input$do_translate,
      flip_type = input$flip_type,
      rotate_type = input$rotate_type
    )
    
    preprocessed_data(processed_data)
    preprocess_done(TRUE)
  })
})

# Plot: Original Transcript
output$pre_transcript_plot <- renderPlot({
  data <- transcript_data_reg()
  req(data)
  
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  ggplot(data, aes(x = X, y = Y, color = TotalCounts)) +
    geom_point(size = 1) +
    scale_color_gradientn(colours = heatmap_Palette(100), name = "Expression") +
    labs(title = "Original Transcriptomics", x = "X", y = "Y") +
    coord_equal() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
})

# Plot: Original Metabolic
output$pre_metabolic_plot <- renderPlot({
  data <- metabolic_data_reg()
  req(data)
  mz_cols <- mz_columns_reg()
  
  if(input$metab_display == "TIC") {
    plot_data <- data
    plot_data$Intensity <- rowSums(data[, mz_cols, drop = FALSE])
    title <- "Original Metabolomics"
  } else {
    req(input$selected_mz_reg)
    plot_data <- data
    plot_data$Intensity <- data[[input$selected_mz_reg]]
    title <- paste("Original Metabolomics (m/z:", input$selected_mz_reg, ")")
  }
  
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  ggplot(plot_data, aes(x = X, y = Y, color = Intensity)) +
    geom_point(size = 1) +
    scale_color_gradientn(colours = heatmap_Palette(100), name = "Intensity") +
    labs(title = title, x = "X", y = "Y") +
    coord_equal() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
})

# Plot: Preprocessed Metabolic
output$preprocessed_metabolic_plot <- renderPlot({
  if(!preprocess_done()) return(NULL)
  data <- preprocessed_data()
  mz_cols <- mz_columns_reg()
  
  if(input$metab_display == "TIC") {
    plot_data <- data
    plot_data$Intensity <- rowSums(data[, mz_cols, drop = FALSE])
    title <- "Preprocessed Metabolomics (Total Counts)"
  } else {
    req(input$selected_mz_reg)
    plot_data <- data
    plot_data$Intensity <- data[[input$selected_mz_reg]]
    title <- paste("Preprocessed Metabolomics (m/z:", input$selected_mz_reg, ")")
  }
  
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  ggplot(plot_data, aes(x = X, y = Y, color = Intensity)) +
    geom_point(size = 1) +
    scale_color_gradientn(colours = heatmap_Palette(100), name = "Intensity") +
    labs(title = title, x = "X", y = "Y") +
    coord_equal() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
})

# Registration Logic
#' @title Perform Image-Based Registration
#' @description Executes the core registration logic using the NiftyReg algorithm.
#' 1. Calculates resolution and handles differences.
#' 2. Normalizes images (Metabolomics TIC vs Transcriptomics TotalCounts).
#' 3. Runs Rigid or Affine registration based on user settings.
#' 4. Applies the resulting transformation to all metabolite channels.
#' 5. Interpolates values for the new grid.
#' @return A list containing:
#' - rgb_overlay: Normalized image data for visualization.
#' - aligned_metabolic_data: The transformed metabolomics data.
#' - mz_columns: List of metabolite names.
#' - quality: Registration quality metrics (MSE, Correlation).
registration_result <- eventReactive(input$run_registration, {
  if(preprocess_done()) {
    metabolic_plot_data <- preprocessed_data()
  } else {
    metabolic_plot_data <- metabolic_data_reg()
  }
  
  req(metabolic_plot_data, transcript_data_reg())
  
  withProgress(message = "Processing Registration...", value = 0.2, {
    transcript_data <- transcript_data_reg()
    transcript_plot_data <- data.frame(
      X = transcript_data$X,
      Y = transcript_data$Y,
      expression_value = transcript_data$TotalCounts
    )
    
    mz_cols <- mz_columns_reg()
    
    # Calculate Resolution
    calc_res <- function(d) {
      x_s <- sort(unique(d$X)); y_s <- sort(unique(d$Y))
      xr <- if(length(x_s)>1) min(diff(x_s)) else 1
      yr <- if(length(y_s)>1) min(diff(y_s)) else 1
      c(xr, yr)
    }
    
    metab_res <- calc_res(metabolic_plot_data)
    trans_res <- calc_res(transcript_plot_data)
    resolution_same <- isTRUE(all.equal(metab_res, trans_res, tolerance = 1e-6))
    
    # Prepare matrices for NiftyReg
    transcript_wide <- acast(transcript_plot_data, Y ~ X, value.var = "expression_value")
    transcript_final <- (transcript_wide - min(transcript_wide, na.rm=TRUE)) / (max(transcript_wide, na.rm=TRUE) - min(transcript_wide, na.rm=TRUE))
    transcript_final[is.na(transcript_final)] <- 0
    
    metabolic_plot_data$TIC <- rowSums(metabolic_plot_data[, mz_cols])
    metabolic_wide_ref <- acast(metabolic_plot_data, Y ~ X, value.var = "TIC")
    metabolic_norm_ref <- (metabolic_wide_ref - min(metabolic_wide_ref, na.rm=TRUE)) / (max(metabolic_wide_ref, na.rm=TRUE) - min(metabolic_wide_ref, na.rm=TRUE))
    metabolic_norm_ref[is.na(metabolic_wide_ref)] <- 0
    
    # Registration Strategy
    if(input$constrain_transform) {
      max_rotation_rad <- input$max_rotation * pi / 180
      
      # Rigid step
      rigid_result <- niftyreg(source=metabolic_norm_ref, target=transcript_final, scope="rigid", nLevels=3, maxIterations=300, symmetric=TRUE, interpolation=0)
      
      # Check rotation
      rigid_transform <- forward(rigid_result)
      rot_angle <- 0
      if(nrow(rigid_transform)==3 && ncol(rigid_transform)==3) {
        rot_angle <- atan2(rigid_transform[2,1], rigid_transform[1,1])
      }
      
      if(abs(rot_angle) > max_rotation_rad) {
        showNotification("Large rotation detected, using constrained affine.", type="warning")
        affine_result <- niftyreg(source=metabolic_norm_ref, target=transcript_final, scope="affine", nLevels=2, maxIterations=100, symmetric=TRUE, interpolation=0)
      } else {
        affine_result <- niftyreg(source=metabolic_norm_ref, target=transcript_final, scope="affine", nLevels=3, maxIterations=200, symmetric=TRUE, interpolation=0)
      }
    } else {
      affine_result <- niftyreg(source=metabolic_norm_ref, target=transcript_final, scope="affine", nLevels=3, maxIterations=500, symmetric=TRUE, interpolation=0)
    }
    
    # Apply Transform to Image
    aligned_m_resize <- applyTransform(forward(affine_result), metabolic_norm_ref)
    aligned_m_resize[is.na(aligned_m_resize)] <- 0
    
    # Overlay
    overlay <- array(0, dim=c(dim(transcript_final), 3))
    overlay[,,1] <- transcript_final
    overlay[,,3] <- aligned_m_resize * 0.5
    overlay[,,2] <- 0
    overlay_normalized <- (overlay - min(overlay)) / (max(overlay) - min(overlay))
    
    # Apply Transform to Data
    message("Applying transform to data...")
    aligned_mz_data <- list()
    interp_method <- if(resolution_same) 0 else 1
    
    for(mz_col in mz_cols) {
      metabolic_wide <- acast(metabolic_plot_data, Y ~ X, value.var = mz_col)
      aligned_mz <- applyTransform(forward(affine_result), metabolic_wide, interpolation=interp_method)
      aligned_mz_data[[mz_col]] <- aligned_mz
    }
    
    # Reconstruct DF
    n_rows <- nrow(aligned_mz_data[[1]])
    n_cols <- ncol(aligned_mz_data[[1]])
    trans_x_range <- range(transcript_plot_data$X)
    trans_y_range <- range(transcript_plot_data$Y)
    
    aligned_X <- seq(trans_x_range[1], trans_x_range[2], length.out = n_cols)
    aligned_Y <- seq(trans_y_range[1], trans_y_range[2], length.out = n_rows)
    
    aligned_grid <- expand.grid(X = aligned_X, Y = aligned_Y)
    
    for(mz_col in mz_cols) {
      aligned_grid[[mz_col]] <- as.vector(t(aligned_mz_data[[mz_col]]))
    }
    
    aligned_metabolic_data <- aligned_grid
    
    # Filter empty rows
    aligned_metabolic_data <- aligned_metabolic_data[rowSums(aligned_metabolic_data[, 3:ncol(aligned_metabolic_data)]!=0, na.rm=TRUE)>0, ]
    
    # Handle negative values
    aligned_metabolic_data[, 3:ncol(aligned_metabolic_data)] <- apply(
      aligned_metabolic_data[, 3:ncol(aligned_metabolic_data)], 2,
      function(x) { x[x<0] <- runif(sum(x<0)) * 1e-5; return(x) }
    )
    
    calc_mse <- mean((aligned_m_resize - transcript_final)^2, na.rm=TRUE)
    calc_cor <- cor(as.vector(aligned_m_resize), as.vector(transcript_final), use="complete.obs")
    
    list(
      rgb_overlay = overlay_normalized,
      aligned_metabolic_data = aligned_metabolic_data,
      mz_columns = mz_cols,
      quality = list(mse = calc_mse, correlation = calc_cor)
    )
  })
})

# Display Overlay
output$rgb_overlay <- renderPlot({
  result <- registration_result()
  grid.raster(result$rgb_overlay)
})

# Final Result Processing (Intersection)
#' @title Finalize Registered Data
#' @description Post-processes the registered metabolomics data.
#' Intersects the aligned metabolomics coordinates with the transcriptomics spots
#' to ensure a 1:1 spatial match. Creates the final Seurat object.
#' @return A Seurat object of registered metabolomics data.
metab_reg <- eventReactive(input$run_registration, {
  transcript_data_reg <- transcript_data_reg()
  result <- registration_result()
  final_data <- result$aligned_metabolic_data
  source("./source/preprocessing/runxy.R")
  # Intersection logic based on coordinates
  # Rounding might be needed if floats are slightly off
  final_data$X_Y <- paste0("sample:",final_data$X, "_", final_data$Y)
  transcript_data_reg$X_Y <- paste0("sample:",transcript_data_reg$X, "_", transcript_data_reg$Y)
  intersect_xy <- intersect(transcript_data_reg$X_Y, final_data$X_Y)
  final_data <- final_data %>%
    dplyr::filter(X_Y %in% intersect_xy) %>%
    dplyr::arrange(match(X_Y, intersect_xy))
  final_data$X_Y <- NULL
  final_data[is.na(final_data)] <- 0
  mz_mat <- create_sparse_from_wide(final_data) 
  coords <- data.frame(x=final_data$X,y=final_data$Y)
  metab_reg <- create_seurat(mz_mat, coords, binsize=1)
  metab_reg$x_y=paste0("x",metab_reg$x,"_y",metab_reg$y)
  return(metab_reg)
})

# Data Info Output
output$data_info_reg <- renderText({
  trans <- transcript_data_reg()
  metab <- metabolic_data_reg()
  req(trans, metab)
  
  info <- paste0(
    "=== Data Information ===\n",
    "Transcriptomics Spots: ", nrow(trans), "\n",
    "Metabolomics Spots: ", nrow(metab), "\n"
  )
  
  if(preprocess_done()) {
    info <- paste0(info, "Preprocessing: Done\n")
  }
  
  if(input$run_registration > 0) {
    res <- registration_result()
    metab_reg<-metab_reg()
    info <- paste0(info, 
                   "Registration: Done\n",
                   "Correlation: ", round(res$quality$correlation, 3), "\n",
                   "MSE: ", round(res$quality$mse, 5), "\n",
                   "Aligned Spots (Raw): ", nrow(res$aligned_metabolic_data), "\n",
                   "Intersected Spots: ", nrow(metab_reg@meta.data), "\n"
    )
  }
  return(info)
})

# Plot: Aligned Metabolic
#' @title Plot Registered Metabolomics
#' @description Renders a spatial plot of the final registered metabolomics data.
#' Displays either Total Ion Current (TIC) or a specific m/z channel for validation.
output$aligned_metabolic <- renderPlot({
  result <- registration_result()
  data <- result$aligned_metabolic_data
  mz_cols <- result$mz_columns
  
  if(input$metab_display == "TIC") {
    plot_data <- data
    plot_data$Intensity <- rowSums(data[, mz_cols])
    title <- "Registered Metabolomics (Total Counts)"
  } else {
    req(input$selected_mz_reg)
    plot_data <- data
    plot_data$Intensity <- data[[input$selected_mz_reg]]
    title <- paste("Registered Metabolomics (m/z:", input$selected_mz_reg, ")")
  }
  
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  ggplot(plot_data, aes(x=X, y=Y, color=Intensity)) +
    geom_point(size=1) +
    scale_color_gradientn(colours=heatmap_Palette(100), name="Intensity") +
    labs(title=title, x="X", y="Y") +
    coord_equal() +
    theme_minimal() +
    theme(plot.title = element_text(hjust=0.5))
})

# Compare Plot Generation
#' @title Generate Comparison Plots for Registration
#' @description Creates a list of four plots to visually assess registration quality:
#' 1. Original Transcriptomics.
#' 2. Original Metabolomics.
#' 3. Registered Metabolomics (Intersected).
#' 4. Transcriptomics (Intersected).
#' @return A list of ggplot objects.
aligned_compareplot_reg <- eventReactive(input$run_registration, {
  # Re-generate plots for download/view
  # Logic similar to original smr but adapted
  
  # P1: Original Trans
  trans <- transcript_data_reg()
  heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  p1 <- ggplot(trans, aes(x=X, y=Y, color=TotalCounts)) +
    geom_point(size=1) + scale_color_gradientn(colours=heatmap_Palette(100)) +
    labs(title="Original Transcript", x="X", y="Y") + coord_equal() + theme_minimal()
  
  # P2: Original Metab
  metab <- metabolic_data_reg()
  mz_cols <- mz_columns_reg()
  val <- rowSums(metab[, mz_cols])
  p2 <- ggplot(metab, aes(x=X, y=Y, color=val)) +
    geom_point(size=1) + scale_color_gradientn(colours=heatmap_Palette(100)) +
    labs(title="Original Metabolomics", x="X", y="Y") + coord_equal() + theme_minimal()
  
  # P3: Registered Metab (Intersected)
  metab_reg <- metab_reg()
  reg_metab <- data.frame(X=metab_reg$x,Y=metab_reg$y,val_reg=metab_reg$nCount_Spatial, val_reg=metab_reg$nCount_Spatial )
  p3 <- ggplot(reg_metab, aes(x=X, y=Y, color=val_reg)) +
    geom_point(size=1) + scale_color_gradientn(colours=heatmap_Palette(100)) +
    labs(title="Registered Metabolomics", x="X", y="Y") + coord_equal() + theme_minimal()
  
  # P4: Registered Trans (Intersected - same as original but filtered)
  # Actually transcript doesn't move. We just show it for comparison of coverage.
  trans$X_Y <- paste0(round(trans$X, 2), "_", round(trans$Y, 2))
  reg_metab$X_Y <- paste0(round(reg_metab$X, 2), "_", round(reg_metab$Y, 2))
  intersect_xy <- intersect(trans$X_Y, reg_metab$X_Y)
  trans_filt <- trans %>% dplyr::filter(X_Y %in% intersect_xy)
  
  p4 <- ggplot(trans_filt, aes(x=X, y=Y, color=TotalCounts)) +
    geom_point(size=1) + scale_color_gradientn(colours=heatmap_Palette(100)) +
    labs(title="Transcript (Intersected)", x="X", y="Y") + coord_equal() + theme_minimal()
  
  list(p2, p1, p3, p4)
})

# Render Comparison Grid
#' @title Render Registration Comparison Grid
#' @description Displays the four comparison plots generated by `aligned_compareplot_reg`
#' in a 2x2 grid structure.
output$aligned_compare <- renderPlot({
  plots <- aligned_compareplot_reg()
  gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol=2)
})

# Downloads
#' @title Download RGB Overlay
#' @description Enables downloading the registration comparison grid as a PNG image.
output$download_rgb_png <- downloadHandler(
  filename = function() { "registration_comparison.png" },
  content = function(file) {
    plots <- aligned_compareplot_reg()
    png(file, width=800, height=800)
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol=2)
    dev.off()
  }
)

#' @title Download All Registration Plots
#' @description Generates a ZIP archive containing all individual plots related to registration:
#' - Original Transcriptomics/Metabolomics
#' - Preprocessed Metabolomics
#' - Registered Metabolomics
#' - RGB Overlay
#' - Comparison Matrix
output$download_reg_plots_zip <- downloadHandler(
  filename = function() {
    paste0("registration_plots_", Sys.Date(), ".zip")
  },
  content = function(file) {
    req(transcript_data_reg(), metabolic_data_reg(), registration_result())
    
    # Create a temporary directory
    temp_dir <- tempdir()
    plot_files <- c()
    
    withProgress(message = 'Generating plots...', value = 0, {
      
      # 1. Original Transcriptomics
      incProgress(0.1, detail = "Original Transcriptomics")
      trans_data <- transcript_data_reg()
      heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
      p1 <- ggplot(trans_data, aes(x = X, y = Y, color = TotalCounts)) +
        geom_point(size = 1) +
        scale_color_gradientn(colours = heatmap_Palette(100), name = "Expression") +
        labs(title = "Original Transcriptomics", x = "X", y = "Y") +
        coord_equal() +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
      
      f1 <- file.path(temp_dir, "Original_Transcriptomics.png")
      ggsave(f1, plot = p1, width = 6, height = 5, dpi = 300)
      plot_files <- c(plot_files, f1)
      
      # 2. Original Metabolomics
      incProgress(0.2, detail = "Original Metabolomics")
      metab_data <- metabolic_data_reg()
      mz_cols <- mz_columns_reg()
      
      if(input$metab_display == "TIC") {
        plot_data <- metab_data
        plot_data$Intensity <- rowSums(metab_data[, mz_cols, drop = FALSE])
        title <- "Original Metabolomics"
      } else {
        if(!is.null(input$selected_mz_reg)) {
          plot_data <- metab_data
          plot_data$Intensity <- metab_data[[input$selected_mz_reg]]
          title <- paste("Original Metabolomics (m/z:", input$selected_mz_reg, ")")
        } else {
          # Fallback if specific not selected
          plot_data <- metab_data
          plot_data$Intensity <- rowSums(metab_data[, mz_cols, drop = FALSE])
          title <- "Original Metabolomics (TIC Fallback)"
        }
      }
      
      p2 <- ggplot(plot_data, aes(x = X, y = Y, color = Intensity)) +
        geom_point(size = 1) +
        scale_color_gradientn(colours = heatmap_Palette(100), name = "Intensity") +
        labs(title = title, x = "X", y = "Y") +
        coord_equal() +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
      
      f2 <- file.path(temp_dir, "Original_Metabolomics.png")
      ggsave(f2, plot = p2, width = 6, height = 5, dpi = 300)
      plot_files <- c(plot_files, f2)
      
      # 3. Preprocessed Metabolomics
      incProgress(0.2, detail = "Preprocessed Metabolomics")
      if(preprocess_done()) {
        pre_data <- preprocessed_data()
        
        if(input$metab_display == "TIC") {
          plot_data_pre <- pre_data
          plot_data_pre$Intensity <- rowSums(pre_data[, mz_cols, drop = FALSE])
          title_pre <- "Preprocessed Metabolomics (Total Counts)"
        } else {
          if(!is.null(input$selected_mz_reg)) {
            plot_data_pre <- pre_data
            plot_data_pre$Intensity <- pre_data[[input$selected_mz_reg]]
            title_pre <- paste("Preprocessed Metabolomics (m/z:", input$selected_mz_reg, ")")
          } else {
            plot_data_pre <- pre_data
            plot_data_pre$Intensity <- rowSums(pre_data[, mz_cols, drop = FALSE])
            title_pre <- "Preprocessed Metabolomics (TIC Fallback)"
          }
        }
        
        p3 <- ggplot(plot_data_pre, aes(x = X, y = Y, color = Intensity)) +
          geom_point(size = 1) +
          scale_color_gradientn(colours = heatmap_Palette(100), name = "Intensity") +
          labs(title = title_pre, x = "X", y = "Y") +
          coord_equal() +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
        
        f3 <- file.path(temp_dir, "Preprocessed_Metabolomics.png")
        ggsave(f3, plot = p3, width = 6, height = 5, dpi = 300)
        plot_files <- c(plot_files, f3)
      }
      
      # 4. Registered Metabolomics
      incProgress(0.2, detail = "Registered Metabolomics")
      res <- registration_result()
      reg_data <- res$aligned_metabolic_data
      mz_cols_res <- res$mz_columns
      
      if(input$metab_display == "TIC") {
        plot_data_reg <- reg_data
        plot_data_reg$Intensity <- rowSums(reg_data[, mz_cols_res])
        title_reg <- "Registered Metabolomics (Total Counts)"
      } else {
        if(!is.null(input$selected_mz_reg)) {
          plot_data_reg <- reg_data
          plot_data_reg$Intensity <- reg_data[[input$selected_mz_reg]]
          title_reg <- paste("Registered Metabolomics (m/z:", input$selected_mz_reg, ")")
        } else {
          plot_data_reg <- reg_data
          plot_data_reg$Intensity <- rowSums(reg_data[, mz_cols_res])
          title_reg <- "Registered Metabolomics (TIC Fallback)"
        }
      }
      
      p4 <- ggplot(plot_data_reg, aes(x=X, y=Y, color=Intensity)) +
        geom_point(size=1) +
        scale_color_gradientn(colours=heatmap_Palette(100), name="Intensity") +
        labs(title=title_reg, x="X", y="Y") +
        coord_equal() +
        theme_minimal() +
        theme(plot.title = element_text(hjust=0.5))
      
      f4 <- file.path(temp_dir, "Registered_Metabolomics.png")
      ggsave(f4, plot = p4, width = 6, height = 5, dpi = 300)
      plot_files <- c(plot_files, f4)
      
      # 5. RGB Overlay
      incProgress(0.1, detail = "RGB Overlay")
      f5 <- file.path(temp_dir, "RGB_Overlay.png")
      png(f5, width = 2000, height = 2000, res=300)
      grid.raster(res$rgb_overlay)
      dev.off()
      plot_files <- c(plot_files, f5)
      
      # 6. Comparison Matrix
      incProgress(0.2, detail = "Comparison Matrix")
      plots_comp <- aligned_compareplot_reg()
      f6 <- file.path(temp_dir, "Comparison_Matrix.png")
      png(f6, width = 3000, height = 3000, res=300)
      gridExtra::grid.arrange(plots_comp[[1]], plots_comp[[2]], plots_comp[[3]], plots_comp[[4]], ncol=2)
      dev.off()
      plot_files <- c(plot_files, f6)
      
      # Zip files
      current_wd <- getwd()
      setwd(temp_dir)
      utils::zip(file, files = basename(plot_files))
      setwd(current_wd)
    })
  }
)




#############################
#' @title Create Shared Sample List
#' @description Identifies the intersection of spots between the final metabolomics and transcriptomics datasets.
#' Creates a mapping dataframe for subsequent processing.
#' @return A data frame with "gene" and "metabolomics" columns containing shared spot IDs.
samplelist0<-reactive({
  # Safeguard: Ensure we have explicit user action
  if (!input$do_registration && input$Submit == 0) return(NULL)
  if (input$do_registration && input$run_registration == 0) return(NULL)
  
  req(pre_genelist(), final_metabolomics_list())
  
  withProgress(message = "Processing data...",value=0.8,{
    
    metabolomicslist<-final_metabolomics_list()
    genelist<-pre_genelist()
    
    intersection <- intersect(metabolomicslist$x_y, genelist$x_y)
    
    samplelist<-data.frame(gene=intersection,metabolomics=intersection) |>
      dplyr::mutate(group="0")
    
    return(samplelist)
  })
})

#' @title Initialize Seurat Objects
#' @description Constructs the final Seurat objects for both modalities using the shared sample list.
#' Calculates basic QC metrics (nFeature_Spatial, nCount_Spatial) and ensures metadata consistency.
#' @return A list containing:
#' - data_mrds: Metabolomics Seurat object.
#' - data_trds: Transcriptomics Seurat object.
data_rds<-reactive({
  req(samplelist0())
  
  withProgress(message = "Processing data...",value=0.8,{
    genelist<-pre_genelist()
    metabolomicslist<-final_metabolomics_list()
    
    samplelist<-samplelist0()
    source("./source/preprocessing/run_prerds.R")
    
    samplelist<-samplelist %>%
      dplyr::select(-metabolomics)
    colnames(samplelist)=c("x_y","group")
    data_mrds<-runrds(rds=metabolomicslist,samplelist)
    data_trds<-runrds(rds=genelist,samplelist)
    
    # Revert to standard processing for initial view (Raw/Basic Distribution)
    # The sophisticated preprocessing (TIC/RMS) is now configured in the Clustering step
    data_mrds <- RUNSCT(data_mrds)
    data_trds <- RUNSCT(data_trds)
    
    overall_test<-function(rds){
      
      spatial_matrix <- GetAssayData(rds[["Spatial"]], layer ="counts")
      nCount_Spatial <- colSums(spatial_matrix)
      rds[["nCount_Spatial"]] <- nCount_Spatial
      nFeature_Spatial <- colSums(spatial_matrix > 0)
      rds[["nFeature_Spatial"]] <- nFeature_Spatial
      rds@meta.data$cell <- paste0("sample:",rds$x,"_",rds$y)
      return(rds)
    }
    data_mrds<-overall_test(data_mrds)
    data_trds<-overall_test(data_trds)
    
    
    data<-list(data_mrds=data_mrds,data_trds=data_trds)
    
    return(data)
  })
})

#' @title Display Basic Dataset Info
#' @description Calculates and renders a summary table of the processed datasets, including:
#' - Number of spots (points).
#' - Dimensions (rows/cols).
#' - Feature counts (genes/metabolites).
output$basic_info <- renderTable({
  data_rds<-data_rds()
  run_basic_info<-function(data){
    ions <- rownames(data@assays$Spatial$counts)
    point_number <- length(unique(colnames(data@assays$Spatial$counts)))
    peak_number <- length(ions)
    number_of_rows <- length(unique(data@meta.data$x))
    number_of_cols <- length(unique(data@meta.data$y))
    summary<-c(point_number,number_of_rows,number_of_cols,peak_number)
    return(summary)
  }
  data_msummary<-run_basic_info(data_rds[[1]])
  data_tsummary<-run_basic_info(data_rds[[2]])
  
  data.frame(omics=c("Metabolomics","transcriptomics"),
             point_number=c(data_msummary[1],data_tsummary[1]),
             number_of_rows=c(data_msummary[2],data_tsummary[2]),
             number_of_cols=c(data_msummary[3],data_tsummary[3]),
             feature_number=c(data_msummary[4],data_tsummary[4]))
})

#' @title Generate Overview Quality Plots
#' @description Creates spatial visualizations of total counts (nCount) for both modalities.
#' Used to verify the spatial distribution of signal intensity after processing.
#' @return A list of two ggplot objects (Metabolite and Gene).
overall_plotlist <-reactive({
  req(data_rds())
  data_rds<-data_rds()
  source("./source/OverallAnalysisFunction/Clustering/clusterplot.R")
  data_mrds=data_rds[[1]]
  data_trds=data_rds[[2]]
  
  overallplot <- function(obj,type,pointSize,breakseq) {
    
    
    plot1 <- Preview(obj, 'nCount_Spatial',pointSize,breakseq)+
      labs(title=paste0(type))+
      theme(plot.title = element_text(hjust = 0.5))
    
    return(plot1)
  }
  mplot<-overallplot(data_mrds,type="Metabolite",pointSize=1,breakseq=50)
  tplot<-overallplot(data_trds,type="Gene",pointSize=1,breakseq=50)
  plotlist<-list(mplot,tplot)
  return(plotlist)
})


#' @title Render Total Counts Plot
#' @description Displays the overview quality plots (Metabolite and Gene total counts) side-by-side.
output$totalcounts_plot <- renderPlot({
  req(overall_plotlist())
  withProgress(message = "Plotting...",value=0.8,{
    overall_plotlist<-overall_plotlist()
    mplot<-overall_plotlist[[1]]
    tplot<-overall_plotlist[[2]]
    p<-gridExtra::grid.arrange(mplot,tplot,ncol=2)
    p
    assign("totalcounts_plot", p, envir = .GlobalEnv)
  })
})



#' @title Dynamic File Upload UI
#' @description Renders the appropriate file input buttons based on the selected data type (TXT or RDS).
#' Returns NULL if using demo data.
output$file_button_container <- renderUI({
  if (input$demo_select=="Upload txt data") {
    
    tagList(
      fileInput("metabfile", "Upload spatial metabolomics txt file",
                accept = ".txt"),
      fileInput("transfile", "Upload spatial transcriptomics txt file",
                accept = ".txt")
      
    )
  } else if (input$demo_select=="Upload rds data") {
    tagList(
      fileInput("metabfile", "Upload spatial metabolomics rds file",
                accept = ".rds"),
      fileInput("transfile", "Upload spatial transcriptomics rds file",
                accept = ".rds")
      
      
    )
  }else{
    NULL
  }
})

#' @title Update Demo Defaults
#' @description Automatically sets default parameters (species: mmu, mode: pos) when "Use demo data" is selected.
observe({
  if (input$demo_select == "Use demo data") {
    updateSelectizeInput(session=getDefaultReactiveDomain(), "speciesname_select", selected = "mmu")
    updateSelectizeInput(session=getDefaultReactiveDomain(), "metab_mode", selected = "pos")
  }
})

