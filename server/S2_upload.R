diff_omics <- reactiveVal(NULL)
spatial_pattern_ionlist <- reactiveVal(NULL)



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
      validate(
        need(is.numeric(data$x), "x column must be numeric"),
        need(is.numeric(data$y), "y column must be numeric")
      )
      data<-colname_change(data)
    }
  }else if(input$demo_select == "Upload rds data"){
    if(!is.null(input$transfile) && input$transfile$name != ""){
      data <- readRDS(input$transfile$datapath)
      validate(
        need( data@active.assay =="Spatial" && "Spatial" %in% names(data@assays), 
              "Active assay must be named 'Spatial'!"),
        
        need( c("x","y") %in% colnames(data@meta.data),
              "The meta.data of Seurat object must contain two columns named 'x' and 'y'."),
        
        need(is.numeric(data@meta.data$x), "x column must be numeric"),
        need(is.numeric(data@meta.data$y), "y column must be numeric")
      )
    }
  }
    return(data)
  }) 
}) 

metabdata<-eventReactive(c(input$Submit), {
  source("./source/preprocessing/runxy.R")
  withProgress(message = "Loading data ...",value=0.8,{

  if(input$demo_select == "Use demo data"){
    data<-demometab_rds()
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
      validate(
        need(is.numeric(data$x), "x column must be numeric"),
        need(is.numeric(data$y), "y column must be numeric")
      )
      data<-colname_change(data)
    }
    }else if(input$demo_select == "Upload rds data"){
      if(!is.null(input$metabfile) && input$metabfile$name != ""){
        data <- readRDS(input$metabfile$datapath)
        validate(
        need( data@active.assay =="Spatial" && "Spatial" %in% names(data@assays), 
             "Active assay must be named 'Spatial'!"),
        
        need( c("x","y") %in% colnames(data@meta.data),
             "The meta.data of Seurat object must contain two columns named 'x' and 'y'."),
         
          need(is.numeric(data@meta.data$x), "x column must be numeric"),
          need(is.numeric(data@meta.data$y), "y column must be numeric")
        )
      }
  }
  
  return(data)
  })
})


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




pre_metabolomicslist<-eventReactive(c(input$Submit), {
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
pre_genelist<-eventReactive(c(input$Submit), {
  req(pre_metabolomicslist())
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


samplelist0<-eventReactive(input$Submit, {
  req(pre_genelist(),pre_metabolomicslist())
  withProgress(message = "Processing data...",value=0.8,{
      
    metabolomicslist<-pre_metabolomicslist()
    genelist<-pre_genelist()

  intersection <- intersect(metabolomicslist$x_y, genelist$x_y)

  samplelist<-data.frame(gene=intersection,metabolomics=intersection) |>
    dplyr::mutate(group="0")

  return(samplelist)
  })
})
data_rds<-eventReactive(input$Submit, {
  req(samplelist0())
  withProgress(message = "Processing data...",value=0.8,{
    genelist<-pre_genelist()
    metabolomicslist<-pre_metabolomicslist()

    samplelist<-samplelist0()
    source("./source/preprocessing/run_prerds.R")

    
    samplelist<-samplelist %>%
      dplyr::select(-metabolomics)
    colnames(samplelist)=c("x_y","group")
    data_mrds<-runrds(rds=metabolomicslist,samplelist)
    data_trds<-runrds(rds=genelist,samplelist)

    overall_test<-function(rds){

      spatial_matrix <- GetAssayData(rds[["Spatial"]], layer ="counts")
      nCount_Spatial <- colSums(spatial_matrix)
      rds[["nCount_Spatial"]] <- nCount_Spatial
      nFeature_Spatial <- colSums(spatial_matrix > 0)
      rds[["nFeature_Spatial"]] <- nFeature_Spatial
      return(rds)
    }
    data_mrds<-overall_test(data_mrds)
    data_trds<-overall_test(data_trds)

    data<-list(data_mrds=data_mrds,data_trds=data_trds)

    return(data)
    })
})

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

overall_plotlist <-eventReactive(input$Submit, {
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

observe({
          if (input$demo_select == "Use demo data") {
            updateSelectizeInput(session=getDefaultReactiveDomain(), "speciesname_select", selected = "mmu")
            updateSelectizeInput(session=getDefaultReactiveDomain(), "metab_mode", selected = "pos")
          }
        })
        
