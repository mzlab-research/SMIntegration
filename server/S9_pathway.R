# ==============================================================================
# S9_pathway.R
# Server logic for Step 5: Functional Association and Annotation
# Handles pathway mapping, annotation statistics, and visualization (DotPlot, Venn).
# ==============================================================================

# ------------------------------------------------------------------------------
# Download Handlers
# ------------------------------------------------------------------------------

#' @title Download Demo Transcriptomics Annotation
output$download_Gene_pathway_annotation_demo <- downloadHandler(
  filename = function() {
    paste0("gene_identi_demo.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData=read_delim("./example_data/trans_identi.xls") |>
        dplyr::select(gene,KEGG.ID)
      write_delim(DemoData,file, delim  = '\t')
    })
  }
)

#' @title Download Demo Metabolomics Annotation
output$download_metab_pathway_annotation_demo <- downloadHandler(
  filename = function() {
    paste0("metabolite_identi_demo.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
      DemoData=read_delim("./example_data/metab_identi.xls") |>
        dplyr::select(mz,metabolite,KEGG.ID)
     # colnames(DemoData)[2]<-"Name"
      write_delim(DemoData,file, delim = '\t')
    })
  }
)

#' @title Download Full Annotation Results
#' @description Exports the annotation count table (DotPlot data).
output$download_annotation_data_all <- downloadHandler(
  filename = function() {
    paste0("annotation_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    DemoData=dotplot_data()[[2]]
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })

#' @title Download Pathway Image (KEGG Map)
output$download_pathway_annotation_plot <- downloadHandler(
  filename = function() {
    DemoData=dotplot_data()[[2]] %>%
      dplyr::filter(Pathway==input$Pathway_select)
    paste0(DemoData$PathwayID, ".png")
  },
  content = function(file) {
    req(temp_pathwayannotation_file())
    withProgress(message = 'Downloading file...', value = 0.7, {
     file.copy(temp_pathwayannotation_file(), file, overwrite = TRUE)
  })
})

#' @title Download Pathway Node Data
#' @description Exports features mapped to the selected pathway (Metab + Gene).
output$download_pathway_annotation_data <- downloadHandler(
  filename = function() {

    paste0(input$Pathway_select,"_data.zip")
  },
  content = function(file) {
    withProgress(message = 'Downloading files...', value = 0.7, {
    annotation_plotdata<-annotation_plotdata()
    annotation_data<-list(annotation_plotdata[[3]],annotation_plotdata[[4]])
    tempdir <- setwd(tempdir())
    on.exit(setwd(tempdir))
    fi=c("Metabolite_data.csv","Gene_data.csv")
    for (i in 1:length(fi)) {
      data<-annotation_data[[i]]
      write.csv(data, fi[i], row.names = FALSE)
    }
    zip(file,fi)
    })
  })



#' @title Download Spatial Feature Plot (Metabolite)
output$download_pathway_annotation_m_plot <- downloadHandler(
  filename = function() {
    "Metabolite_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- pathway_annotation_m_plot_save()

    ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)
  })
})
output$download_pathway_annotation_m_data <- downloadHandler(
  filename = function() {
    paste0("Metabolite_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    plotdatalist <- pathway_annotation_ion_data()
    DemoData<-plotdatalist[[1]]
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })
output$download_pathway_annotation_t_plot <- downloadHandler(
  filename = function() {
    "Gene_plot.png"
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- pathway_annotation_t_plot_save()

    ggsave(file, plot = p, device = "png", width = 6, height = 6,bg = "#FFFFFF", dpi = 300)
  })
})
output$download_pathway_annotation_t_data <- downloadHandler(
  filename = function() {
    paste0("Gene_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    plotdatalist <- pathway_annotation_ion_data()
    DemoData<-plotdatalist[[2]]
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })

#' @title Download Annotation Dot Plot Data
output$download_dotplot_data <- downloadHandler(
  filename = function() {

    paste0("annotation_dotplot_data.txt")
  },
  content = function(file) {
    withProgress(message = 'Downloading file...', value = 0.7, {
    DemoData=dotplot_data_save()
    write.table(DemoData,file,row.names = F, quote = F, sep = '\t')
    })
  })

#' @title Download Annotation Dot Plot Image
output$download_dotplot <- downloadHandler(
  filename = function() {
    "annotation_dotplot.png"
  },
  content = function(file) {

    withProgress(message = 'Downloading file...', value = 0.7, {
    p <- dotplot_save()

    ggsave(file, plot = p, device = "png", width = 8, height = 6,bg = "#FFFFFF", dpi = 300)
  })
})

#' @title Download Venn Plot
output$download_venn_plot <- downloadHandler(
  filename = function() {
    "venn_plot.png"
  },
  content = function(file) {
   req(temp_venn_plot_file())
    withProgress(message = 'Downloading file...', value = 0.7, {
    file.copy(temp_venn_plot_file(), file, overwrite = TRUE)
    })
  })




# ------------------------------------------------------------------------------
# ID Mapping and Database Setup
# ------------------------------------------------------------------------------

#' @title Metabolite Identification Mapping
#' @description Maps input metabolite names or m/z values to KEGG IDs.
#' Uses built-in database (level 3 identification) or user-uploaded annotation file.
metabidenti<-eventReactive(c(input$start_annotation), {
  if(input$demo_select == "Use demo data"){
    metabidenti=read_delim("./example_data/metab_identi.xls")
    if("mz" %in% colnames(metabidenti)){
      metabidenti$mz<-as.character(metabidenti$mz)
    }
    print("metabidenti done")
    return(metabidenti)
  }else{
    if (input$annotation_select == "use built-in database") {
      req(data_rds())
      withProgress(message = "Processing data...",value=0.8,{
        print("metabidenti start")
        source("./source/OverallAnalysisFunction/Annotation/metab_identi_new.R")
        data<-data_rds()[[1]]
        mode<-input$metab_mode
        sp<-run_sp(data)
        input_type <- detect_input_type(sp)
        if(input_type == "metab_name") {
          metabidenti <- run_metab_name_processing(sp,mode)
        }else{
          # MS-based identification logic
          combine<-run_isotope(sp,mode)
          end_deal_data<-run_add(combine,mode)
          more_to_one_identi<-run_metab_identi(end_deal_data,mode)
          metabidenti<-run_more_to_one(more_to_one_identi)
          metabidenti$metabolite<-as.character(metabidenti$mz)
        }
        
        if("KEGG ID" %in% colnames(metabidenti)){
          metabidenti<-metabidenti %>%
            rename(KEGG.ID=`KEGG ID`)
        }
        print("metabidenti done")
        return(metabidenti)
      })
    }else{
      # User uploaded annotation
      if(!is.null(input$metabannotationfile) && input$metabannotationfile$name != ""){
        metabidenti <- read_delim(input$metabannotationfile$datapath)
        if("KEGG ID" %in% colnames(metabidenti)){
          metabidenti<-metabidenti %>%
            rename(KEGG.ID=`KEGG ID`)
        }
        validate(
          need(ncol(metabidenti) >= 2, 
               "Data must have at least 2 columns (metabolite and KEGG.ID)"),
          
          need("metabolite" %in% colnames(metabidenti) && "KEGG.ID" %in% colnames(metabidenti),
               "First two columns must be named 'metabolite' and 'KEGG.ID'")       

        )
        metabidenti$metabolite<-as.character(metabidenti$metabolite)
        return(metabidenti)
      }else{
        # Fallback to built-in logic if no file uploaded
        req(data_rds())
        withProgress(message = "Processing data...",value=0.8,{
          print("metabidenti start")
          source("./source/OverallAnalysisFunction/Annotation/metab_identi_new.R")
          data<-data_rds()[[1]]
          mode<-input$metab_mode
          sp<-run_sp(data)
          input_type <- detect_input_type(sp)
          if(input_type == "metab_name") {
            metabidenti <- run_metab_name_processing(sp,mode)
          }else{
            combine<-run_isotope(sp,mode)
            end_deal_data<-run_add(combine,mode)
            more_to_one_identi<-run_metab_identi(end_deal_data,mode)
            metabidenti<-run_more_to_one(more_to_one_identi)
            metabidenti$metabolite<-as.character(metabidenti$mz)
          }
          if("KEGG ID" %in% colnames(metabidenti)){
            metabidenti<-metabidenti %>%
              rename(KEGG.ID=`KEGG ID`)
          }
          metabidenti$metabolite<-as.character(metabidenti$metabolite)
          
          print("metabidenti done")
          return(metabidenti)
        })
      }
    }
  }
})

#' @title Transcript ID Mapping
#' @description Maps gene symbols to Entrez IDs and then to KEGG KO IDs.
#' Supports Human (hsa) and Mouse (mmu).
transidenti<-eventReactive(c(input$start_annotation), {
  if(input$demo_select == "Use demo data"){
    transidenti=read_delim("./example_data/trans_identi.xls")
    return(transidenti)
  }else{
  if (input$annotation_select == "use built-in database") {
    req(data_rds())
    withProgress(message = "Processing data...",value=0.8,{

      data<-data_rds()[[2]]
      genename<-rownames(data@assays$Spatial$counts)
      speciesname=input$speciesname_select
      
      # Convert Symbol to Entrez ID
      if(speciesname=="mmu"){
        id <- bitr(genename,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Mm.eg.db)
      }else if(speciesname=="hsa"){
        id <- bitr(genename,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)    
      }
      source("./source/OverallAnalysisFunction/Annotation/gene_identi.R")
      transidenti <-trans_identi_processing(id,speciesname)
      transidenti$gene<-as.character(transidenti$gene)

      return(transidenti)
    })
  }else{
    # User uploaded annotation
    if(!is.null(input$transannotationfile) && input$transannotationfile$name != ""){
      transidenti <- read_delim(input$transannotationfile$datapath,delim = "\t")
      validate(
        need(ncol(transidenti) >= 2, 
             "Data must have at least 2 columns (gene and KEGG.ID)"),
        
          need("gene" %in% colnames(transidenti) && "KEGG.ID" %in% colnames(transidenti),
               "First two columns must be named 'gene' and 'KEGG.ID'")

      )
      transidenti$gene<-as.character(transidenti$gene)
      return(transidenti)
    }else{
      # Fallback
      req(data_rds())
      withProgress(message = "Processing data...",value=0.8,{
        data<-data_rds()[[2]]
        genename<-rownames(data@assays$Spatial$counts)
        speciesname=input$speciesname_select

        if(speciesname=="mmu"){
          id <- bitr(genename,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Mm.eg.db)
        }else if(speciesname=="hsa"){
          id <- bitr(genename,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)    
        }
        source("./source/OverallAnalysisFunction/Annotation/gene_identi.R")
        transidenti <-trans_identi_processing(id,speciesname)
        transidenti$gene<-as.character(transidenti$gene)

        return(transidenti)
      })
    }
  }
  }
})

#' @title Initialize Annotation Database
#' @description Loads KEGG pathway definitions for the selected species.
annotation_database<-eventReactive(c(input$start_annotation), {
  species=input$speciesname_select
  print(species)
  req(species)
    withProgress(message = "Processing data...",value=0.8,{
      source("./source/Pathway/annotation_stats.R")
      metabdatabse<-annotationdatabase(omics="metab",species=species)
      transdatabse<-annotationdatabase(omics="trans",species=species)
      annotation_database<-list(metabdatabse,transdatabse)
      return(annotation_database)
    })

  
})

metabannotationdata<-eventReactive(c(input$start_annotation), {

  if(!is.null(input$functional_data) && input$functional_data!=""){
    if(input$functional_data=="Differential features"){
      req(diff_omics(),metabidenti())
      withProgress(message = "Processing data...",value=0.8,{
        print("metabannotationdata start")
        diff_omics<-diff_omics()

        metabidenti<-metabidenti()
        diff_m<-diff_omics[[1]] %>%
          dplyr::filter(State!="Non-significant") %>%
          dplyr::select(metabolite,p_val_adj,`log2(Fold Change)`,State)

        #up down all
        if(!is.null(input$functional_diffdata) && input$functional_diffdata!=""){
        if(input$functional_diffdata=="Up-regulated features"){
          diff_m<-diff_m %>%
            dplyr::filter(State=="Up")
          validate(
            need(nrow(diff_m)>0,"The number of up-regulated metabolites was 0")
          )
        }else if(input$functional_diffdata=="Down-regulated features"){
          diff_m<-diff_m %>%
            dplyr::filter(State=="Down")
          validate(
            need(nrow(diff_m)>0,"The number of down-regulated metabolites was 0")
          )
        }else{
          validate(
            need(nrow(diff_m)>0,"The number of differential metabolites was 0")
          )
        }}
        if(nrow(diff_m) > 300){
          diff_m %<>% arrange(p_val_adj)
          diff_m <- diff_m[1:300,]
        }
          diff_m$metabolite<-as.character(diff_m$metabolite)
          mzid <-unique(diff_m$metabolite)
          data<-data.frame(metabolite=mzid) %>%
            left_join(metabidenti,by=c("metabolite"="metabolite"))
          data$metabolite<-as.character(data$metabolite)
          data <-data[!is.na(data$KEGG.ID), ]
          metabannotationdata<-data %>%
            left_join(diff_m,by=c("metabolite"="metabolite"))
          print("metabannotationdata done")
          return(metabannotationdata)
        
        
      })
    }else{
      req(spatial_pattern_ionlist(),metabidenti())
      withProgress(message = "Processing data...",value=0.8,{
        print("metabannotationdata start")
        metabidenti<-metabidenti()
        spatial_pattern_ionlist<-spatial_pattern_ionlist()
        mz<-spatial_pattern_ionlist[[1]]
        validate(
          need(length(mz)>0,"The number of pattern-specific metabolites was 0")
        )
          mzid <-unique(as.character(mz))
          data<-data.frame(metabolite=mzid) %>%
            left_join(metabidenti,by=c("metabolite"="metabolite"))
          data$metabolite<-as.character(data$metabolite)
          metabannotationdata <-data[!is.na(data$KEGG.ID), ]
          print("metabannotationdata done")
          return(metabannotationdata)

        
      })
    }
  }


})
#' @title Prepare Transcriptomics Annotation Data
#' @description Prepares the gene list for pathway analysis based on user selection:
#' 1. Differential features: Filters by Up/Down regulation and P-value (top 300).
#' 2. Pattern-specific features: Uses genes identified in spatial pattern analysis.
#' Maps genes to KEGG IDs using the loaded identifier mapping.
#' @return A data frame containing genes, KEGG IDs, and statistical metrics.
transannotationdata<-eventReactive(c(input$start_annotation), {
  
  if(!is.null(input$functional_data) && input$functional_data!=""){
    if(input$functional_data=="Differential features"){
      req(diff_omics(),transidenti())
      withProgress(message = "Processing data...",value=0.8,{
        print("transannotationdata start")
        diff_omics<-diff_omics()
        transidenti<-transidenti()
        diff_t<-diff_omics[[2]] %>%
          dplyr::filter(State!="Non-significant") %>%
          dplyr::select(gene,p_val_adj,`log2(Fold Change)`,State)
        
        #up down all
        if(!is.null(input$functional_diffdata) && input$functional_diffdata!=""){
          if(input$functional_diffdata=="Up-regulated features"){
            diff_t<-diff_t %>%
              dplyr::filter(State=="Up")
            validate(
              need(nrow(diff_t)>0,"The number of up-regulated genes was 0")
            )
          }else if(input$functional_diffdata=="Down-regulated features"){
            diff_t<-diff_t %>%
              dplyr::filter(State=="Down")
            validate(
              need(nrow(diff_t)>0,"The number of down-regulated genes was 0")
            )
          }else{
            validate(
              need(nrow(diff_t)>0,"The number of differential genes was 0")
            )
          }}
        if(nrow(diff_t) > 300){
          diff_t %<>% arrange(p_val_adj)
          diff_t <- diff_t[1:300,]
        }
        genename=as.character(unique(diff_t$gene))
        data<-data.frame(gene=genename) %>%
          left_join(transidenti,by=c("gene"="gene"))
        data <- data[!is.na(data$KEGG.ID), ]
        transannotationdata<-data %>%
          left_join(diff_t,by=c("gene"="gene"))
        print("transannotationdata done")
        return(transannotationdata)
        
      })
    }else{
      req(spatial_pattern_ionlist(),transidenti())
      withProgress(message = "Processing data...",value=0.8,{
        print("transannotationdata start")
        transidenti<-transidenti()
        spatial_pattern_ionlist<-spatial_pattern_ionlist()
        
        gene<-spatial_pattern_ionlist[[2]] 
        validate(
          need(length(gene)>0,"The number of pattern-specific genes was 0")
        )
        genename=as.character(unique(gene))
        data<-data.frame(gene=genename) %>%
          left_join(transidenti,by=c("gene"="gene"))
        transannotationdata <- data[!is.na(data$KEGG.ID), ]
        
        print("transannotationdata done")
        
        return(transannotationdata)
        
      })
    }
  }
  
})

#' @title Run Pathway annotation Analysis
#' @description Performs KEGG pathway annotation analysis for both Metabolomics and Transcriptomics.
#' Uses the `run_annotation_count` function (likely hypergeometric test or count-based) against the loaded species database.
#' @return A list containing annotation results for Metabolomics (1) and Transcriptomics (2).
annotation_results<-eventReactive(c(input$start_annotation), {
  req(transannotationdata(),metabannotationdata(),annotation_database())
  if(nrow(transannotationdata())>0 && nrow(metabannotationdata())>0){
    withProgress(message = "Processing data...",value=0.8,{
      print("annotation_results start")
      source("./source/Pathway/annotation_stats.R")
      metabannotationdata<-metabannotationdata()
      transannotationdata<-transannotationdata()
      annotation_database<-annotation_database()
      KEGG.ID_m<-data.frame(node=metabannotationdata$KEGG.ID,KEGG.ID=metabannotationdata$KEGG.ID)
      KEGG.ID_t<-data.frame(node=transannotationdata$KEGG.ID,KEGG.ID=transannotationdata$KEGG.ID)
      metabdatabse<-annotation_database[[1]]
      transdatabse<-annotation_database[[2]]
      result_m<-run_annotation_count(selectCid_data=KEGG.ID_m, database=metabdatabse,omics_type="metab")
      result_t<-run_annotation_count(selectCid_data=KEGG.ID_t, database=transdatabse,omics_type="trans")
      annotation_results<-list(result_m,result_t)
      print("annotation_results done")
      return(annotation_results)
    })
  }
})


#' @title Prepare Dot Plot Data
#' @description Processes annotation results for visualization.
#' Categorizes pathways into:
#' 1. Shared (Gene-Metabolite co-annotated).
#' 2. Metabolite-only.
#' 3. Gene-only.
#' Selects top 20 pathways for each category.
#' @return A list containing processed plot data, combined annotation table, and pathway category lists.
dotplot_data<-eventReactive(c(input$start_annotation), {
  req(annotation_results())
  withProgress(message = "Processing data...",value=0.8,{
    print("dotplot_data start")
    annotation_results<-annotation_results()
    annotation_m<-annotation_results[[1]]
    annotation_t<-annotation_results[[2]]
    if(nrow(annotation_m) > 0 & nrow(annotation_t) > 0){
      annot_metab<-annotation_m %>%
        arrange(desc(Count))
      annot_ko<-annotation_t %>%
        arrange(desc(Count))
      annot_ko$Types<- "Gene"
      annot_metab$Types<- "Metabolite"
      all_annot<- rbind(annot_metab,annot_ko)
      
      same_path<- intersect(annot_metab$Pathway,annot_ko$Pathway)
      unique_metab <- setdiff(annot_metab$Pathway, annot_ko$Pathway)
      unique_ko <- setdiff(annot_ko$Pathway, annot_metab$Pathway)
      maplist<-list(same_path=same_path,unique_metab=unique_metab,unique_ko=unique_ko)
      process_annot<-function(pathtype){
        topnumber=20
        annot_metab<- annot_metab[which(annot_metab$Pathway %in% pathtype),]
        if(nrow(annot_metab) > topnumber){
          annot_metab <- annot_metab[1:topnumber,]
        }else{
          annot_metab <- annot_metab
        }
        annot_ko<- annot_ko[which(annot_ko$PathwayID %in% annot_metab$PathwayID),]
        plotdata<- rbind(annot_metab,annot_ko) 
        plotdata <- plotdata[complete.cases(plotdata), ]
        return(plotdata)
      }
      plotdata<- list(process_annot(same_path),process_annot(unique_metab),process_annot(unique_ko))
      dotplot_data<-list(plotdata,all_annot,maplist)
      
      return(dotplot_data)
      print("dotplot_data done")
      
      
    }else{
      showNotification("The annotation table is empty.")
    }
  })
})

#table--------------------------------

#' @title Display Annotation Statistics
#' @description Renders a summary table showing the number of annotated pathways and the range of feature counts.
output$annotation_info <- renderTable({
  dotplot_data<-dotplot_data()
  data <- dotplot_data[[2]]
  if(nrow(data)>0){
    data.frame(number_of_annotated_pathways=length(unique(data$Pathway)),
               Count_max=max(data$Count),
               Count_min=min(data$Count))
  }
  
  
  
})

##annotation plot
#' @title Prepare Pathway Map Data
#' @description Prepares data for coloring the KEGG pathway map nodes.
#' Identifies which metabolites and genes from the dataset are present in the selected pathway.
#' extracting their Fold Change values (if differential) or presence status (if pattern-based).
#' @return A list containing mapping vectors (ID -> Value) and details of mapped features.
annotation_plotdata<-eventReactive(c(input$start_annotation,input$Pathway_select), {
  req(dotplot_data())
  req(input$Pathway_select)
  withProgress(message = "Processing data...",value=0.8,{
    print("annotation_plotdata start")

    Pathwaydata<-dotplot_data()[[2]]
    metabannotationdata<-metabannotationdata()
    transannotationdata<-transannotationdata()
    if (!is.null(input$Pathway_select) && input$Pathway_select!="") {
      Pathwayselect<-input$Pathway_select
    }else{
      Pathwayselect<-sort(Pathwaydata$Pathway)[1]
    }
    Pathwaydata_m<-Pathwaydata %>%
      filter(Pathway %in% Pathwayselect)  %>%
      filter(Types=="Metabolite")
    Pathwaydata_t<-Pathwaydata %>%
      filter(Pathway %in% Pathwayselect)  %>%
      filter(Types=="Gene")
    
    if(nrow(Pathwaydata_t)>0){
      pathway_id<-Pathwaydata_t$PathwayID
      
      keggt<-strsplit(as.character(Pathwaydata_t$`KEGG IDs`), "\\+")[[1]]
      t_f<-transannotationdata %>%
        filter(KEGG.ID %in% keggt)
      t_gene<-t_f$gene
    }else{
      t_f<-NA
      t_gene<-NA
    }
    
    if(nrow(Pathwaydata_m)>0){
      pathway_id<-Pathwaydata_m$PathwayID
      
      keggm<-strsplit(as.character(Pathwaydata_m$`KEGG IDs`), "\\+")[[1]]
      m_f<-metabannotationdata %>%
        filter(KEGG.ID %in% keggm)
      m_mz<-m_f$metabolite
    }else{
      m_f<-NA
      m_mz<-NA
    }
    
    if(input$functional_data=="Differential features"){
      if(nrow(Pathwaydata_m)>0){
        m_fplot<- setNames(m_f$`log2(Fold Change)`, m_f$KEGG.ID)
      }else{
        m_fplot<-NA
      }
      if(nrow(Pathwaydata_t)>0){
        t_fplot<- setNames(t_f$`log2(Fold Change)`, t_f$KEGG.ID) 
      }else{
        t_fplot<-NA
      }
    }else{
      if(nrow(Pathwaydata_m)>0){
        
        m_fplot<- setNames(rep(1,length(m_f$KEGG.ID)), m_f$KEGG.ID)
      }else{
        m_fplot<-NA
      }
      if(nrow(Pathwaydata_t)>0){
        
        t_fplot<- setNames(rep(1,length(t_f$KEGG.ID)), t_f$KEGG.ID) 
      }else{
        t_fplot<-NA
      }
    }
    
    
    annotation_plotdata<-list(m_fplot,t_fplot,m_f,t_f,pathway_id,m_mz,t_gene)
    print("annotation_plotdata done")
    return(annotation_plotdata)
  })
  
})

temp_pathwayannotation_file <- reactiveVal(NULL)

#' @title Generate KEGG Pathway Map
#' @description Draws the KEGG pathway map with colored nodes.
#' 1. Loads the base pathway image and coordinate configuration (.conf).
#' 2. Draws circles for Metabolites (colored by FC: Orange=Up, Blue=Down).
#' 3. Draws rectangles for Genes (colored by FC: Red=Up, Green=Down).
#' @return The composite image object (magick image).
annotation_plotsave <- eventReactive(c(input$start_annotation, input$Pathway_select), {
  req(annotation_plotdata())
  req(input$Pathway_select)
  req(input$speciesname_select)
  
  withProgress(message = "Processing data...", value = 0.8, {
    annotation_plotdata <- annotation_plotdata()
    
    m_fplot <- annotation_plotdata[[1]]
    t_fplot <- annotation_plotdata[[2]]
    pathway_id <- annotation_plotdata[[5]]
    species <- input$speciesname_select
    
    local({
      devs <- dev.list()
      if (!is.null(devs)) {
        for (d in devs) {
          tryCatch({
            dev.off(d)
          }, error = function(e) {})
        }
      }
      
      mappath <- "./source/Database/KEGG/map"
      img <- image_read(file.path(mappath, paste0(pathway_id, ".png")))
      img <- image_convert(img, colorspace = 'gray')
      
      lines <- suppressWarnings(readLines(file.path(mappath, paste0(pathway_id, ".conf"))))
      dim <- image_info(img)
      width <- dim$width
      height <- dim$height

      canvas <- image_blank(width, height, "none")

      img_draw <- image_draw(canvas)
      on.exit(dev.off(), add = TRUE)
      
      circle <- function(x, y, r, border = "black", lty = "solid", lwd = 1) {
        angle <- seq(0, 2 * pi, length.out = 100)
        x_circle <- x + r * cos(angle)
        y_circle <- y + r * sin(angle)
        lines(x_circle, y_circle, col = border, lty = lty, lwd = lwd)
      }

      for (line in lines) {
        line <- trimws(line)
        
        if (grepl("^circ", line) && any(!is.na(m_fplot))) {
          coords <- as.numeric(str_extract_all(line, "\\d+")[[1]])
          id_match <- str_match(line, "(C\\d{5})")
          matched_ids <- intersect(id_match, names(m_fplot))
          if (length(matched_ids) > 0 && length(coords) >= 3) {
            x <- coords[1]
            y <- coords[2]
            r <- coords[3]
            
            values <- unlist(m_fplot[matched_ids])
            median_value <- median(values)
            
            if (!is.na(median_value)) {
              color <- ifelse(median_value > 0, "orange",
                              ifelse(median_value < 0, "blue", "gray"))
              circle(x, y, r, border = color, lty = "solid", lwd = 3)
            }
          }
        } else if (grepl("^rect", line) && any(!is.na(t_fplot))) {      
          coords <- as.numeric(str_extract_all(line, "\\d+")[[1]])
          ids <- str_extract_all(line, "K\\d{5}")[[1]]
          matched_ids <- intersect(ids, names(t_fplot))
          
          if (length(matched_ids) > 0 && length(coords) >= 4) {
            xleft <- coords[1]
            ybottom <- coords[2]
            xright <- coords[3]
            ytop <- coords[4]
            
            values <- unlist(t_fplot[matched_ids])
            median_value <- median(values)
            if (!is.na(median_value)) {
              color <- ifelse(median_value > 0, "red",
                              ifelse(median_value < 0, "green", "gray"))
              rect(xleft, ybottom, xright, ytop, border = color, lty = "solid", lwd = 3)
            }
          }
        }
      }

      img <- image_composite(img, img_draw, offset = "+0+0")
      
      temp_dir <- tempdir()
      output_file <- file.path(temp_dir, paste0(pathway_id, ".new.png"))
      
      image_write(img, path = output_file, format = "png")
      temp_pathwayannotation_file(output_file)
      
      return(img)
    })
  })
})

#' @title Render KEGG Pathway Map
#' @description Displays the generated KEGG pathway image in the UI.
observeEvent(c(input$start_annotation, input$Pathway_select), {
  output$pathway_annotation_plot <- renderImage({
    req(annotation_plotsave())
    temp_path <- temp_pathwayannotation_file()
    list(src = temp_path)
  }, deleteFile = FALSE)
})



#' @title Display Mapped Metabolites Table
#' @description Shows a table of metabolites mapped to the currently selected pathway,
#' including their statistics (P-value, FC) if available.
output$annotation_namemap_info_m <- renderTable({
  annotation_plotdata<-annotation_plotdata()
  req(annotation_plotdata)
  if(any(!is.na(annotation_plotdata[[3]]))){
    data<-annotation_plotdata[[3]]
  }else{
    data<-NA
  }
  if(any(!is.na(data))){
    if("State" %in% colnames(data)){
      d= data.frame(metabolite=data$metabolite,
                    KEGG.ID=data$KEGG.ID,
                    State=data$State,
                    p_val_adj=data$p_val_adj,
                    `log2(Fold Change)`=data$`log2(Fold Change)`
      )
    }else{
      d=data.frame(metabolite=data$metabolite,
                   KEGG.ID=data$KEGG.ID
                   
      )
    }
    if("mz" %in% colnames(data) && "Name" %in% colnames(data)){
      d$Name=data$Name
    }
    d
  }
  
})
#' @title Display Mapped Genes Table
#' @description Shows a table of genes mapped to the currently selected pathway.
output$annotation_namemap_info_t <- renderTable({
  annotation_plotdata<-annotation_plotdata()
  req(annotation_plotdata)
  if(any(!is.na(annotation_plotdata[[4]]))){
    data<-annotation_plotdata[[4]]
  }else{
    data<-NA
  }
  if(any(!is.na(data))){
    if("State" %in% colnames(data)){
      data.frame(gene=data$gene,
                 KEGG.ID=data$KEGG.ID,
                 State=data$State,
                 p_val_adj=data$p_val_adj,
                 `log2(Fold Change)`=data$`log2(Fold Change)`
      )
    }else{
      data.frame(gene=data$gene,
                 KEGG.ID=data$KEGG.ID
      )
    }
  }
})

#' @title Prepare Feature Spatial Data
#' @description Extracts spatial intensity data for a specific metabolite or gene selected from the pathway map.
#' Used to verify the spatial distribution of interesting features found in the pathway.
pathway_annotation_ion_data<- eventReactive(c(input$start_annotation,input$Pathway_select,input$Pathway_m_select,input$Pathway_t_select), {#
  req(data_rds())
  req(annotation_plotdata())
  withProgress(message = "Processing data...",value=0.8,{
    data_rds <- data_rds()
    #m
    data_m=data_rds[[1]]
    plotdata_m <- data_m@meta.data
    if(!is.null(input$Pathway_m_select) && input$Pathway_m_select!=""){
      ion <- as.character(input$Pathway_m_select)
    }else{
      ion <- NA
    }
    if(!is.na(ion)){
      plotdata_m$intensity<- data_m@assays$Spatial$counts[ion,]
      plotdata_m$norm_intensity <-100*(plotdata_m$intensity)/max(plotdata_m$intensity)
    }else{
      plotdata_m<-NA
    }
    
    ####
    data_t=data_rds[[2]]
    plotdata_t <- data_t@meta.data
    if(!is.null(input$Pathway_t_select) && input$Pathway_t_select!=""){
      ion_t <- as.character(input$Pathway_t_select)
    }else{
      ion_t <-NA 
    }
    
    if(!is.na(ion_t)){
      plotdata_t$intensity<- data_t@assays$Spatial$counts[ion_t,]
      plotdata_t$norm_intensity <-100*(plotdata_t$intensity)/max(plotdata_t$intensity)
    }else{
      plotdata_t<-NA
    }
    
    if (exists("plotdata_t") && exists("plotdata_m")) {
      pathway_annotation_ion_data<-list(plotdata_m,plotdata_t)
      return(pathway_annotation_ion_data)
    }else{
      return(NULL)
    }
  })
})
#' @title Generate Metabolite Spatial Plot (Pathway)
#' @description Creates a spatial heatmap for the selected pathway metabolite.
pathway_annotation_m_plot_save <- eventReactive(c(input$start_annotation,input$Pathway_select,input$Pathway_m_select), {#
  req(pathway_annotation_ion_data()[[1]])
  
  plotdatalist <- pathway_annotation_ion_data()
  data<-plotdatalist[[1]]
  if(any(!is.na(data))){
    heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
    
    p <- ggplot(data, aes(x = x, y = y)) +
      geom_point(aes(color = norm_intensity), size = 1) +
      scale_color_gradientn(colours = heatmap_Palette(100)) +
      xlim(min(data$x),max(data$x))+
      ylim(min(data$y),max(data$y))+
      coord_equal() +
      xlab("")+ylab("")+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank()) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
  }else{
    p<-NULL
  }
  
  return(p)
  
})
#' @title Generate Gene Spatial Plot (Pathway)
#' @description Creates a spatial heatmap for the selected pathway gene.
pathway_annotation_t_plot_save <- eventReactive(c(input$start_annotation,input$Pathway_select,input$Pathway_t_select), {#
  req(pathway_annotation_ion_data()[[2]])
  
  plotdatalist <- pathway_annotation_ion_data()
  data<-plotdatalist[[2]]
  if(any(!is.na(data))){
    heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
    p <- ggplot(data, aes(x = x, y = y)) +
      geom_point(aes(color = norm_intensity), size = 1) +
      scale_color_gradientn(colours = heatmap_Palette(100)) +
      xlim(min(data$x),max(data$x))+
      ylim(min(data$y),max(data$y))+
      coord_equal() +
      xlab("")+ylab("")+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank()) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
  }else{
    p<-NULL
  }
  return(p)
})

#' @title Render Metabolite Spatial Plot (Pathway)
#' @description Displays the spatial plot for the selected metabolite.
output$pathway_annotation_m_plot <- renderPlot({
  if(any(!is.null(pathway_annotation_m_plot_save()))){
    withProgress(message = "Plotting...",value=0.8,{
      pathway_annotation_m_plot_save()
    })
  }else{
    NULL
  }
})

#' @title Render Gene Spatial Plot (Pathway)
#' @description Displays the spatial plot for the selected gene.
output$pathway_annotation_t_plot <- renderPlot({
  if(any(!is.null(pathway_annotation_t_plot_save()))){
    withProgress(message = "Plotting...",value=0.8,{
      pathway_annotation_t_plot_save()
    })
  }else{
    NULL
  }
})


##dotplot
#' @title Filter Dot Plot Data
#' @description Selects the subset of pathway results (Co-annotated, Metab-only, or Gene-only)
#' to be displayed in the dot plot based on user choice.
dotplot_data_save<-eventReactive(c(input$start_annotation,input$bubble_pathway_types), {
  req(input$bubble_pathway_types)
  req(dotplot_data()[[1]])
  
  if(input$bubble_pathway_types=="Gene-metabolite co-annotated pathways"){
    df<-dotplot_data()[[1]][[1]]
  }else if(input$bubble_pathway_types=="Metabolite-annotated pathways"){
    df<-dotplot_data()[[1]][[2]]
  }else if(input$bubble_pathway_types=="Gene-annotated pathways"){
    df<-dotplot_data()[[1]][[3]]
  }else{
    df<-NULL
    message("No annotation results!")
  }
  return(df)
})
#' @title Generate Dot Plot
#' @description Creates the pathway annotation dot plot (bubble chart) using `AnnotationDotPlot`.
#' Visualizes Pathway Name vs Rich Factor/Count, sized by Count, colored by P-value.
dotplot_save<-eventReactive(c(input$start_annotation,input$bubble_pathway_types), {
  req(dotplot_data_save())
  df<-dotplot_data_save()
  
  source("./source/Pathway/AnnotationDotPlot.R")
  if(any(!is.na(df))){
    
    plot<-AnnotationDotPlot(df)
  }else{
    plot<-NULL
  }
  return(plot)
})
#' @title Render Dot Plot
#' @description Displays the annotation dot plot in the UI.
observeEvent(c(input$start_annotation,input$bubble_pathway_types),{
  req(dotplot_save())
  output$dotplot_run <- renderPlot({
    withProgress(message = "Plotting...",value=0.8,{
      devs <- dev.list()
      current_dev <- dev.cur()
      if (!is.null(devs) && length(devs) > 1) {
        for (d in names(devs)[names(devs) != names(current_dev)]) {
          tryCatch({
            dev.off(dev.list()[[d]])
          }, error = function(e) {})
        }
      }
      plot_obj<- dotplot_save()
      
      grid::grid.newpage()
      grid::grid.draw(plot_obj)
      
    })
  })
})

#overlap

#' @title Prepare Venn Diagram Data
#' @description Extracts pathway lists from Metabolomics and Transcriptomics results
#' to visualize the overlap of annotated pathways.
vennplot_data<-eventReactive(c(input$start_annotation), {
  req(annotation_results())
  annotation_results<-annotation_results()
  annotation_m<-annotation_results[[1]]
  annotation_t<-annotation_results[[2]]
  venn_list <- list(Metabolite = annotation_m$Pathway, Gene = annotation_t$Pathway)
  vennplot_data<-venn_list
  return(vennplot_data)
})
temp_venn_plot_file <- reactiveVal(NULL)
#' @title Render Venn Diagram
#' @description Displays the Venn diagram showing shared and unique pathways between omics.
output$venn_plot <- renderPlot({
  req(vennplot_data())
  vennplot_data<-vennplot_data()
  
  withProgress(message = "Plotting...",value=0.8,{
    temp_file <- tempfile(fileext = ".png")
    
    devs <- dev.list()
    current_dev <- dev.cur()
    if (!is.null(devs) && length(devs) > 1) {
      for (d in names(devs)[names(devs) != names(current_dev)]) {
        tryCatch({
          dev.off(dev.list()[[d]])
        }, error = function(e) {})
      }
    }
    
    p<- venn.diagram(vennplot_data, filename = temp_file,
                     height = 600, width = 600,
                     resolution =300,
                     imagetype="png",
                     col = "transparent",
                     fill = c('yellow', 'skyblue'),
                     alpha = 0.5, cex = 0.5,
                     fontfamily = "serif",
                     fontface = "bold",
                     cat.cex = 0.5,cat.pos = 0,
                     cat.fontfamily = "serif",
                     rotation.degree = 0)
    
    print("venn done")
    temp_venn_plot_file(temp_file)
    grid::grid.raster(png::readPNG(temp_file))
    
  })
})




#' @title Dynamic File Upload UI
#' @description Renders file input buttons for custom annotation tables (Metabolomics/Transcriptomics)
#' and download buttons for demo data.
output$annotationfile_button_container <- renderUI({
  if (input$annotation_select=="upload your annotation data") {
    tagList(
      p("The metabolic annotation file (.txt) must contain the following columns: metabolite, KEGG.ID. The gene annotation file must contain the following columns: gene, KEGG.ID. You can click on the download button below to get demo annotation files"),
      downloadButton("download_metab_pathway_annotation_demo", "Download demo metabolic annotation file"),
      downloadButton("download_Gene_pathway_annotation_demo", "Download demo gene annotation file"),
      fileInput("metabannotationfile", "Upload the spatial metabolomics annotation file",
                accept = ".txt"),
      fileInput("transannotationfile", "Upload the spatial transcriptomics annotation file",
                accept = ".txt")
    )
    
  }else {
    NULL
  }
  
})

#' @title Update Gene Selector
#' @description Populates the gene selection dropdown with genes present in the currently displayed pathway.
observeEvent(annotation_plotdata(), {
  annotation_plotdata<-annotation_plotdata()
  t_gene<-as.character(annotation_plotdata[[7]])
  updateSelectizeInput(
    session = getDefaultReactiveDomain(),
    inputId = "Pathway_t_select",
    choices =t_gene,
    server = TRUE,
    selected = if (length(t_gene) > 0) t_gene[1] else NULL,
    options = list(
      maxOptions = 1000,          
      searchConjunction = 'and',   
      render = I("{                
      option: function(item, escape) {
        return '<div>' + escape(item.label) + '</div>';
      }
    }"),
      loadThrottle = 300,          
      placeholder = "Search by feature name...",
      closeAfterSelect = TRUE      
    )
  )
  
})
#' @title Update Metabolite Selector
#' @description Populates the metabolite selection dropdown with metabolites present in the currently displayed pathway.
observeEvent(annotation_plotdata(), {
  
  annotation_plotdata<-annotation_plotdata()
  m_mz<-as.character(annotation_plotdata[[6]])
  updateSelectizeInput(
    session = getDefaultReactiveDomain(),
    inputId = "Pathway_m_select",
    choices =m_mz,
    server = TRUE,
    selected = if (length(m_mz) > 0) m_mz[1] else NULL,
    options = list(
      maxOptions = 1000,          
      searchConjunction = 'and',   
      render = I("{                
      option: function(item, escape) {
        return '<div>' + escape(item.label) + '</div>';
      }
    }"),
      loadThrottle = 300,          
      placeholder = "Search by feature name...",
      closeAfterSelect = TRUE      
    )
  )
  
})


#' @title Dynamic Analysis Type UI
#' @description Allows user to choose between "Differential features" and "Pattern-specific features"
#' based on which upstream analyses have been completed.
output$functional_data_button_container <- renderUI({
  if(!is.null(diff_omics()) && !is.null(spatial_pattern_ionlist())){
    selectizeInput("functional_data", "Select data for functional association analysis:", choices = c("Differential features","Pattern-specific features"), selected = "Differential features")
  }else if(!is.null(diff_omics()) && is.null(spatial_pattern_ionlist())){
    selectizeInput("functional_data", "Select data for functional association analysis:", choices = c("Differential features"), selected = "Differential features")
    
  }else if(is.null(diff_omics()) && !is.null(spatial_pattern_ionlist())){
    selectizeInput("functional_data", "Select data for functional association analysis:", choices = c("Pattern-specific features"), selected = "Pattern-specific features")
    
  }else{
    p("Note: Please complete differential analysis or spatial pattern analysis first.")
  }
})
####pathway type
#' @title Update Pathway Selector
#' @description Updates the list of available pathways for the map view based on the selected category (Shared, Metab-only, Gene-only).
observeEvent(c(input$start_annotation,input$Pathway_types), {
  dotplot_data<-dotplot_data()
  maplist<-dotplot_data[[3]]
  req(maplist)
  str(maplist)
  if(input$Pathway_types=="gm"){
    Pathway<-maplist$same_path
  }else if(input$Pathway_types=="onlym"){
    Pathway<-maplist$unique_metab
  }else{
    Pathway<-maplist$unique_ko
  }
  Pathway<-sort(Pathway)
  updateSelectInput(session = getDefaultReactiveDomain(), "Pathway_select", choices = Pathway, 
                    selected = if (length(Pathway) > 0) Pathway[1] else NULL)
  
  
})


#' @title Update Bubble Plot Type
#' @description Updates the available bubble plot categories based on the results (e.g., if only metabolites annotated, only show that option).
observeEvent(c(input$start_annotation), {
  dotplot_data<-dotplot_data()
  datasets<-dotplot_data[[1]]
  namec=c("Gene-metabolite co-annotated pathways","Metabolite-annotated pathways","Gene-annotated pathways")
  x <- character(0)
  for (i in 1:3) {
    dataset <- datasets[[i]]
    name<-namec[i]
    
    if (!is.null(dataset) && 
        any(!is.na(dataset))) {
      
      if (nrow(dataset) > 0 && ncol(dataset) > 0) {
        x <- c(x, namec[i])
      }
    }
  }
  updateSelectInput(session = getDefaultReactiveDomain(), "bubble_pathway_types", choices = x, 
                    selected = if (length(x) > 0) x[1] else NULL)
  
})

#' @title Dynamic Direction Filter UI
#' @description Renders options to filter differential features by direction (Up/Down/Both).
output$functional_diffdata_button_container <- renderUI({
  if(!is.null(input$functional_data)){
    if(input$functional_data=="Differential features"){
      selectizeInput("functional_diffdata", "Select up/down-regulated differential features for analysis",
                     choices = c("Up-regulated features","Down-regulated features","up-regulated and down-regulated features"), selected = "up-regulated and down-regulated features")
      
    }else{
      NULL
    }
  }
  
})


