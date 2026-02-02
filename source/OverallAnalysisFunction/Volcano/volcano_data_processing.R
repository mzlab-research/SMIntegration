#' @title Prepare Data for Volcano Plot
#' @description Formats differential analysis results for volcano plotting.
#' Handles numeric conversions, factor levels for 'State' (Up/Down/NS), and optional VIP scores.
#' @param data Data frame containing differential analysis results.
#' @return A list containing:
#'   - [1] Processed data frame ready for plotting.
#'   - [2] Simplified data frame for export/download.
volcano_data_processing<-function(data){
  Pvalue_type="p_val_adj"
  
  # Ensure numeric types
  data$log_test<-as.numeric(data$log_test)
  data$`log2(Fold Change)` <-as.numeric(data$`log2(Fold Change)`)
  
  # Set factor levels for consistent plotting order/colors
  data$State <- factor(data$State, levels = c('Down', 'Non-significant', 'Up')) ##New add
  
  # Handle VIP scores if present (e.g., for PLS-DA results)
  if("VIP" %in% toupper(colnames(data))){
    
    vipsize<-1
    vipmax<-paste0("VIP>=",vipsize)
    vipmin<-paste0("VIP<",vipsize)
    data$VIP<-as.numeric(data$VIP)
    data$VIP_state<-ifelse(data$VIP >= vipsize, vipmax, vipmin)
    data$VIP_state <- factor(data$VIP_state, levels = unique(data$VIP_state))
  }
  
  # Create simplified export version
  Volcano_data_temp<-data
  Volcano_data_temp<-Volcano_data_temp %>%
    dplyr::rename(log_p_val_adj=log_test)  %>%
    dplyr::select(-avg_log2FC,-p_val)
  
  data$absfc <- abs(data$`log2(Fold Change)`)
  result<-list(data,Volcano_data_temp)
 return(result)
}

#' @title Generate Volcano Plot
#' @description Creates a volcano plot using ggplot2.
#' Supports standard Differential Analysis (P-value vs LogFC) and OPLS-DA (with VIP shape).
#' @param data Data frame returned by volcano_data_processing().
#' @param type Character. Plot title/type label.
#' @param group Vector. Comparison groups.
#' @param pvalue Numeric. P-value threshold for horizontal line.
#' @param FC_Threshold Numeric. Fold change threshold for vertical lines.
#' @return A ggplot object.
volcano_plot_processing<-function(data,type,group,pvalue,FC_Threshold){
  if(nrow(data) == 0) return(NULL) # Return NULL if data is empty

  Pvalue_type="p_val_adj"
  comparename=paste0(group[1],":",group[2])
  xlims <- ceiling(max(data$absfc))
  State_value <- unique(data$State)
  State_len <- length(State_value)
  
  ### Color Arguments
  # Robust color mapping using named vector
  # ggplot2 scale_color_manual handles named vectors automatically
  # Ensures 'Up' is always HotPink, 'Down' is LightSkyBlue, regardless of presence
  scale_color_map <- c("Down" = "LightSkyBlue", 
                       "Non-significant" = "grey", 
                       "Up" = "HotPink")
  
    # Plotting Logic
    if("VIP" %in% toupper(colnames(data))){
    # OPLS-DA Style: Shape maps to VIP score
    p<- data %>%
      ggplot(aes(`log2(Fold Change)`,log_test))+
      theme_classic()+
      labs(title=paste(type))+
      geom_point(alpha= I(1/2),aes(color = State,shape = VIP_state),size = 2.5)+
      scale_color_manual(values = scale_color_map)+
      scale_shape_manual(values = c(17, 16))+
      geom_hline(yintercept = -log10(pvalue),linetype=6,size = .3,color = "black")+
      geom_vline(xintercept=c(-log2(FC_Threshold),log2(FC_Threshold)),linetype=6,size = .3,color = "black")+
      scale_fill_discrete(label = c("VIP < 1","VIP >= 1"))+
      xlim(-xlims,xlims)+
      xlab("log2(Fold Change)")+
      ylab(paste("-log10(",Pvalue_type,")",sep = ""))+
      theme(panel.grid=element_blank())+
      theme(axis.line.x = element_line(color="black", size = 0.5),
            axis.line.y = element_line(color="black", size = 0.5),
            plot.title = element_text(hjust = 0.5))

    
  }else{
    # Standard Style
    p<- data %>%
      ggplot(aes(`log2(Fold Change)`,log_test))+
      theme_classic()+
      labs(title=paste(type))+
      geom_point(alpha= I(1/2),aes(color = State),size = 2.5)+
      scale_color_manual(values = scale_color_map)+
      scale_shape_manual(values = c(17, 16))+
      geom_hline(yintercept = -log10(pvalue),linetype=6,size = .3,color = "black")+
      geom_vline(xintercept=c(-log2(FC_Threshold),log2(FC_Threshold)),linetype=6,size = .3,color = "black")+
      xlim(-xlims,xlims)+
      xlab("log2(Fold Change)")+
      ylab(paste("-log10(",Pvalue_type,")",sep = ""))+
      theme(panel.grid=element_blank())+
      theme(axis.line.x = element_line(color="black", size = 0.5),
            axis.line.y = element_line(color="black", size = 0.5),
            plot.title = element_text(hjust = 0.5))

  }
return(p)
}