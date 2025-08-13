volcano_data_processing<-function(data){
  Pvalue_type="p_val_adj"
  data$log_test<-as.numeric(data$log_test)
  data$`log2(Fold Change)` <-as.numeric(data$`log2(Fold Change)`)
  data$State <- factor(data$State, levels = c('Down', 'Non-significant', 'Up')) ##New add
  if("VIP" %in% toupper(colnames(data))){
    
    vipsize<-1
    vipmax<-paste0("VIP>=",vipsize)
    vipmin<-paste0("VIP<",vipsize)
    data$VIP<-as.numeric(data$VIP)
    data$VIP_state<-ifelse(data$VIP >= vipsize, vipmax, vipmin)
    data$VIP_state <- factor(data$VIP_state, levels = unique(data$VIP_state))
  }
  
  Volcano_data_temp<-data
  Volcano_data_temp<-Volcano_data_temp %>%
    dplyr::rename(log_p_val_adj=log_test)  %>%
    dplyr::select(-avg_log2FC,-p_val)
  
  data$absfc <- abs(data$`log2(Fold Change)`)
  result<-list(data,Volcano_data_temp)
 return(result)
}
#plot
volcano_plot_processing<-function(data,type,group,pvalue,FC_Threshold){
  Pvalue_type="p_val_adj"
  comparename=paste0(group[1],":",group[2])
  xlims <- ceiling(max(data$absfc))
  State_value <- unique(data$State)
  State_len <- length(State_value)
  ###color args
  if(State_len == 1){
    if(State_value == 'Non-significant'){
      scale_color <- "grey"
    }else if(State_value == "Down"){
      scale_color <- "LightSkyBlue"
    }else if(State_value == "Up"){
      scale_color <- "HotPink"
    }
  }else if(State_len == 2){
    if('Down' %in% State_value & 'Non-significant' %in% State_value){
      scale_color <- c("LightSkyBlue","grey")
    }else if('Non-significant' %in% State_value & 'Up' %in% State_value){
      scale_color <- c("grey","HotPink")
    }else if('Down' %in% State_value & 'Up' %in% State_value){
      scale_color <- c("LightSkyBlue","HotPink")
    }
  }else if(State_len == 3){
    scale_color <- c("LightSkyBlue","grey","HotPink")
  }
    if("VIP" %in% toupper(colnames(data))){
    p<- data %>%
      ggplot(aes(`log2(Fold Change)`,log_test))+
      theme_classic()+
      labs(title=paste(type))+
      geom_point(alpha= I(1/2),aes(color = State,shape = VIP_state),size = 2.5)+
      scale_color_manual(values = scale_color)+
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
    p<- data %>%
      ggplot(aes(`log2(Fold Change)`,log_test))+
      theme_classic()+
      labs(title=paste(type))+
      geom_point(alpha= I(1/2),aes(color = State),size = 2.5)+
      scale_color_manual(values = scale_color)+
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