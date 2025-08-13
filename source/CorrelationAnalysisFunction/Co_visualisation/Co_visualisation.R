spatial_coord_processing<-function(combined_matrix,meta.data,rescale=TRUE){

  decon_mtrx = t(combined_matrix)
  if(rescale){
    decon_mtrx<-as.data.frame(apply(decon_mtrx, 2, function(x) {
      x/max(x)
      }))
  }else{
    decon_mtrx<-as.data.frame(decon_mtrx)
  }

 
  ion_types_all <- colnames(decon_mtrx)
  decon_mtrx$x<-meta.data$x
  decon_mtrx$y<-meta.data$y
  spatial_coord <- decon_mtrx %>%
    dplyr::select(x,y,all_of(ion_types_all))
 
  return(spatial_coord)
}


multiple_ions_plot<-function(spatial_coord,pair){
  pt.size=1
  pairname=c("feature1","feature2","feature3")
  LRpair =pair[!is.na(pair)]
  data = spatial_coord[,c('x','y',LRpair)] 
  if(!is.na(pair[1])){
    f1=data[,pair[1]]
  }else{
    f1=0
  }
  if(!is.na(pair[2])){
    f2=data[,pair[2]]
  }else{
    f2=0
  }
  if(!is.na(pair[3])){
    f3=data[,pair[3]]
  }else{
    f3=0
  }


  data$color <- rgb(f1,f2,f3, maxColorValue = 1)

  
  p1<-ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = color), size = 1) +
    xlim(min(data$x), max(data$x)*1.3) +
    ylim(min(data$y), max(data$y)) +
    coord_equal() +
    scale_color_identity() +  
  xlab("")+ylab("")+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank()) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())
  if(!is.na(pair[1])){

  p1<-p1 + annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)/2, ymax = max(data$y)/2+4, 
           fill = "red", color = "black") +  
    annotate("text", x = max(data$x)+ 5, y = max(data$y)/2+2, label = pairname[1], 
             hjust = 0, vjust = 0.5, size = 3, color = "black")   
    
  }
  if(!is.na(pair[2])){

    p1<-p1 +annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)*1.2/2, ymax = (max(data$y)*1.2/2)+4, 
           fill = "green", color = "black") +  
    annotate("text", x = max(data$x)+ 5, y = (max(data$y)*1.2/2)+2, label = pairname[2], 
             hjust = 0, vjust = 0.5, size = 3, color = "black")   
  }
  if(!is.na(pair[3])){

    p1<-p1 +annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)*1.4/2, ymax = (max(data$y)*1.4/2)+4, 
                     fill = "blue", color = "black") +  
      annotate("text", x = max(data$x)+ 5, y =  (max(data$y)*1.4/2)+2, label =pairname[3], 
               hjust = 0, vjust = 0.5, size = 3, color = "black")  
  }
  result<-list(p1,data)
  return(result)
}



Co_visualisation_plot<-function(spatial_coord,pair){
  pairname=c("feature1","feature2","feature3")
  f1_col=NULL
  f2_col=NULL
  f3_col=NULL
  pt.size=1
  Fpair =pair[!is.na(pair)]
  
  location = spatial_coord[,c('x','y')]
  topn=floor(0.2*dim(location)[1])
  expr =spatial_coord[,Fpair,drop=FALSE]
  ncell<-dim(expr)[1]
  
  expstatus<-function(f,fname){
    f_q25 = quantile(f, probs = 0.25)
    f_q75 = quantile(f, probs = 0.75)
    n1 = which(f > f_q75)
    n2 = which(f < f_q25)
    f_col<-rep(paste0(fname,"_medium"),ncell)
    f_col[n1]<-paste0(fname,"_high")
    f_col[n2]<-paste0(fname,"_low")
    return(f_col)
  }
  if(!is.na(pair[1])){
    f1<-expr[,pair[1]]
    f1_col<-expstatus(f1,fname=pairname[1])
  }
  if(!is.na(pair[2])){
    f2<-expr[,pair[2]]
    f2_col<-expstatus(f2,fname=pairname[2])
  }
  if(!is.na(pair[3])){
    f3<-expr[,pair[3]]
    f3_col<-expstatus(f3,fname=pairname[3])
  }
  
  
  expcol <- paste0(
    ifelse(nzchar(f1_col), paste0(f1_col, "_"), ""),
    ifelse(nzchar(f2_col), paste0(f2_col, "_"), ""),
    ifelse(nzchar(f3_col), f3_col, "")
  )

  expcol <- gsub("_$", "", expcol)
  
  tmp<-data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol))
  
  p1<-ggplot(tmp,aes(x=x,y=y,col=Exp))+
    geom_point(size=pt.size)+
    labs(color="") + 
    xlim(min(tmp$x),max(tmp$x))+
    ylim(min(tmp$y),max(tmp$y))+
    coord_equal()+
    xlab("")+ylab("")+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank()) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())
  
  data<-list(p1,tmp)
  return(data)
}
