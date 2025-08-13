
#ID
sankyplot_ID=function(links){
  nodes <- data.frame(
    name=c(as.character(links$source),
           as.character(links$target)) %>% unique()
  )
  links$IDsource <- match(links$source, nodes$name)-1
  links$IDtarget <- match(links$target, nodes$name)-1
  p<-sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",
                   LinkGroup = "color", 
                   sinksRight=FALSE,nodeWidth = 10, 
                   nodePadding = 4) 
  return(p)
}



