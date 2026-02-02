#' @title Create Sankey Diagram
#' @description Generates a Sankey diagram to visualize flows or relationships between nodes.
#' Uses the `networkD3` package.
#' @param links Data frame containing 'source', 'target', 'value', and optionally 'color' columns.
#' @return A Sankey diagram widget (networkD3 object).
sankyplot_ID = function(links){
  # Create nodes dataframe from unique source and target values
  nodes <- data.frame(
    name=c(as.character(links$source),
           as.character(links$target)) %>% unique()
  )
  
  # Map source nodes to IDs (starting from 0)
  links$IDsource <- match(links$source, nodes$name)-1
  
  # Map target nodes to IDs (starting from 0)
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # Create Sankey diagram using networkD3 package
  p<-sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",
                   LinkGroup = "color", 
                   sinksRight=FALSE,nodeWidth = 10, 
                   nodePadding = 4) 
  
  # Return the plot object
  return(p)
}