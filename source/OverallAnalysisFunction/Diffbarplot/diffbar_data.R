#' @title Prepare Differential Count Data for Bar Plot
#' @description Summarizes the number of Up/Down regulated features per comparison group.
#' Formats data for stacked bar plot visualization.
#' @param data Data frame of differential results (with 'State' and 'Sample' columns).
#' @param group Vector. Comparison groups (e.g., c("Treatment", "Control")).
#' @param type Character. Feature type (e.g., "Metabolite").
#' @return A list containing:
#'   - [1] data_count_long: Long format data for ggplot2.
#'   - [2] diff_count: Wide format summary table.
diffbar_data<-function(data,group,type){
  # Handle empty input gracefully
  if(nrow(data) == 0){
      # Return empty structure if input data is empty
      data_count_long <- data.frame(Sample=character(), type=character(), State=character(), value=numeric())
      diff_count <- data.frame(Sample=character(), State=character(), count=numeric(), type=character()) # Minimal structure
      return(list(data_count_long, diff_count))
  }

  compare_group<-paste0(group[1],":",group[2])
  
  # Summarize counts by State (Up/Down) and Sample
  diff_count <- data %>%
    filter(Sample!="") %>%
    group_by(State, Sample) %>%
    dplyr::summarise(count = n()) %>%
    reshape2::dcast(Sample~State) # Pivot to wide format
  
  diff_count$type=type
  diff_count[is.na(diff_count)] <- 0
  
  # Handle missing groups or states (e.g., if only Up regulated found)
  no_diff_compare_group <- setdiff(compare_group, diff_count$Sample)
  names_diff_count <- names(diff_count)
  
  if(!('Up' %in% names_diff_count)){   #only down
    diff_count <- diff_count %>% mutate(Up = 0)
    
  }
  if(!('Down' %in% names_diff_count)){
    diff_count <- diff_count %>% mutate(Down = 0)
  }
  
  # Add total count and handle groups with no differences
  diff_count %>%
    mutate(`Diff of total` = Down + Up) %>%
    dplyr::select(Sample, `Diff of total`, Up, Down) %>%
    dplyr::rename(Group = Sample) %>%
    bind_rows(tibble(Group = no_diff_compare_group,`Diff of total` = 0, Up = 0, Down = 0)) #%>%

  # Convert back to long format for plotting
  data_count_long <- diff_count  %>%
    pivot_longer(cols = c("Down", "Non-significant", "Up"), 
                 names_to = "State",
                 values_to = "value") %>%
    filter(State != "Non-significant") %>%
    dplyr::select(Sample,type,State,value)
  
  result<-list(data_count_long,diff_count)
  return(result)
}



