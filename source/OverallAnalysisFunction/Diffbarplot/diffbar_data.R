diffbar_data<-function(data,group,type){
  compare_group<-paste0(group[1],":",group[2])
  diff_count <- data %>%
    filter(Sample!="") %>%
    group_by(State, Sample) %>%
    dplyr::summarise(count = n()) %>%
    reshape2::dcast(Sample~State)
  diff_count$type=type
  diff_count[is.na(diff_count)] <- 0
  no_diff_compare_group <- setdiff(compare_group, diff_count$Sample)
  names_diff_count <- names(diff_count)
  
  if(!('Up' %in% names_diff_count)){   #only down
    diff_count <- diff_count %>% mutate(Up = 0)
    
  }
  if(!('Down' %in% names_diff_count)){
    diff_count <- diff_count %>% mutate(Down = 0)
  }
  diff_count %>%
    mutate(`Diff of total` = Down + Up) %>%
    dplyr::select(Sample, `Diff of total`, Up, Down) %>%
    dplyr::rename(Group = Sample) %>%
    bind_rows(tibble(Group = no_diff_compare_group,`Diff of total` = 0, Up = 0, Down = 0)) #%>%

  
  data_count_long <- diff_count  %>%
    pivot_longer(cols = c("Down", "Non-significant", "Up"), 
                 names_to = "State",
                 values_to = "value") %>%
    filter(State != "Non-significant") %>%
    dplyr::select(Sample,type,State,value)
  result<-list(data_count_long,diff_count)
  return(result)
}



