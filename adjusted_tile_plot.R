#test different function to display frost risk

library(chillR)
library(tidyverse)

#read frost data
all_frost <- read.csv(paste0("data/predictions/frost_summary/",frost_summary_files[1])) %>%
  rename(Scenario = Time,
         Year = year) %>%
  mutate(End_year = Year,
         Scenario = as.character(Scenario))

for (i in 2:length(frost_summary_files))
{current <- read.csv(paste0("data/predictions/frost_summary/",frost_summary_files[i])) %>%
  rename(Scenario = Time,
         Year = year) %>%
  mutate(End_year = Year,
         Scenario = as.character(Scenario))
all_frost <- bind_rows(all_frost,current)}

all_frost$bloom[all_frost$bloom>200] <- all_frost$bloom[all_frost$bloom>200] -365
all_frost$bloomdate <- as.Date(ISOdate(2000,12,31)+all_frost$bloom*3600*24)
all_frost$bloomdate[which(all_frost$bloomdate=="2000-12-31")]<-NA
all_frost$SSP[all_frost$SSP=="ssp126"] <- "SSP1"
all_frost$SSP[all_frost$SSP=="ssp245"] <- "SSP2"
all_frost$SSP[all_frost$SSP=="ssp370"] <- "SSP3"
all_frost$SSP[all_frost$SSP=="ssp585"] <- "SSP5"


#function to make range plot
#this function may look frightening, but I think it is great
#it can handle a lot of things
#it can add isotopic lines and it can exclude lines that are too close to another 
#or to the border of the plot (and you can control how far apart the labels should be)
#you set the seperation of the color strips flexibly. you can add borders
#you can control the inversion point of the color scale
#you can adjust the labelling of the color legend
#the function automatically formats the input data for the plot
#it assumes that you have facets for the rows and columns, I don't know if it could
#handle missing facets, though
#you can also have no isotopic lines. Or you could add lines dependent on quantiles of the distribution of the baseline 
get_range_plot <- function(df, colours, target_col, group_var,
                           col_facet, row_facet, x_col, legend_name,
                           order_y_axis, h_levels =NULL, h_level_baseline = NULL,
                           critical_value  = NULL,
                           start = 1,
                           end = 100,
                           tile_width = 1,
                           add_border_color_stripes = FALSE,
                           rng = NULL,
                           x_label = NULL,
                           y_label = NULL,
                           fill_breaks = NULL,
                           fill_labels = NULL,
                           d_levels = NULL,
                           d_label = NULL,
                           d_label_text_size = NULL,
                           d_label_mindiff = 2,
                           seg_length = NULL,
                           int_d_levels = NULL, 
                           base_size = 15){
  
  
  df_plot <- prepare_id(df, 
                        target_col = target_col, 
                        group_var = group_var, 
                        order_y_axis = order_y_axis, 
                        h_levels = h_levels,
                        start = start,
                        end = end) %>% 
    pluck('df_plot') %>% 
    as.data.frame()
  
  df_plot$target_col_i <- df_plot[,target_col]
  df_plot$col_facet_i <- df_plot[,col_facet]
  df_plot$row_facet_i <- df_plot[,row_facet]
  df_plot$x_col_i <- df_plot[,x_col]
  
  if(is.factor(df_plot$col_facet_i) == FALSE) df_plot$col_facet_i <- as.factor(df_plot$col_facet_i )
  if(is.factor(df_plot$row_facet_i) == FALSE) df_plot$row_facet_i <- as.factor(df_plot$row_facet_i )
  if(is.factor(df_plot$x_col_i) == FALSE) df_plot$x_col_i <- as.factor(df_plot$x_col_i )
  
  #range of the legend
  if(is.null(rng))   rng <- range(df_plot[,target_col], na.rm = TRUE)
  if(is.null(fill_breaks)) fill_breaks <- waiver()
  if(is.null(fill_labels)) fill_labels <- waiver()
  if(is.null(x_label)) x_label <- waiver()
  if(is.null(y_label)) y_label <- waiver()
  if(is.null(d_label_text_size)) d_label_text_size  <- 3.88
  if(is.null(critical_value)){
    fill_values <- c(0,1)
  } else {
    rng_col <- range(df_plot$target_col_i, na.rm = TRUE)
    mx <- critical_value/rng_col[2]
    fill_values <- c(0,0.1,0.2,0.3,0.6,0.7,0.8,1,1/mx)*mx
  }
  #make the ssp numeric, control the breaks later on
  p <- df_plot %>% 
    mutate(x_col_i = as.integer(as.factor(x_col_i))) %>% 
    ggplot(aes(x = x_col_i, y = id, fill = target_col_i)) + 
    geom_tile(aes(height = step), width = tile_width) +
    facet_grid(row_facet_i~col_facet_i, scales = 'free_x', space = 'free_x') +
    scale_fill_gradientn(colours = colours,
                         values = fill_values,
                         breaks = fill_breaks, 
                         labels = fill_labels,
                         name = legend_name,
                         limits = rng) +
    xlab(x_label) +
    ylab(y_label) + 
    scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0,0)) +
    scale_x_continuous(breaks = 1:length(unique(df_plot$x_col_i)),
                       labels = sort(unique(df_plot$x_col_i)), expand = c(0,0)) +
    theme_bw(base_size = base_size) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if(add_border_color_stripes){
    
    id_unique <- df_plot %>% 
      mutate(x_col_i = as.integer(as.factor(x_col_i)),
             id = paste(x_col_i, row_facet_i, col_facet_i, sep = '--')) %>% 
      pull(id) %>% 
      unique()
    
    x_col_i <- id_unique %>% str_split(pattern = '--') %>% purrr::map_chr(1) %>% as.integer()
    row_facet_i <- id_unique %>% str_split(pattern = '--') %>% purrr::map_chr(2)
    col_facet_i <- id_unique %>% str_split(pattern = '--') %>% purrr::map_chr(3) 
    
    df_hline <- data.frame(x_col_i = x_col_i,
                           row_facet_i = factor(row_facet_i,levels = levels(df_plot$row_facet_i)),
                           col_facet_i = factor(col_facet_i, levels = levels(df_plot$col_facet_i))) %>% 
      mutate(left_border = x_col_i - (tile_width/2),
             right_border = x_col_i + (tile_width/2))
    
    p <- p +
      geom_vline(data = df_hline, aes(xintercept = left_border))+
      geom_vline(data = df_hline, aes(xintercept = right_border))
  }
  
  
  #contour lines based on fixed levels
  if(is.null(d_levels) == FALSE){
    
    if(is.null(d_label)) d_label <- as.character(d_levels)
    #in case we want day ranges instead, find them out
    df_contour2 <- df_plot 
    df_contour2$group_id <-  apply(df_contour2[group_var],MARGIN = 1, function(y) paste0(y, collapse = '_'))
    df_contour2['target_col'] <- df_contour2[target_col]
    d_df <- data.frame(d_level = d_levels, d_label = d_label)
    df_contour2 <- expand.grid(unique(df_contour2$group_id), d_levels) %>% 
      setNames(c('group_id', 'd_level')) %>% 
      merge(df_contour2, by = c('group_id')) %>% 
      mutate(diff = abs(d_level - target_col)) %>% 
      arrange(diff) %>% 
      group_by(group_id, d_level) %>% 
      slice(1) %>% 
      ungroup() %>% 
      as.data.frame() %>% 
      merge(d_df, by = 'd_level')
    df_contour2$target_col_i <- df_contour2[,target_col]
    df_contour2$col_facet_i <- df_contour2[,col_facet]
    df_contour2$row_facet_i <- df_contour2[,row_facet]
    df_contour2$x_col_i <- df_contour2[,x_col]
    
    #make sure that there is a minimum difference between the d_labels
    if(is.null(d_label_mindiff) == FALSE){
      
      #add for each facet a ceiling and bottm. Then calculate the difference
      #if there is a conflict (difference too small) --> remove one label and then re-calculate
      
      df_contour2 <- df_contour2 %>% 
        mutate(extra_row = FALSE) %>% 
        group_by(col_facet_i, row_facet_i, x_col_i) %>% 
        arrange(id) %>% 
        group_modify(~ add_row(.x, id = end, extra_row = TRUE, .after = nrow(.x))) %>% 
        group_modify(~ add_row(.x, id = start, extra_row = TRUE, .before = 1)) %>% 
        ungroup() %>% 
        mutate(iteration_id = paste(col_facet_i, row_facet_i, x_col_i, sep = '--'))
      
      df_contour2_cleaned <- data.frame()
      #iterate over each iteration id
      for(it_id in unique(df_contour2$iteration_id)){
        #it_id <- unique(df_contour2$iteration_id)[9]
        sub <- df_contour2[df_contour2$iteration_id ==it_id,]
        #only run the cleaning if there is at least one entry apart from start and end
        if(nrow(sub)< 3){
          df_contour2_cleaned <- rbind(df_contour2_cleaned, sub)
          next()
        } 
        
        #falg to stop the while loop
        clean_sub <- TRUE
        count <- 1
        while(clean_sub){
          #calculate difference to next row
          diff_vec <- c()
          for(i in 1:(nrow(sub)-1)){
            diff_vec[i] <- abs(sub$id[i] - sub$id[i+1])
          }
          
          #check if any below the threshold for min differences
          if(any(diff_vec < d_label_mindiff) == FALSE){
            clean_sub <- FALSE
          } else {
            #if first or last entry has smallest difference, remove the second or second-last entry
            #otherwise kick the entry before
            i_kick <- which.min(diff_vec)
            if(i_kick == 1){
              sub <- sub[-2, ]
            } else {
              sub <- sub[-i_kick,]
            }
          }
          
          #if there is no entry apart from start or stop, then stop the while loop
          if(nrow(sub) < 3){
            clean_sub <- FALSE
          }
        }
        df_contour2_cleaned <- rbind(df_contour2_cleaned, sub)
      }
      
      df_contour2 <- df_contour2_cleaned %>% 
        filter(extra_row == FALSE) %>% 
        select(-iteration_id, -extra_row)
      
      
    }
    
    
    
    if(is.null(seg_length)) seg_length <- 0.35
    
    p <- p +
      geom_segment(data = df_contour2,
                   aes(x = as.integer(as.factor(x_col_i)) - (tile_width / 2),
                       xend = as.integer(as.factor(x_col_i)) - ((tile_width /2) - seg_length),
                       y = id, yend = id)) +
      geom_segment(data = df_contour2,
                   aes(x = as.integer(as.factor(x_col_i)) +(tile_width/2),
                       xend = as.integer(as.factor(x_col_i)) + ((tile_width/2) - seg_length),
                       y = id, yend = id)) +
      geom_text(data = df_contour2,
                aes(x = as.integer(as.factor(x_col_i)), 
                    y = id, label = d_label), size = d_label_text_size) 
    
    
    
    if(is.null(int_d_levels) == FALSE){
      df_contour2 <- df_plot 
      df_contour2$group_id <-  apply(df_contour2[group_var],MARGIN = 1, function(y) paste0(y, collapse = '_'))
      df_contour2['target_col'] <- df_contour2[target_col]
      d_df <- data.frame(d_level = int_d_levels)
      df_contour2 <- expand.grid(unique(df_contour2$group_id), int_d_levels) %>% 
        setNames(c('group_id', 'd_level')) %>% 
        merge(df_contour2, by = c('group_id')) %>% 
        mutate(diff = abs(d_level - target_col)) %>% 
        arrange(diff) %>% 
        group_by(group_id, d_level) %>% 
        slice(1) %>% 
        ungroup() %>% 
        as.data.frame() %>% 
        merge(d_df, by = 'd_level')
      df_contour2$target_col_i <- df_contour2[,target_col]
      df_contour2$col_facet_i <- df_contour2[,col_facet]
      df_contour2$row_facet_i <- df_contour2[,row_facet]
      df_contour2$x_col_i <- df_contour2[,x_col]
      
      p <- p +
        geom_segment(data = df_contour2,
                     aes(x = as.integer(as.factor(x_col_i)) - (tile_width/2),
                         xend = as.integer(as.factor(x_col_i)) + (tile_width/2),
                         y = id, yend = id), col = 'grey20', linetype = 'dashed') 
      
    }
    
    
    #force that h_levels are not shown when d_levels specified
    h_levels <- NULL
    
    
    #contour line defined on the quantiles of the baseline year
  } else if(is.null(h_levels) == FALSE){
    
    df_contour <- prepare_id(df, 
                             target_col = target_col, 
                             group_var = group_var, 
                             order_y_axis = order_y_axis, 
                             h_levels = h_levels,
                             h_level_baseline = h_level_baseline) %>% 
      pluck('df_contour') %>% 
      as.data.frame()
    
    df_contour$target_col_i <- df_contour[,target_col]
    df_contour$col_facet_i <- df_contour[,col_facet]
    df_contour$row_facet_i <- df_contour[,row_facet]
    df_contour$x_col_i <- df_contour[,x_col]
    
    p <- p +
      geom_segment(data = df_contour,
                   aes(x = as.integer(as.factor(x_col_i)) - (tile_width/2),
                       xend = as.integer(as.factor(x_col_i)) + (tile_width/2),
                       y = id, yend = id)) 
    
    
  }
  return(p)
}

#function called inside the get range plot to prepare the input data
#prepare id per ssp and scenario year
prepare_id <- function(x, target_col, group_var, start = 1, end = 100, order_y_axis,
                       h_levels = c(0.1, 0.5, 0.9), h_level_baseline){
  x$group_id <-  apply(x[group_var],MARGIN = 1, function(y) paste0(y, collapse = '_'))
  x['target_col'] <- x[target_col]
  
  n_obs <- x %>% 
    group_by(group_id) %>% 
    summarise(n = n(),
              step = end / n)
  #id_i <- x$group_id[1]
  df_out <-  purrr::map(unique(x$group_id), function(id_i){
    step <- n_obs %>% filter(group_id == id_i) %>% pull(step)
    n <- n_obs %>% filter(group_id == id_i) %>% pull(n)
    
    df_int <-       x %>% 
      filter(group_id == id_i)
    
    if(order_y_axis == 'descending'){
      df_int <- df_int %>% 
        group_by(group_id) %>% 
        arrange(desc(target_col),.by_group = TRUE) 
      
    } else {
      df_int <- df_int %>% 
        group_by(group_id) %>% 
        arrange((target_col),.by_group = TRUE)
      
    }  
    df_out <- df_int %>% 
      mutate(id = seq(from = start, to = end, length.out = n),
             step = step) %>% 
      ungroup() %>% 
      select(-group_id, -target_col)
    
    return(df_out)
  }) %>% 
    bind_rows() 
  
  if(is.null(h_levels)){
    quan_df_all <- NULL
  } else {
    #prepare contour df
    df_contour <- df_plot
    df_contour$group_id <-  apply(df_contour[group_var],MARGIN = 1, function(y) paste0(y, collapse = '_'))
    df_contour['target_col'] <- df_contour[target_col]
    
    #get quantiles of baseline
    baseline_id <- h_level_baseline
    
    df_baseline <- df_contour %>% 
      filter(grepl(pattern = baseline_id, x = group_id))
    
    quan_df_baseline <- purrr::map(h_levels, function(h) get_h_levels(df_baseline, h)) %>% 
      bind_rows() %>% 
      mutate(id_sub = gsub(pattern = baseline_id, replacement = '', x = group_id)) %>% 
      select(-group_id)
    
    #partly match each entry of quan_df_baseline
    quan_df_all <- purrr::map(1:nrow(quan_df_baseline), function(i){
      df_contour %>% 
        filter(grepl(pattern = quan_df_baseline$id_sub[i], x = group_id)) %>% 
        mutate(val = quan_df_baseline$q[i],
               quantile = quan_df_baseline$quantile[i],
               diff = abs(target_col - val)) %>% 
        arrange(diff) %>% 
        group_by(group_id) %>% 
        slice(1) %>% 
        ungroup() %>% 
        return()
      
    }) %>% 
      bind_rows()
  }
  
  return(list(df_plot = df_out, df_contour = quan_df_all))
}

library(colorRamps)
all_frost %>% 
  filter(location == 'CKA',
         Scenario != 'historical') %>% 
  mutate(scenario_type = factor(Scenario,
                                levels = c(1960, 1970, 1980, 1990, 2000, 2010, 2020, 2050, 2085),
                                labels = c(rep('Historical', 7), '2050', '2085')),
         scenario_x_label = ifelse(SSP == 'none', yes = Scenario, no = SSP)) %>% 
  #critical value habe ich von deinem code entnommen, markiert den Wendepunkt in der Farbskala (wenn ich es richtigverstehe)
  get_range_plot(critical_value = 20, 
                 #bestimmt den farbverlauf im plot
                 colours = rev(blue2yellow(15)), 
                 #welche spalte wird in dem plot dargestellt
                 target_col = 'frostH0',
                 #bestimmt die gruppen, die bei der Vorbereitung der Daten berücksichtig wird
                 group_var = c('scenario_type','scenario_x_label', 'SSP', 'cultivar'), 
                 #facet spalten
                 col_facet = 'scenario_type', 
                 #facet reihen
                 row_facet = 'cultivar',
                 #x-axis-ticks
                 x_col = 'scenario_x_label', 
                 #labeling of plot
                 x_label = 'Climate Scenario',
                 y_label = 'Cumulative Distribution of Frost Hours',
                 legend_name = 'Frost Hours',
                 #highest value at the bottom (can be reversed, usefull for plottomg bloom dates)
                 order_y_axis = 'descending', 
                 #which values get labelled isotopic lines 
                 d_levels = c(0, 10, 20, 30), 
                 #here you could specify the label the isotopic lines get (usefull for bloom dates)
                 #if null, label is the same as level
                 d_label = NULL,
                 #add vertical borders, to distinguish the colorstripes
                 add_border_color_stripes = TRUE,
                 #lower tile widths would make them thinner, so that coloms do not 'touch'
                 tile_width = 1,
                 #you can also draw intermediate, non-labelled isotopic lines
                 int_d_levels = c(15))
ggsave('plots/rangeplot_frostrisk_forstH0.jpeg',
       height = 20, width = 25, units = 'cm', device = 'jpeg')
