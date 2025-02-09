# generate temperature scenarios for Bonn
library(chillR)
library(tidyverse)
CKA_weather <- read_tab("data/cka_clean.csv")
fix_weather(CKA_weather)$QC
hist_scenarios <- temperature_scenario_from_records(weather = CKA_weather,
                                                   year = c(1960,1970,1980,1990,2000,2010,2020))

median(c(1958,2021))
baseline_scenario <- temperature_scenario_from_records(weather = CKA_weather,
                                                       year = 1989.5)

hist_adjusted <- temperature_scenario_baseline_adjustment(baseline_temperature_scenario = baseline_scenario,
                                                          temperature_scenario = hist_scenarios)

for (scen in hist_adjusted)
  {temps <- temperature_generation(weather = CKA_weather,
                                   years = c(1958, 2021),
                                   sim_years = c(3000, 3100),
                                   temperature_scenario = scen[[1]])
  write.csv(temps, paste0("weather/CKA_past_",scen$scenario_year,".csv"))
}

location <- c(7.1, 50.8)

area <- c(52, 6, 50, 8)

download_cmip6_ecmwfr(
  scenarios = c("ssp126", "ssp245", "ssp370", "ssp585"),
  area = area,
  user = '4cdedeb2-d3db-4c9e-8b7b-7ede68be142a',
  key = '48fe5079-b173-4339-af59-55abb2a936d4',
  model = 'all',
  frequency = 'monthly',
  variable = c('Tmin', 'Tmax'),
  year_start = 2015,
  year_end = 2100)


download_baseline_cmip6_ecmwfr(
  area = area,
  user = '4cdedeb2-d3db-4c9e-8b7b-7ede68be142a',
  key = '48fe5079-b173-4339-af59-55abb2a936d4',
  model = 'match_downloaded',
  frequency = 'monthly',
  variable = c('Tmin', 'Tmax'),
  year_start = 1986,
  year_end = 2014,
  month = 1:12)

station <- data.frame(
  station_name = c("Bonn"),
  longitude = c(7.1),
  latitude = c(50.8))

extracted <- extract_cmip6_data(stations = station)

change_scenarios <- 
  gen_rel_change_scenario(extracted,
                          scenarios = c(2050, 2085),
                          reference_period = c(1986:2014),
                          future_window_width = 30)

write.csv(change_scenarios, "data/all_change_scenarios.csv", row.names = FALSE)

change_scenarios <- read.csv("data/all_change_scenarios.csv")

scen_list <- convert_scen_information(change_scenarios,
                                      give_structure = FALSE)


baseline_scenario <- temperature_scenario_from_records(weather = CKA_weather,
                                                       year = 1989.5)
scen_reference_scenario <- temperature_scenario_from_records(weather = CKA_weather,
                                                             year = 2000)
baseline_adjustment_scenario <- temperature_scenario_baseline_adjustment(baseline_temperature_scenario = baseline_scenario,
                                                                         temperature_scenario = scen_reference_scenario)


adjusted_list <- temperature_scenario_baseline_adjustment(baseline_temperature_scenario = baseline_adjustment_scenario,
                                                          temperature_scenario = scen_list,
                                                          temperature_check_args = 
                                                            list(scenario_check_thresholds = c(-5, 15)))

for(scen in 1:length(adjusted_list))
{
  if(!file.exists(paste0("data/future_climate/CKA_future_",
                         scen,"_",
                         names(adjusted_list)[scen],".csv")) )
  {temp_temp <- temperature_generation(CKA_weather,
                                       years = c(1958, 2021),
                                       sim_years = c(3000, 3100),
                                       adjusted_list[scen],  
                                       temperature_check_args = 
                                         list( scenario_check_thresholds = c(-5, 20)))
  write.csv(temp_temp[[1]],paste0("data/future_climate/CKA_future_",scen,"_",names(adjusted_list)[scen],".csv"),
            row.names=FALSE)
  print(paste("Processed object",scen,"of", length(adjusted_list)))
  
  
  }
  
}



# now make projections
# par_all <- read.csv('data/parameter_cultivars.csv') %>% 
#   filter(species == 'Sweet Cherry') %>% 
#   mutate(Tc = 36, theta_star = 279)

par_best <- read.csv('data/par_cherry_best.csv')

#convert parameters from "new" format to old format (makes running PhenoFlex easier)
out <- par_best[,LarsChill::phenoflex_parnames_new] %>% 
  apply(MARGIN = 1, FUN = LarsChill::convert_parameters) %>% 
  t() 
out <- out[,5:8]
colnames(out) <- c('E0', 'E1', 'A0', 'A1')

#attach E0-A1 to the parameter table
par_best <- cbind(par_best, out)


#read the calibration / validation data
obs <- read.csv('data/master_phenology_repeated_splits.csv') %>% 
  filter(species == 'Sweet Cherry')

#read the weather data, make seasonlist object


cka_season <- read.csv('data/cka_clean.csv') %>% 
  chillR::stack_hourly_temps(latitude = 50.61) %>% 
  purrr::pluck('hourtemps') %>% 
  chillR::genSeasonList(years = 1959:2022) %>% 
  setNames(paste0('CKA--', 1959:2022))

seasonlist <- cka_season

#----------------------------#
#predict the bloom dates
#----------------------------#

#iterate over the rows of the parameter data.frame
#i <- 1
historical_predictions <- purrr::map(1:nrow(par_best), function(i){
  #extract the model parameters
  par <- par_best[i, LarsChill::phenoflex_parnames_old] %>% unlist()
  #iterate over model parameters
  
  #function iterates over the seasonlist and returns bloom dates
  pred <- LarsChill::return_predicted_days(par, modelfn = chillR::PhenoFlex_GDHwrapper, SeasonList = seasonlist, na_penalty = 365)
  
  loc <- names(pred) %>% str_split('--') %>% purrr::map_chr(1)
  yr <- names(pred) %>% str_split('--') %>% purrr::map_chr(2)
  
  return(data.frame(pred = pred, location = loc, year = yr, cultivar = par_all$cultivar[i], repetition = par_all$repetition[i]))
}, .progress = TRUE) %>% 
  bind_rows() %>%
  mutate(GCM="none",
         SSP="none",
         Time="historical")

write.csv(historical_predictions, "data/predictions/historical/cka_predictions_historical.csv",row.names = FALSE)

clim_scens <- list.files("data/future_climate/")

for(clim in clim_scens)
{
  seasonlist <- read.csv(paste0("data/future_climate/",clim)) %>% 
    chillR::stack_hourly_temps(latitude = 50.61) %>% 
    purrr::pluck('hourtemps') %>% 
    chillR::genSeasonList(years = 3001:3100) %>% 
    setNames(paste0('CKA--', 3001:3100))
  
  predictions <- purrr::map(1:nrow(par_best), function(i){
    #extract the model parameters
    par <- par_best[i, LarsChill::phenoflex_parnames_old] %>% unlist()
    #iterate over model parameters
    
    #function iterates over the seasonlist and returns bloom dates
    pred <- LarsChill::return_predicted_days(par, modelfn = chillR::PhenoFlex_GDHwrapper, SeasonList = seasonlist, na_penalty = 365)
    
    loc <- names(pred) %>% str_split('--') %>% purrr::map_chr(1)
    yr <- names(pred) %>% str_split('--') %>% purrr::map_chr(2)
    
    return(data.frame(pred = pred, location = loc, year = yr, cultivar = par_all$cultivar[i], repetition = par_all$repetition[i]))
  }, .progress = TRUE) %>% 
    bind_rows() %>%
    mutate(GCM=strsplit(clim,"\\.")[[1]][3],
           SSP=strsplit(clim,"\\.")[[1]][2],
           Time=strsplit(clim,"\\.")[[1]][4])
  write.csv(predictions, paste0("data/predictions/future_scenarios/cka_predictions_",clim),row.names = FALSE)
  
}  

clim_scens <- list.files("weather")

for(clim in clim_scens)
{
  seasonlist <- read.csv(paste0("weather/",clim)) %>% 
    chillR::stack_hourly_temps(latitude = 50.61) %>% 
    purrr::pluck('hourtemps') %>% 
    chillR::genSeasonList(years = 3001:3100) %>% 
    setNames(paste0('CKA--', 3001:3100))
  
  predictions <- purrr::map(1:nrow(par_best), function(i){
    #extract the model parameters
    par <- par_best[i, LarsChill::phenoflex_parnames_old] %>% unlist()
    #iterate over model parameters
    
    #function iterates over the seasonlist and returns bloom dates
    pred <- LarsChill::return_predicted_days(par, modelfn = chillR::PhenoFlex_GDHwrapper, SeasonList = seasonlist, na_penalty = 365)
    
    loc <- names(pred) %>% str_split('--') %>% purrr::map_chr(1)
    yr <- names(pred) %>% str_split('--') %>% purrr::map_chr(2)
    
    return(data.frame(pred = pred, location = loc, year = yr, cultivar = par_all$cultivar[i], repetition = par_all$repetition[i]))
  }, .progress = TRUE) %>% 
    bind_rows() %>%
    mutate(GCM="none",
           SSP="none",
           Time=strsplit(strsplit(clim,"_")[[1]][3],"\\.")[[1]][1])
  write.csv(predictions, paste0("data/predictions/past_scenarios/cka_predictions_",clim),row.names = FALSE)
  
}  


# characterize_frost

hist_preds <- list.files("data/predictions/historical")

hist_scen <- read.csv(paste0("data/predictions/historical/",hist_preds)) %>%
  mutate(bloom=round(pred))

weather <- read.csv("data/cka_clean.csv")

select_scen <- hist_scen %>%
  filter(cultivar %in% c("Bing", "Lapins", "Regina", "Schneiders", "Van"))

for(yy in unique(select_scen$year))  
{
 ww <-
   weather %>%
   mutate(Season=Year+(Month>7)) %>% 
   filter(Season==yy) %>% 
   make_JDay() %>%
   stack_hourly_temps(latitude = 50.61) %>% 
   purrr::pluck('hourtemps')
 
for (i in which(select_scen$year==yy))
  select_scen[i,10:273] <- ww[(which(ww$JDay==select_scen$bloom[i])[1]-120):(which(ww$JDay==select_scen$bloom[i])[24]+120),] %>%
   select(Temp) %>% unlist() %>% round(2)
}

write.csv(select_scen,"data/predictions/frost/historical.csv",row.names = FALSE)

# frost for past scenarios

hist_scens <- list.files("data/predictions/past_scenarios")

for(hs in hist_scens)
{
  hist_scen <- read.csv(paste0("data/predictions/past_scenarios/",hs)) %>%
    mutate(bloom=round(pred))
  select_scen <- hist_scen %>%
    filter(cultivar %in% c("Bing", "Lapins", "Regina", "Schneiders", "Van")) %>%
    filter(year %in% 3002:3100)
  
  weather <- read.csv(paste0("weather/CKA_past_", strsplit(hs,"_")[[1]][5]))
  
  for(yy in unique(select_scen$year))  
  {
    ww <-
      weather %>%
      mutate(Season=Year+(Month>7)) %>% 
      filter(Season==yy) %>% 
      make_JDay() %>%
      stack_hourly_temps(latitude = 50.61) %>% 
      purrr::pluck('hourtemps')
    
    for (i in which(select_scen$year==yy))
      select_scen[i,10:273] <- ww[(which(ww$JDay==select_scen$bloom[i])[1]-120):(which(ww$JDay==select_scen$bloom[i])[24]+120),] %>%
        select(Temp) %>% unlist() %>% round(2)
  }
  
  write.csv(select_scen,paste0("data/predictions/frost/pastscen_",strsplit(hs,"_")[[1]][5]),row.names = FALSE)
  
  
}

# frost for future scenarios

futu_scens <- list.files("data/predictions/future_scenarios")

for(fs in futu_scens)
{
  scen <- read.csv(paste0("data/predictions/future_scenarios/",fs)) %>%
    mutate(bloom=round(pred))
  select_scen <- scen %>%
    filter(cultivar %in% c("Bing", "Lapins", "Regina", "Schneiders", "Van")) %>%
    filter(year %in% 3002:3100)
  
  weather <- read.csv(paste0("data/future_climate/", strsplit(fs,"cka_predictions_")[[1]][2]))
  
  for(yy in unique(select_scen$year))  
  {
    ww <-
      weather %>%
      mutate(Season=Year+(Month>7)) %>% 
      filter(Season==yy) %>% 
      make_JDay() %>%
      stack_hourly_temps(latitude = 50.61) %>% 
      purrr::pluck('hourtemps')
    
    for (i in which(select_scen$year==yy))
      select_scen[i,10:273] <- ww[(which(ww$JDay==select_scen$bloom[i])[1]-120):(which(ww$JDay==select_scen$bloom[i])[24]+120),] %>%
        select(Temp) %>% unlist() %>% round(2)
  }
  
  write.csv(select_scen,paste0("data/predictions/frost/",strsplit(fs,"cka_predictions_CKA_")[[1]][2]),row.names = FALSE)
  
  
}


#### now evaluate the temperature windows (bloom date +- 5 days, at hourly basis)

frost_files <- list.files("data/predictions/frost")

for (ff in frost_files)
{
  fff <- read.csv(paste0("data/predictions/frost/",ff))
  meantemp <- apply(fff[,10:273],1,mean)
  mintemp <- apply(fff[,10:273],1,min)
  maxtemp <- apply(fff[,10:273],1,max)
  frost_hours_0 <- apply(fff[,10:273],1, function (x) sum(x<0))
  frost_hours_2.2 <- apply(fff[,10:273],1, function (x) sum(x<(-2.2))) 
  write.csv(data.frame(cbind(fff[,1:9],
                             Tmean = meantemp,
                             Tmin = mintemp,
                             Tmax = maxtemp,
                             frostH0 = frost_hours_0,
                             frostH2.2 = frost_hours_2.2)),
    paste0("data/predictions/frost_summary/",ff),row.names = FALSE)
}

frost_summary_files <- list.files("data/predictions/frost_summary")

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

ggplot(all_frost, aes(Scenario,Tmean)) + geom_boxplot() + facet_wrap(~cultivar)
ggplot(all_frost, aes(Scenario,bloomdate)) + geom_boxplot() + facet_wrap(~cultivar)

#historic plot
hist_data <- all_frost %>% filter(Scenario == "historical")

ggplot(hist_data, aes(Year,bloomdate)) + 
  geom_smooth() + 
  geom_point() +
  facet_wrap(~cultivar)

#past scenarios plot

past_data <- all_frost %>% filter(Scenario %in% c("1960", "1970", "1980", "1990",
                                  "2000", "2010", "2020"))
#%>%
#  mutate(Time = as.numeric(Time))

ggplot(past_data, aes(Scenario,bloomdate, group = Scenario)) + 
  geom_smooth() + 
  geom_boxplot() +
  facet_wrap(~cultivar)

#future scenarios plot
future_data <- all_frost %>% filter(Scenario %in% c("2050", "2085"))

ggplot(past_data, aes(Scenario,bloomdate, group = Scenario)) + 
  geom_boxplot() +
  facet_wrap(~cultivar)





plot_scenarios_gg <- function(past_observed,
                              past_simulated,
                              future_data,
                              metric,
                              axis_label,
                              time_points)
{
  rng <- range(past_observed[[metric]],
               past_simulated[[metric]],
               future_data[[metric]])  
  past_plot <- ggplot() +
    geom_boxplot(data = past_simulated %>% mutate(Scen = "Historical"),
                 aes_string("as.numeric(Scenario)",
                            metric,
                            group="Scenario"),
                 fill="skyblue") +
    scale_y_continuous(limits = c(0, 
                                  round(round(1.1*rng[2])))) +
    labs(x = "Year", y = axis_label) +
    facet_grid(~Scen) +
    theme_bw(base_size = 15) +  
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle=45, 
                                     hjust=1)) +
    geom_point(data = past_observed,
               aes_string("End_year",
                          metric),
               col="blue")
  
  future_plot_list <- list()
  
  for(y in time_points)
  {
    future_plot_list[[which(y == time_points)]] <-
      ggplot(data = future_data[which(future_data$Scenario==y),]) +
      geom_boxplot(aes_string("GCM", 
                              metric, 
                              fill="GCM")) +
      facet_wrap(vars(SSP), nrow = 1) +
      scale_x_discrete(labels = NULL,
                       expand = expansion(add = 1)) +
      scale_y_continuous(limits = c(0, 
                                    round(round(1.1*rng[2])))) +
      geom_text_npc(aes(npcx = "center",
                        npcy = "top",
                        label = Scenario),
                    size = 5) +
      theme_bw(base_size = 15) +
      theme(axis.ticks.y = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.margin = margin(0,
                                   0, 
                                   0, 
                                   0, 
                                   "cm"),
            legend.background = element_rect(),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            legend.box.spacing = unit(0, "cm"),
            plot.subtitle = element_text(hjust = 0.5,
                                         vjust = -1,
                                         size = 15 * 1.05,
                                         face = "bold")) 
  }
  
  plot <- (past_plot +
             future_plot_list +
             plot_layout(guides = "collect",
                         widths = c(1,rep(1.8,length(future_plot_list))))
  ) & theme(legend.position = "bottom",
            legend.key.spacing.y = unit(0,"mm"),
            legend.key.height = unit(4,"mm"),
            legend.text = element_text(size=8),
            legend.title = element_text(size=10),
            axis.title.x = element_blank(),
            legend.box.spacing = unit(0,"mm"))
  plot
  
}










plot_scenarios_gg_dates <- function(past_observed,
                              past_simulated,
                              future_data,
                              metric,
                              axis_label,
                              time_points)
{
  rng <- range(c(past_observed[[metric]],
               past_simulated[[metric]],
               future_data[[metric]]),na.rm=TRUE)  
  past_plot <- ggplot() +
    geom_boxplot(data = past_simulated %>% mutate(Scen = "Historical"),
                 aes_string("as.numeric(Scenario)",
                            metric,
                            group = "Scenario"),
                 fill = "skyblue",
                 na.rm = TRUE) +
    scale_y_date(limits=c(rng[1], rng[2]+(rng[2]-rng[1])*0.1)) +
    labs(x = "Year", y = axis_label) +
    facet_grid(~Scen) +
    theme_bw(base_size = 15) +  
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle=45, 
                                     hjust=1)) +
    geom_point(data = past_observed,
               aes_string("End_year",
                          metric),
               col="blue")
  
  future_plot_list <- list()
  
  for(y in time_points)
  {
    future_plot_list[[which(y == time_points)]] <-
      ggplot(data = future_data[which(future_data$Scenario==y),]) +
      geom_boxplot(aes_string("GCM", 
                              metric, 
                              fill = "GCM"),
                   na.rm = TRUE) +
      facet_wrap(vars(SSP), nrow = 1) +
      scale_x_discrete(labels = NULL,
                       expand = expansion(add = 1)) +
      scale_y_date(limits=c(rng[1], rng[2]+(rng[2]-rng[1])*0.1)) +
      geom_text_npc(aes(npcx = "center",
                        npcy = "top",
                        label = Scenario),
                    size = 5) +
      theme_bw(base_size = 15) +
      theme(axis.ticks.y = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.margin = margin(0,
                                   0, 
                                   0, 
                                   0, 
                                   "cm"),
            legend.background = element_rect(),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            legend.box.spacing = unit(0, "cm"),
            plot.subtitle = element_text(hjust = 0.5,
                                         vjust = -1,
                                         size = 15 * 1.05,
                                         face = "bold")) 
  }
  
  plot <- (past_plot +
             future_plot_list +
             plot_layout(guides = "collect",
                         widths = c(1,rep(1.8,length(future_plot_list))))
  ) & theme(legend.position = "bottom",
            legend.key.spacing.y = unit(0,"mm"),
            legend.key.height = unit(4,"mm"),
            legend.text = element_text(size=8),
            legend.title = element_text(size=10),
            axis.title.x = element_blank(),
            legend.box.spacing = unit(0,"mm"))
  plot
  
}


library(ggpp)
library(patchwork)

cultivars <- c("Bing","Van","Regina","Lapins", "Schneiders")

for (cultivar in cultivars)
  {plot_scenarios_gg_dates(past_observed = hist_data[hist_data$cultivar==cultivar,],
                  past_simulated = past_data[past_data$cultivar==cultivar,],
                  future_data = future_data[future_data$cultivar==cultivar,],
                  metric = "bloomdate",
                  axis_label = "Predicted bloom date",
                  time_points = c(2050, 2085))
  ggsave(paste0("plots/bloom_",cultivar,".png"),height=8,width=14,dpi=400)
}



library(colorRamps)

frost_sorting <- function(past_data, future_data, cultivar)
  {
  past_simulated <- past_data[past_data$cultivar==cultivar,]
  past_sorted <- past_simulated
  for(cult in unique(past_sorted$cultivar))
    for(Scenario in unique(past_sorted$Scenario))
      {past_sorted$frostH0[past_sorted$cultivar==cult &
                             past_sorted$Scenario==Scenario] <- 
        sort(past_sorted$frostH0[past_sorted$cultivar==cult &
                                   past_sorted$Scenario==as.character(Scenario)],
             decreasing = TRUE)
      past_sorted$frostH2.2[past_sorted$cultivar==cult &
                              past_sorted$Scenario==Scenario] <- 
        sort(past_sorted$frostH2.2[past_sorted$cultivar==cult &
                                     past_sorted$Scenario==as.character(Scenario)],
             decreasing = TRUE)
      past_sorted$Tmean[past_sorted$cultivar==cult &
                          past_sorted$Scenario==Scenario] <- 
        sort(past_sorted$Tmean[past_sorted$cultivar==cult &
                                 past_sorted$Scenario==as.character(Scenario)],
             decreasing = FALSE)
      past_sorted$Tmin[past_sorted$cultivar==cult &
                         past_sorted$Scenario==Scenario] <- 
        sort(past_sorted$Tmin[past_sorted$cultivar==cult &
                                past_sorted$Scenario==as.character(Scenario)],
             decreasing = FALSE)
      past_sorted$Tmax[past_sorted$cultivar==cult &
                         past_sorted$Scenario==Scenario] <- 
        sort(past_sorted$Tmax[past_sorted$cultivar==cult &
                                past_sorted$Scenario==as.character(Scenario)],
             decreasing = FALSE)
      
      }
  
  future_simulated <- future_data[future_data$cultivar==cultivar,]
  future_sorted <- future_simulated
  for(cult in unique(future_sorted$cultivar))
    for(Scenario in unique(future_sorted$Scenario))
      for(GCM in unique(future_sorted$GCM))
        for(SSP in unique(future_sorted$SSP))
          {future_sorted$frostH0[future_sorted$cultivar==cult &
                                   future_sorted$Scenario==Scenario &
                                   future_sorted$SSP==SSP &
                                   future_sorted$GCM==GCM] <- 
            sort(future_sorted$frostH0[future_sorted$cultivar==cult &
                                         future_sorted$Scenario==as.character(Scenario) &
                                         future_sorted$SSP==SSP &
                                         future_sorted$GCM==GCM],
                 decreasing = TRUE)
          future_sorted$frostH2.2[future_sorted$cultivar==cult &
                                    future_sorted$Scenario==Scenario &
                                    future_sorted$SSP==SSP &
                                    future_sorted$GCM==GCM] <- 
            sort(future_sorted$frostH2.2[future_sorted$cultivar==cult &
                                           future_sorted$Scenario==as.character(Scenario) &
                                           future_sorted$SSP==SSP &
                                           future_sorted$GCM==GCM],
                 decreasing = TRUE)
          future_sorted$Tmean[future_sorted$cultivar==cult &
                                future_sorted$Scenario==Scenario &
                                future_sorted$SSP==SSP &
                                future_sorted$GCM==GCM] <- 
            sort(future_sorted$Tmean[future_sorted$cultivar==cult &
                                       future_sorted$Scenario==as.character(Scenario) &
                                       future_sorted$SSP==SSP &
                                       future_sorted$GCM==GCM],
                 decreasing = FALSE)
          future_sorted$Tmin[future_sorted$cultivar==cult &
                               future_sorted$Scenario==Scenario &
                               future_sorted$SSP==SSP &
                               future_sorted$GCM==GCM] <- 
            sort(future_sorted$Tmin[future_sorted$cultivar==cult &
                                      future_sorted$Scenario==as.character(Scenario) &
                                      future_sorted$SSP==SSP &
                                      future_sorted$GCM==GCM],
                 decreasing = FALSE)
          future_sorted$Tmax[future_sorted$cultivar==cult &
                               future_sorted$Scenario==Scenario &
                               future_sorted$SSP==SSP &
                               future_sorted$GCM==GCM] <- 
            sort(future_sorted$Tmax[future_sorted$cultivar==cult &
                                      future_sorted$Scenario==as.character(Scenario) &
                                      future_sorted$SSP==SSP &
                                      future_sorted$GCM==GCM],
                 decreasing = FALSE)
        }
return(list(past=past_sorted,future=future_sorted))
  }




library(RColorBrewer)

time_points <- c(2050,2085)

tile_plots <- function(sorted_data, metric, critical_value, colours = rev(blue2yellow(15)))
{
  past_sorted <- sorted_data$past
  future_sorted <- sorted_data$future
  rng <- range(past_sorted[,metric],future_sorted[,metric])
  mx <- critical_value/rng[2]
  
  past_sorted$Scen <- "Historical"
  
  past_plot <- ggplot(past_sorted, aes(Scenario, Year, fill = .data[[metric]])) + 
    geom_tile() +
    facet_wrap(~Scen) +
    scale_fill_gradientn(colours = colours,
                         values = c(0,0.1,0.2,0.3,0.6,0.7,0.8,1,1/mx)*mx,"Frost hours",
                         limits = rng) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=1))
  
  future_plot_list <- list()
  
  for(y in time_points)
    {
    future_plot_list[[which(y == time_points)]] <-
      ggplot(data = future_sorted[which(future_sorted$Scenario == y),], aes(GCM, Year, fill = .data[[metric]])) +
      geom_tile() +
      facet_wrap(~SSP, nrow=1) +
      scale_fill_gradientn(colours = colours,
                           values = c(0,0.1,0.2,0.3,0.6,0.7,0.8,1,1/mx)*mx,"Frost hours",
                           limits = rng) +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0,
            size=6)) +
      ggtitle(y) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  plot <- (past_plot +
             future_plot_list +
             plot_layout(guides = "collect",
                         widths = c(1,rep(1.8,length(future_plot_list))))
  ) & theme(legend.position = "bottom",
            legend.key.spacing.y = unit(0,"mm"),
            legend.key.height = unit(4,"mm"),
            legend.text = element_text(size=8),
            legend.title = element_text(size=10),
            axis.title.x = element_blank(),
            legend.box.spacing = unit(0,"mm"))
  plot
  
  
}

for (cult in cultivars)
{
  sorted <- frost_sorting(past_data, future_data, cult)
  tile_plots(sorted, "frostH2.2", critical_value = 10, colours = ygobb(15))
  ggsave(paste0("plots/",cult,"_frostH2.2.png"),height=8,width=14,dpi=400)
  tile_plots(sorted, "frostH0", critical_value = 20, colours = ygobb(15))
  ggsave(paste0("plots/",cult,"_frostH0.png"),height=8,width=14,dpi=400)
  tile_plots(sorted, "Tmean", critical_value = 20, colours = matlab.like(15))
  ggsave(paste0("plots/",cult,"_Tmean.png"),height=8,width=14,dpi=400)
  tile_plots(sorted, "Tmin", critical_value = 10, colours = matlab.like(15))
  ggsave(paste0("plots/",cult,"_Tmin.png"),height=8,width=14,dpi=400)
  tile_plots(sorted, "Tmax", critical_value = 30, colours = matlab.like(15))
  ggsave(paste0("plots/",cult,"_Tmax.png"),height=8,width=14,dpi=400)
  
}



# compute some statistics

for (cult in cultivars)
  {
  cult_past <- past_data[which(past_data$cultivar==cult),]
  
  past <- cbind(Cultivar = cult,
                GCM = "none",
                SSP = "none",
                aggregate(cult_past["bloom"],
                          list(cult_past$Scenario), min),
                aggregate(cult_past["bloom"],
                          list(cult_past$Scenario), median)[,2],
                aggregate(cult_past["bloom"],
                          list(cult_past$Scenario), max)[,2],
                aggregate(cult_past[c("frostH0","frostH2.2","Tmean","Tmin","Tmax")],
                          list(cult_past$Scenario), mean)[,2:6],
                aggregate(cult_past[,c("frostH0","frostH2.2")],
                          list(cult_past$Scenario), 
                          function(x) sum(x >= 5))[,2:3],
                aggregate(cult_past[,c("frostH0","frostH2.2")],
                          list(cult_past$Scenario),
                          function(x) sum(x>0))[,2:3])
  
  colnames(past)[4:16] <- c("Time","earliest_bloom","median_bloom","latest_bloom",
                            "mean_frost0","mean_frost2.2",
                            "mean_Tmean","mean_Tmin","mean_Tmax",
                            "frost0_5plus","frost2.2_5plus",
                            "frost0_1plus","frost2.2_1plus")
  
  cult_future <- future_data[which(future_data$cultivar==cult),]
  
  future <- cbind(Cultivar = cult,
                  aggregate(cult_future["bloom"],
                            list(cult_future$GCM,cult_future$SSP,cult_future$Scenario), min),
                  aggregate(cult_future["bloom"],
                            list(cult_future$GCM,cult_future$SSP,cult_future$Scenario), median)[,4],
                  aggregate(cult_future["bloom"],
                            list(cult_future$GCM,cult_future$SSP,cult_future$Scenario), max)[,4],
                  aggregate(cult_future[c("frostH0","frostH2.2","Tmean","Tmin","Tmax")],
                            list(cult_future$GCM,cult_future$SSP,cult_future$Scenario), mean)[,4:8],
                  aggregate(cult_future[,c("frostH0","frostH2.2")],
                            list(cult_future$GCM,cult_future$SSP,cult_future$Scenario),
                            function(x) sum(x >= 5))[,4:5],
                  aggregate(cult_future[,c("frostH0","frostH2.2")],
                            list(cult_future$GCM,cult_future$SSP,cult_future$Scenario),
                            function(x) sum(x>0))[,4:5])
  
  colnames(future)<-c("Cultivar","GCM","SSP","Time","earliest_bloom","median_bloom","latest_bloom",
                      "mean_frost0","mean_frost2.2",
                      "mean_Tmean","mean_Tmin","mean_Tmax",
                      "frost0_5plus","frost2.2_5plus",
                      "frost0_1plus","frost2.2_1plus")                  
  
  cult_stats <- rbind(past,future)
  
  if(cult == cultivars[1]) all_stats <- cult_stats else
    all_stats <- rbind(all_stats,cult_stats)
  
  write.csv(all_stats, "data/all_frost_results.csv",row.names = FALSE)

}


summa <- data.frame(Metric = c("bloom_past_slope",
                               "bloom_1960_2020",
                               "bloom_2020_2050_SSP1",
                               "bloom_2020_2050_SSP2",
                               "bloom_2020_2050_SSP3",
                               "bloom_2020_2050_SSP5",
                               "bloom_2020_2085_SSP1",
                               "bloom_2020_2085_SSP2",
                               "bloom_2020_2085_SSP3",
                               "bloom_2020_2085_SSP5",
                               "frost2.2_5Hplus_1960",
                               "frost2.2_5Hplus_2020",
                               "frost2.2_5Hplus_2050_SSP1",
                               "frost2.2_5Hplus_2050_SSP2",
                               "frost2.2_5Hplus_2050_SSP3",
                               "frost2.2_5Hplus_2050_SSP5",
                               "frost2.2_5Hplus_2085_SSP1",
                               "frost2.2_5Hplus_2085_SSP2",
                               "frost2.2_5Hplus_2085_SSP3",
                               "frost2.2_5Hplus_2085_SSP5"))


for (cult in cultivars)
{
  bloom_medians <- all_stats %>% filter(Cultivar == cult) %>%
    group_by(Time,SSP) %>% 
    summarize(mean = mean(median_bloom),.groups = "drop_last")
  
  summa[summa$Metric == "bloom_past_slope",cult] <- 
    (lm(
      bloom_medians %>% filter(SSP=="none") %>% pluck("mean") ~
        as.numeric((bloom_medians %>% filter(SSP=="none") %>% pluck("Time"))))$coefficients[2]*10) %>% round(2)

  summa[summa$Metric == "bloom_1960_2020",cult] <- bloom_medians %>% filter(Time == 2020) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 1960) %>% pluck("mean")
  
  summa[summa$Metric == "bloom_2020_2050_SSP1",cult] <- (bloom_medians %>% filter(SSP == "SSP1", Time == 2050) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 2020) %>% pluck("mean")) %>% round(2)
  summa[summa$Metric == "bloom_2020_2050_SSP2",cult] <- (bloom_medians %>% filter(SSP == "SSP2", Time == 2050) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 2020) %>% pluck("mean")) %>% round(2)
  summa[summa$Metric == "bloom_2020_2050_SSP3",cult] <- (bloom_medians %>% filter(SSP == "SSP3", Time == 2050) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 2020) %>% pluck("mean")) %>% round(2)
  summa[summa$Metric == "bloom_2020_2050_SSP5",cult] <- (bloom_medians %>% filter(SSP == "SSP5", Time == 2050) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 2020) %>% pluck("mean")) %>% round(2)

  summa[summa$Metric == "bloom_2020_2085_SSP1",cult] <- (bloom_medians %>% filter(SSP == "SSP1", Time == 2085) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 2020) %>% pluck("mean")) %>% round(2)
  summa[summa$Metric == "bloom_2020_2085_SSP2",cult] <- (bloom_medians %>% filter(SSP == "SSP2", Time == 2085) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 2020) %>% pluck("mean")) %>% round(2)
  summa[summa$Metric == "bloom_2020_2085_SSP3",cult] <- (bloom_medians %>% filter(SSP == "SSP3", Time == 2085) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 2020) %>% pluck("mean")) %>% round(2)
  summa[summa$Metric == "bloom_2020_2085_SSP5",cult] <- (bloom_medians %>% filter(SSP == "SSP5", Time == 2085) %>% pluck("mean") -
    bloom_medians %>% filter(Time == 2020) %>% pluck("mean")) %>% round(2)
  
  frost2.2_5Hplus_medians <- all_stats %>% filter(Cultivar == cult) %>%
    group_by(Time,SSP) %>% 
    summarize(mean = mean(frost2.2_5plus),.groups = "drop_last")
  
  summa[summa$Metric == "frost2.2_5Hplus_1960",cult] <- frost2.2_5Hplus_medians %>% filter(Time == 1960) %>% pluck("mean") %>% round(2)
  summa[summa$Metric == "frost2.2_5Hplus_2020",cult] <- frost2.2_5Hplus_medians %>% filter(Time == 2020) %>% pluck("mean") %>% round(2)

  summa[summa$Metric == "frost2.2_5Hplus_2050_SSP1",cult] <- frost2.2_5Hplus_medians %>% filter(SSP == "SSP1", Time == 2050) %>% pluck("mean") %>% round(2)
  summa[summa$Metric == "frost2.2_5Hplus_2050_SSP2",cult] <- frost2.2_5Hplus_medians %>% filter(SSP == "SSP2", Time == 2050) %>% pluck("mean") %>% round(2)
  summa[summa$Metric == "frost2.2_5Hplus_2050_SSP3",cult] <- frost2.2_5Hplus_medians %>% filter(SSP == "SSP3", Time == 2050) %>% pluck("mean") %>% round(2)
  summa[summa$Metric == "frost2.2_5Hplus_2050_SSP5",cult] <- frost2.2_5Hplus_medians %>% filter(SSP == "SSP5", Time == 2050) %>% pluck("mean") %>% round(2)
  
  summa[summa$Metric == "frost2.2_5Hplus_2085_SSP1",cult] <- frost2.2_5Hplus_medians %>% filter(SSP == "SSP1", Time == 2085) %>% pluck("mean") %>% round(2)
  summa[summa$Metric == "frost2.2_5Hplus_2085_SSP2",cult] <- frost2.2_5Hplus_medians %>% filter(SSP == "SSP2", Time == 2085) %>% pluck("mean") %>% round(2)
  summa[summa$Metric == "frost2.2_5Hplus_2085_SSP3",cult] <- frost2.2_5Hplus_medians %>% filter(SSP == "SSP3", Time == 2085) %>% pluck("mean") %>% round(2)
  summa[summa$Metric == "frost2.2_5Hplus_2085_SSP5",cult] <- frost2.2_5Hplus_medians %>% filter(SSP == "SSP5", Time == 2085) %>% pluck("mean") %>% round(2)
  
}

write.csv(summa, "data/all_frost_results_summary.csv", row.names = FALSE)


