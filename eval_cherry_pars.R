library(tidyverse)

devtools::install_github('https://github.com/larscaspersen/eval_phenoflex')
devtools::install_github('https://github.com/larscaspersen/addition_chillR')
library(LarsChill)
library(evalpheno)

par_all <- read.csv('data/parameter_cultivars.csv') %>% 
  filter(species == 'Sweet Cherry') %>% 
  mutate(Tc = 36, theta_star = 279)

#convert parameters from "new" format to old format (makes running PhenoFlex easier)
out <- par_all[,LarsChill::phenoflex_parnames_new] %>% 
  apply(MARGIN = 1, FUN = LarsChill::convert_parameters) %>% 
  t() 
out <- out[,5:8]
colnames(out) <- c('E0', 'E1', 'A0', 'A1')

#attach E0-A1 to the parameter table
par_all <- cbind(par_all, out)


#read the calibration / validation data
obs <- read.csv('data/master_phenology_repeated_splits.csv') %>% 
  filter(species == 'Sweet Cherry')

#read the weather data, make seasonlist object
zaragoza_season <- read.csv('data/zaragoza_clean.csv') %>% 
  chillR::stack_hourly_temps(latitude = 41.65) %>% 
  purrr::pluck('hourtemps') %>% 
  chillR::genSeasonList(years = 1974:2022) %>% 
  setNames(paste0('Zaragoza--', 1974:2022))

cka_season <- read.csv('data/cka_clean.csv') %>% 
  chillR::stack_hourly_temps(latitude = 50.61) %>% 
  purrr::pluck('hourtemps') %>% 
  chillR::genSeasonList(years = 1958:2022) %>% 
  setNames(paste0('Klein-Altendorf--', 1958:2022))

seasonlist <- c(cka_season, zaragoza_season)
rm(zaragoza_season, cka_season)

#----------------------------#
#predict the bloom dates
#----------------------------#

#iterate over the rows of the parameter data.frame
#i <- 1
pred_df <- purrr::map(1:nrow(par_all), function(i){
  #extract the model parameters
  par <- par_all[i, LarsChill::phenoflex_parnames_old] %>% unlist()
  #iterate over model parameters
  
  #function iterates over the seasonlist and returns bloom dates
  pred <- LarsChill::return_predicted_days(par, modelfn = chillR::PhenoFlex_GDHwrapper, SeasonList = seasonlist, na_penalty = 365)
  
  loc <- names(pred) %>% str_split('--') %>% purrr::map_chr(1)
  yr <- names(pred) %>% str_split('--') %>% purrr::map_chr(2)
  
  return(data.frame(pred = pred, location = loc, year = yr, cultivar = par_all$cultivar[i], repetition = par_all$repetition[i]))
}, .progress = TRUE) %>% 
  bind_rows()

#merge with calibration / validation data
str(pred_df)

pred_obs <- pred_df %>% 
  mutate(year = as.numeric(year)) %>% 
  merge(obs, by = c('cultivar', 'repetition', 'location', 'year'))


#calculate the rmse, rpiq (for rpiq we use interquartile range of the data before splitting)
iqr_df <- obs %>% 
  group_by(repetition, cultivar) %>% 
  summarise(n_total = n(),
            iqr = IQR(pheno))

performance <- pred_obs %>% 
  group_by(cultivar, repetition, split) %>% 
  summarise(mean_bias = mean(pred - pheno) %>%  round(digits = 2),
            rmse = chillR::RMSEP(pred, pheno) %>% round(digits = 2)) %>% 
  merge(iqr_df, by = c('repetition', 'cultivar')) %>% 
  mutate(rpiq = round(iqr / rmse, digits = 2))

#take the best performing parameters (validation?) per cultivar
min_rmse <- performance %>% 
  filter(split == 'Validation') %>% 
  group_by(cultivar,) %>% 
  summarise(rmse = min(rmse))


min_rmse_df <- performance %>% 
  filter(split == 'Validation') %>% merge(min_rmse, by = c('cultivar', 'rmse')) %>% 
  mutate(r_cult = paste(repetition, cultivar))



#make a plot of predicted vs observed for each cultivar, add info on calibration / validation performance
performance_sub <- performance %>% 
  mutate(r_cult = paste(repetition, cultivar)) %>% 
  filter(r_cult %in% min_rmse_df$r_cult) %>% 
  pivot_longer(cols = c('rmse', 'rpiq', 'mean_bias')) %>% 
  mutate(measurement_split = paste(name, split, sep = '_'),
         value = round(value, digits = 1)) %>% 
  select(-name, -split) %>% 
  pivot_wider(names_from = measurement_split, values_from = value)

length(unique(performance_sub$cultivar))
#make 3 3x3 plots and one larger one

ypos_text <- 120
xpos_text <- 65
cult_select <- sort(unique(performance_sub$cultivar))[1:9]
pred_obs %>% 
  mutate(r_cult = paste(repetition, cultivar)) %>% 
  filter(r_cult %in% min_rmse_df$r_cult,
         cultivar %in% cult_select) %>% 
  ggplot(aes(x = pheno, y = pred)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_text(data = performance_sub[performance_sub$cultivar %in% cult_select,],
            aes(x = xpos_text, y = ypos_text, label = paste0('Calibration (Validation)\nRMSE: ', 
                                                format(rmse_Calibration, digits = 2), 
                                                ' (', 
                                                format(rmse_Validation, digits = 2),
                                                ')\nRPIQ: ', 
                                                format(rpiq_Calibration, digits = 2),
                                                ' (',
                                                format(rpiq_Validation, digits = 2),
                                                ')\nMean Bias: ',
                                                format(mean_bias_Calibration, digits = 1),
                                                ' (',
                                                format(mean_bias_Validation, digits = 1),')')),
            ,hjust = 0) +
  geom_point(aes(col = split)) +
  facet_wrap(~cultivar) +
  theme_bw(base_size = 15) +
  ylab('Predicted Bloom Date') +
  xlab('Observed Bloom Date') +
  scale_x_continuous(breaks = c(1, 32, 60,91, 121, 152), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'), 
                     limits = c(60, 140)) +
  scale_y_continuous(breaks = c(1, 32, 60,91, 121, 152), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'),
                     limits = c(60, 140)) +
  theme(legend.position = 'bottom') 
ggsave('plots/pred_obs_cherry_1.jpeg', height = 20, width =25, units = 'cm', device = 'jpeg')

cult_select <- sort(unique(performance_sub$cultivar))[10:18]
pred_obs %>% 
  mutate(r_cult = paste(repetition, cultivar)) %>% 
  filter(r_cult %in% min_rmse_df$r_cult,
         cultivar %in% cult_select) %>% 
  ggplot(aes(x = pheno, y = pred)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_text(data = performance_sub[performance_sub$cultivar %in% cult_select,],
            aes(x = xpos_text, y = ypos_text, label = paste0('Calibration (Validation)\nRMSE: ', 
                                                format(rmse_Calibration, digits = 2), 
                                                ' (', 
                                                format(rmse_Validation, digits = 2),
                                                ')\nRPIQ: ', 
                                                format(rpiq_Calibration, digits = 2),
                                                ' (',
                                                format(rpiq_Validation, digits = 2),
                                                ')\nMean Bias: ',
                                                format(mean_bias_Calibration, digits = 1),
                                                ' (',
                                                format(mean_bias_Validation, digits = 1),')')),
            ,hjust = 0) +
  geom_point(aes(col = split)) +
  facet_wrap(~cultivar) +
  theme_bw(base_size = 15) +
  ylab('Predicted Bloom Date') +
  xlab('Observed Bloom Date') +
  scale_x_continuous(breaks = c(1, 32, 60,91, 121, 152), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'), 
                     limits = c(60, 140)) +
  scale_y_continuous(breaks = c(1, 32, 60,91, 121, 152), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'),
                     limits = c(60, 140)) +
  theme(legend.position = 'bottom') 
ggsave('plots/pred_obs_cherry_2.jpeg', height = 20, width =25, units = 'cm', device = 'jpeg')

cult_select <- sort(unique(performance_sub$cultivar))[19:27]
pred_obs %>% 
  mutate(r_cult = paste(repetition, cultivar)) %>% 
  filter(r_cult %in% min_rmse_df$r_cult,
         cultivar %in% cult_select) %>% 
  ggplot(aes(x = pheno, y = pred)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_text(data = performance_sub[performance_sub$cultivar %in% cult_select,],
            aes(x = xpos_text, y = ypos_text, label = paste0('Calibration (Validation)\nRMSE: ', 
                                                format(rmse_Calibration, digits = 2), 
                                                ' (', 
                                                format(rmse_Validation, digits = 2),
                                                ')\nRPIQ: ', 
                                                format(rpiq_Calibration, digits = 2),
                                                ' (',
                                                format(rpiq_Validation, digits = 2),
                                                ')\nMean Bias: ',
                                                format(mean_bias_Calibration, digits = 1),
                                                ' (',
                                                format(mean_bias_Validation, digits = 1),')')),
            ,hjust = 0) +
  geom_point(aes(col = split)) +
  facet_wrap(~cultivar) +
  theme_bw(base_size = 15) +
  ylab('Predicted Bloom Date') +
  xlab('Observed Bloom Date') +
  scale_x_continuous(breaks = c(1, 32, 60,91, 121, 152), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'), 
                     limits = c(60, 140)) +
  scale_y_continuous(breaks = c(1, 32, 60,91, 121, 152), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'),
                     limits = c(60, 140)) +
  theme(legend.position = 'bottom') 
ggsave('plots/pred_obs_cherry_3.jpeg', height = 20, width =25, units = 'cm', device = 'jpeg')

cult_select <- sort(unique(performance_sub$cultivar))[28:37]
pred_obs %>% 
  mutate(r_cult = paste(repetition, cultivar)) %>% 
  filter(r_cult %in% min_rmse_df$r_cult,
         cultivar %in% cult_select) %>% 
  ggplot(aes(x = pheno, y = pred)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_text(data = performance_sub[performance_sub$cultivar %in% cult_select,],
            aes(x = xpos_text, y = ypos_text, label = paste0('Calibration (Validation)\nRMSE: ', 
                                                format(rmse_Calibration, digits = 2), 
                                                ' (', 
                                                format(rmse_Validation, digits = 2),
                                                ')\nRPIQ: ', 
                                                format(rpiq_Calibration, digits = 2),
                                                ' (',
                                                format(rpiq_Validation, digits = 2),
                                                ')\nMean Bias: ',
                                                format(mean_bias_Calibration, digits = 1),
                                                ' (',
                                                format(mean_bias_Validation, digits = 1),')')),
            ,hjust = 0) +
  geom_point(aes(col = split)) +
  facet_wrap(~cultivar) +
  theme_bw(base_size = 15) +
  ylab('Predicted Bloom Date') +
  xlab('Observed Bloom Date') +
  scale_x_continuous(breaks = c(1, 32, 60,91, 121, 152), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'), 
                     limits = c(60, 140)) +
  scale_y_continuous(breaks = c(1, 32, 60,91, 121, 152), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'),
                     limits = c(60, 140)) +
  theme(legend.position = 'bottom') 
ggsave('plots/pred_obs_cherry_4.jpeg', height = 20, width =25, units = 'cm', device = 'jpeg')







