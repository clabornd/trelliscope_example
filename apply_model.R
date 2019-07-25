source('sourced_script/model_fourier.R')
source('sourced_script/helper_functions.R')
source('sourced_script/plot_fn.R')

library(trelliscopejs)
library(plotly)
library(lme4)
library(tidyverse)
library(lubridate)
library(chron)

model_df <- read.csv('data/plotting_df.csv', check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(t_awake = ifelse(DayNight == 'day', 
                          yes = ifelse(hours <= 24 & hours >= 6, hours - 6, 18 + hours),
                          no = ifelse(hours <= 24 & hours >= 18, hours - 18, 6 + hours))) %>%
  mutate(cost = cos(2*pi*hours/24), sint = sin(2*pi*hours/24),
         cost12 = cos(2*pi*hours/12), sint12 = sin(2*pi*hours/12),
         cost8 = cos(2*pi*hours/8), sint8 = sin(2*pi*hours/8))

cats = c("1:30", "4:30",  "7:30", "10:30", "13:30", "16:30", "22:30")

model_df <- model_df %>% group_by(Lipid) %>% nest()

model_df_params <- model_df %>%
  mutate(params = map(data, model_fs))

# try plotting one
df <- model_df_params$data[[1]]
params <- model_df_params$params[[1]]
lipidplot_nfits_facet(df, params, subjects_day = TRUE, subjects_night = FALSE, fit_day_24 = TRUE, fit_night_24 = FALSE)

facet_plots <- model_df_params %>%
  mutate(plot = pmap_plot(list(df = data, params = params, subjects_day = TRUE, subjects_night = FALSE, fit_day_24 = TRUE, fit_night_24 = FALSE,
                               cats = list(cats), period = list(c(24,12,8))), 
                          lipidplot_nfits_facet)) 

facet_plots_cogs <- facet_plots %>% 
  mutate(Lipid = cog(Lipid, desc = "Lipid Identifier", default_label = T)) %>%
  cbind(do.call(rbind.data.frame, facet_plots$params) %>% mutate_each(cog))

facet_plots_cogs %>% trelliscope(name = 'Plots of all lipids', nrow = 1, ncol = 1, path = 'displays/', thumb = TRUE, width = 1000, desc = 'Plots of all lipids')
