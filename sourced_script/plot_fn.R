# functions to plot fitted curves
# just a single sinusoid
cos_fit <- function(x, M, A, phi, period = 24){M + A*cos(2*pi*x/period + phi)}

# multiple sinusoids
cos_fitn <- function(x, M, amp, phi, period){
  sum = M
  for(i in 1:(length(amp))){
    sum <- sum + amp[i]*cos(2*pi*x/period[i] + phi[i])
  }
  return(sum)
}

#' Plot a single panel one period fit 
#' @param df dataframe for a singe lipid - same dataframe passed to model()
#' @param params list output from model()
lipidplot = function(df, params, subjects_day = TRUE, subjects_night = TRUE, fit_day = TRUE, fit_night = TRUE, period = 24,
                            cats = c("1:30", "4:30",  "7:30","10:30", "13:30", "16:30", "19:30", "22:30")){

  # subset by individual, sort by numeric time
  df <- df %>% group_by(SubjectID) %>% 
    arrange(hours) %>% 
    ungroup()
  
  p <- plot_ly(type = 'scatter', mode = 'lines+markers')
  
  # create dataframes that contain fitted values against time
  fit_df_day <- data.frame(x = seq(0, 24, 0.05)) %>% 
    mutate(fit = cos_fit(x, params[["mesor_day"]], params[["amp_day"]], params[["phi_day"]], period = period))
  fit_df_night <- data.frame(x = seq(0, 24, 0.05)) %>% 
    mutate(fit = cos_fit(x, params[["mesor_night"]], params[["amp_night"]], params[["phi_night"]], period = period))
  
  # add traces separately to be able to hide them if necessary
  p <- p %>%
    add_trace(data = df %>% filter(DayNight == "day"), x = ~hours, y = ~Value,
              name = ~SubjectID, colors = "Set2", line = list(dash = 'dot'), marker = list(color = "red"), visible = subjects_day) %>%
    add_trace(data = df %>% filter(DayNight == "night"), x = ~hours, y = ~Value,
              name = ~SubjectID, colors = "Set1", marker = list(color = "blue"), visible = subjects_night) %>%
    add_trace(data = fit_df_day, x = ~x, y = ~fit, line = list(color = 'red', width = 6, dash = 'dot'), mode = 'lines', name = "Day fit", visible = fit_day) %>%
    add_trace(data = fit_df_night, x = ~x, y = ~fit, line = list(color = 'blue', width = 6), mode = 'lines', name = "Night fit", visible = fit_night)
  
  # styling
  p %>% layout(yaxis = list(title = "Normalized Log2 Abundance"), xaxis = list(title = "Collection Time", tickvals = unique(df$hours), ticktext = cats))
}

# single panel of multi-period model (up to three periods)
lipidplot_nfits = function(df, params, subjects_day = TRUE, subjects_night = TRUE, 
                                    fit_day_24 = TRUE, fit_night_24 = TRUE,
                                    fit_day_both = 'legendonly', fit_night_both = 'legendonly',
                                    period = c(24,12,8),
                                    cats = c("1:30", "4:30",  "7:30","10:30", "13:30", "16:30", "19:30", "22:30")){

  # subset by individual, sort by numeric time
  df <- df %>% group_by(SubjectID) %>% 
    arrange(hours) %>% 
    ungroup()
  
  p <- plot_ly(type = 'scatter', mode = 'lines+markers')
  
  # fits from the full model (there is no trace for the model fit with only a 24 hour period)
  
  # just the 24 hour component
  fit_df_day <- data.frame(x = seq(0, 24, 0.05)) %>% 
    mutate(fit = cos_fit(x, params[["mesor_day"]], params[["amp_day24"]], params[["phi_day24"]]))
  fit_df_night <- data.frame(x = seq(0, 24, 0.05)) %>% 
    mutate(fit = cos_fit(x, params[["mesor_night"]], params[["amp_night24"]], params[["phi_night24"]]))
  
  # both the 24 and 12 hour component
  fit_df_day2 <- data.frame(x = seq(0, 24, 0.05)) %>%
    mutate(fit = cos_fitn(x, params[["mesor_day"]], 
                          amp = c(params[["amp_day24"]], params[["amp_day12"]], params[['amp_day8']]), 
                          phi = c(params[["phi_day24"]], params[["phi_day12"]], params[['phi_day8']]),
                          period = period)
    )
           
  fit_df_night2 <- data.frame(x = seq(0, 24, 0.05)) %>%
    mutate(fit = cos_fitn(x, params[["mesor_night"]], 
                          amp = c(params[["amp_night24"]], params[["amp_night12"]], params[['amp_night8']]), 
                          phi = c(params[["phi_night24"]], params[["phi_night12"]], params[['phi_night8']]),
                          period = period)
    )

  # separate traces to control visibility
  p <- p %>%
    add_trace(data = df %>% filter(DayNight == "day"), x = ~hours, y = ~Value,
              name = ~SubjectID, colors = "Set2", line = list(dash = 'dot'), marker = list(color = "red"), visible = subjects_day) %>%
    add_trace(data = df %>% filter(DayNight == "night"), x = ~hours, y = ~Value,
              name = ~SubjectID, colors = "Set1", marker = list(color = "blue"), visible = subjects_night) %>%
    add_trace(data = fit_df_day, x = ~x, y = ~fit, line = list(color = 'red', width = 6, dash = 'dot'), mode = 'lines', name = "24 hour component(day)", visible = fit_day_24) %>%
    add_trace(data = fit_df_night, x = ~x, y = ~fit, line = list(color = 'blue', width = 6), mode = 'lines', name = "24 hour component(night)", visible = fit_night_24) %>%
    add_trace(data = fit_df_day2, x = ~x, y = ~fit, line = list(color = 'red', width = 6, dash = 'dot'), mode = 'lines', name = "Full fit(day)", visible = fit_day_both) %>%
    add_trace(data = fit_df_night2, x = ~x, y = ~fit, line = list(color = 'blue', width = 6), mode = 'lines', name = "Full fit(night)", visible = fit_night_both)
  p %>% layout(yaxis = list(title = "Normalized Log2 Abundance"), xaxis = list(title = "Collection Time", tickvals = unique(df$hours), ticktext = cats))
}

## functions for faceted displays:

# faceted single fit
lipidplot_facet <- function(df, params, subjects_day = TRUE, subjects_night = TRUE, fit_day = TRUE, fit_night = TRUE, period = 24,
                            cats = c("1:30", "4:30",  "7:30","10:30", "13:30", "16:30", "19:30", "22:30")){
  p1 <- lipidplot(df, params, subjects_day, subjects_night, fit_day, fit_night, period = period, cats = cats)
  p2 <- lipidplot(df, params, !subjects_day, !subjects_night, !fit_day, !fit_night, period = period, cats = cats)
  subplot(p1, p2, nrows = 1)
}


# faceted two period fit
lipidplot_nfits_facet <- function(df, params, subjects_day = TRUE, subjects_night = TRUE, 
                                    fit_day_24 = TRUE, fit_night_24 = TRUE,
                                    fit_day_both = 'legendonly', fit_night_both = 'legendonly',
                                    period = c(24,12,8),
                                    cats = c("1:30", "4:30",  "7:30","10:30", "13:30", "16:30", "19:30", "22:30")){
  p1 <- lipidplot_nfits(df, params, subjects_day, subjects_night, fit_day_24, fit_night_24, fit_day_both, fit_night_both,
                                period = period, cats = cats)
  p2 <- lipidplot_nfits(df, params, !subjects_day, !subjects_night, !fit_day_24, !fit_night_24, fit_day_both, fit_night_both,
                                 period = period, cats = cats)
  subplot(p1, p2, nrows = 1)
}


