#'@param beta point estimate for beta from the cosinor or mixed model output
#'@param gamma point estimate for gamma from the cosinor or mixed model output
#'@param sigma variance covariance matrix for the two parameters with the variance for beta and gamma at cells [1,1] and [2,2] and cov(beta, gamma) on the off diagonals

deltamethod_amp = function(beta, gamma, sigma){
  grad1 <- beta/sqrt(beta^2 + gamma^2)
  grad2 <- gamma/sqrt(beta^2 + gamma^2)
  
  grad <- c(grad1, grad2)
  
  return(grad%*%sigma%*%grad)
}

deltamethod_phi = function(beta, gamma, sigma){
  grad1 <- gamma/(beta^2*(1+(-gamma/beta)^2))
  grad2 <- -1/(beta*(1+(-gamma/beta)^2))
  
  grad <- c(grad1, grad2)
  
  return(grad%*%sigma%*%grad)
}

# convert radians to HH:MM string
to_hm <- function(phi, signed = TRUE){
  if(signed){ paste0(gsub(":[0-9]{2}$", "", chron::times(ifelse(sign(phi) <= 0, abs(phi)/2/pi, 1 - abs(phi)/2/pi)))) }
  else paste0(gsub(":[0-9]{2}$", "", chron::times(abs(phi)/2/pi)))
}

# convert radians to decimal hours
to_hours <- function(phi, signed = TRUE){
  if(signed){ifelse(sign(phi) <= 0, abs(phi)%%(2*pi)/2/pi*24, (2*pi - phi%%(2*pi))%%(2*pi)/2/pi*24)}
  else(abs(phi)/2/pi*24)
}

# make cog columns from params column
make_cogs <- function(df, full_model = FALSE){
  adj_p_amp_day = p.adjust(sapply(df$params, function(x) x$p_amp_day), method = 'BH')
  adj_p_amp_night = p.adjust(sapply(df$params, function(x) x$p_amp_night), method = 'BH')
  adj_p_amp_diff = p.adjust(sapply(df$params, function(x) x$p_amp_diff), method = 'BH')
  adj_p_phi_diff = p.adjust(sapply(df$params, function(x) x$p_phi_diff), method = 'BH')
  
  df <- df %>% 
    mutate(
      Lipid = cog(Lipid, desc = "Lipid Identifier"),
      Ionization = cog(Ionization, desc = "Instrument Mode", default_active = FALSE),
      correlations = map_cog(correlations, ~as_tibble(.)),
      additional_cogs = map_cog(params,
            ~tibble(
              icc = cog(.x$icc, desc = 'Intra-Class Correlation Coefficient for single period model'),
              icc_day = cog(.x$icc_day, desc = 'Intra- Class Correlation Coefficient for the day group in the single period model'),
              icc_night = cog(.x$icc_night, desc = 'Intra- Class Correlation Coefficient for the night group in the single period model'),
              # AIC = cog(.x$AIC, desc = 'AIC for single period model'),
              deviance = cog(.x$deviance, desc = 'Deviance for single period model'),
              # mean_cv_day = cog(.x$mean_cv_day, desc = 'Average coefficient of variation across all time points (day group)'),
              # mean_cv_night = cog(.x$mean_cv_night, desc = 'Average coefficient of variation across all time points (night group)'),
              mesor_day = cog(.x$mesor_day, desc = "Midline of lipid abundance for day shift condition", group = "Point Estimates"),
              mesor_night = cog(.x$mesor_night, desc = "Midline of lipid abundance for night shift condition", group = "Point Estimates"),
              # mesor_diff = cog(.x$mesor_diff, desc = "Difference in midline (night vs day shift)", group = "Point Estimates"),
              amp_day = cog(.x$amp_day, desc = "Amplitude for day shift condition", group = "Point Estimates"), 
              amp_night = cog(.x$amp_night, desc = "Amplitude for night shift condition", group = "Point Estimates"), 
              # amp_diff = cog(.x$amp_diff, desc = "Difference in amplitudes (night vs day shift)", group = "Point Estimates"),
              phi_day = cog(.x$phi_day, desc = "Time of maximum abundance in radians for day shift condition", group = "Point Estimates"), 
              phi_day_HM = cog(.x$phi_day_HM, desc = "Time of maximum abundance in decimal hours for day shift condition", group = "Point Estimates"), 
              phi_night = cog(.x$phi_night, desc = "Time of maximum abundance in radians for night shift condition", group = "Point Estimates"), 
              phi_night_HM = cog(.x$phi_night_HM, desc = "Time of maximum abundance in decimal hours for night shift condition", group = "Point Estimates"),
              phi_diff = cog(.x$phi_diff, desc = "Difference in time of maximum abundance in radians", group = "Point Estimates"), 
              phi_diff_HM = cog(.x$phi_diff_HM, desc = "Difference in time of maximum abundance in decimal hours", group = "Point Estimates"), 
              
              #SE_mesor_day = cog(.x$SE_mesor_day, desc = "Standard error for midline estimate (day shift)", group = "Standard Errors"), 
              #SE_mesor_night = cog(.x$SE_mesor_night, desc = "Standard error for midline estimate (night shift)", group = "Standard Errors"), 
              #SE_mesordiff = cog(.x$SE_mesordiff, desc = "Standard error for the difference in midline estimates", group = "Standard Errors"), 
              #SE_amp_day = cog(.x$SE_amp_day, desc = "Standard error of amplitude estimate (day shift)", group = "Standard Errors"), 
              #SE_amp_night = cog(.x$SE_amp_night, desc = "Standard error of amplitude estimate (night shift)", group = "Standard Errors"), 
              #SE_ampdiff = cog(.x$SE_ampdiff, desc = "Standard error for difference in amplitude esimates", group = "Standard Errors"), 
              #SE_phi_day = cog(.x$SE_phi_day, desc = "Standard error for peak time estimate in radians (day shift)", group = "Standard Errors"), 
              #SE_phi_day_HM = cog(.x$SE_phi_day_HM, desc = "Standard error for peak time estimate in decimal hours (day shift)", group = "Standard Errors"),
              #SE_phi_night = cog(.x$SE_phi_night, desc = "Standard error for peak time estimate in radians (night shift)", group = "Standard Errors"), 
              #SE_phi_night_HM = cog(.x$SE_phi_night_HM, desc = "Standard error for peak time estimate in decimal hours (night shift)", group = "Standard Errors"),
              #SE_phidiff = cog(.x$SE_phidiff, desc = "Standard error for difference in peak time estimates in radians", group = "Standard Errors"),
              
              # p_mesordiff = cog(.x$p_mesordiff, desc = "p value for test of no difference between midlines (night vs day shifts)", group = "p values"),   
              p_amp_day = cog(.x$p_amp_day, desc = "p value for test of zero amplitude (day shift)", group = "p values"),
              p_amp_night = cog(.x$p_amp_night, desc = "p value for test of zero amplitude (night shift)", group = "p values"),
              p_amp_diff = cog(.x$p_amp_diff, desc = "p value for test of no difference in amplitudes (night vs day shifts)", group = "p values"),
              p_phi_diff = cog(.x$p_phi_diff, desc = "p value for test of no difference in peak time (night vs day shifts)", group = "p values"),
              
              # t_mesordiff = cog(t_mesordiff, desc = "t statistic for difference in midline estimates", group = "t values"),
              # t_amp_day = cog(t_amp_day, desc = "t statistic for amplitude estimate (day shift)", group = "t values"), 
              # t_amp_night = cog(t_amp_night, desc = "t statistic for amplitude estimate (night shift)", group = "t values"), 
              # t_amp_diff = cog(t_amp_diff, desc = "t statistic for difference in amplitude estimates", group = "t values"), 
              # t_phi_diff = cog(t_phi_diff, desc = "t statistic for difference in phi estimates", group = "t values"),
            )
        ),
      # adj_p_mesordiff = cog(p.adjust(.x$p_mesordiff, method = "BH"), desc = "Benjamini-Hochberg adjusted p value for test of no difference between midlines (night vs day shifts)", group = "p values"),
      adj_p_amp_day = cog(adj_p_amp_day, desc = "Benjamini-Hochberg adjusted p value for test of zero amplitude (day shift)", group = "p values"),
      adj_p_amp_night = cog(adj_p_amp_night, desc = "Benjamini-Hochberg adjusted p value for test of zero amplitude (night shift)", group = "p values"),
      adj_p_amp_diff = cog(adj_p_amp_diff, desc = "Benjamini-Hochberg adjusted p value for test of no difference in amplitudes (night vs day shifts)", group = "p values"),
      adj_p_phi_diff = cog(adj_p_phi_diff, desc = "Benjamini-Hochberg adjusted p value for test of no difference in peak time (night vs day shifts)", group = "p values")
    )
  
  if(full_model){
    adj_p_amp_day24 = p.adjust(sapply(df$params, function(x) x$p_amp_day24), method = 'BH')
    adj_p_amp_night24 = p.adjust(sapply(df$params, function(x) x$p_amp_night24), method = 'BH')
    adj_p_amp_diff24 = p.adjust(sapply(df$params, function(x) x$p_amp_diff24), method = 'BH')
    adj_p_phi_diff24 = p.adjust(sapply(df$params, function(x) x$p_phi_diff24), method = 'BH')
    
    adj_p_amp_day12 = p.adjust(sapply(df$params, function(x) x$p_amp_day12), method = 'BH')
    adj_p_amp_night12 = p.adjust(sapply(df$params, function(x) x$p_amp_night12), method = 'BH')
    adj_p_amp_diff12 = p.adjust(sapply(df$params, function(x) x$p_amp_diff12), method = 'BH')
    adj_p_phi_diff12 = p.adjust(sapply(df$params, function(x) x$p_phi_diff12), method = 'BH')
    
    df <- df %>% 
      mutate(
        full_cogs = map_cog(params,
          ~tibble(
            icc_full = cog(.x$icc_full, desc = 'Intra- Class Correlation Coefficient for full model'),
            icc_day_full = cog(.x$icc_day_full, desc = 'Intra- Class Correlation Coefficient for the day group in the full model'),
            icc_night_full = cog(.x$icc_night_full, desc = 'Intra- Class Correlation Coefficient for the night group in the full model'),
            # AIC_full = cog(.x$AIC_full, desc = 'AIC for two period model'),
            deviance_full = cog(.x$deviance_full, desc = 'Deviance for two period model'),
            p_lratio = cog(.x$p_likelihood_ratio, desc = 'p-value for likelihood ratio test'),
            
            mesor_day_full = cog(.x$mesor_day_full, 'Full model day mesor'),
            mesor_night_full = cog(.x$mesor_night_full, 'Full model night mesor'),
            
            amp_day24 = cog(.x$amp_day24,'Day amplitude estimate for the 24 hour period in full model'), 
            amp_night24 = cog(.x$amp_night24, 'Night amplitude estimate for the 24 hour period in full model'),
            phi_day24 = cog(.x$phi_day24, 'Day phi estimate for the 24 hour period in full model'),
            phi_day_HM24 = cog(.x$phi_day_HM24, 'Day phi estimate (in decimal hours) for the 24 hour period in full model'),
            phi_night24 = cog(.x$phi_night24, 'Night phi estimate for the 24 hour period in full model'), 
            phi_night_HM24 = cog(.x$phi_night_HM24, 'Night phi estimate (in decimal hours) for the 24 hour period in full model'),
            
            amp_day12 = cog(.x$amp_day12,'Day amplitude estimate for the 12 hour period in full model'), 
            amp_night12 = cog(.x$amp_night12, 'Night amplitude estimate for the 12 hour period in full model'),
            phi_day12 = cog(.x$phi_day12, 'Day phi estimate for the 12 hour period in full model'),
            phi_day_HM12 = cog(.x$phi_day_HM12, 'Day phi estimate (in decimal hours) for the 12 hour period in full model'),
            phi_night12 = cog(.x$phi_night12, 'Night phi estimate for the 12 hour period in full model'), 
            phi_night_HM12 = cog(.x$phi_night_HM12, 'Night phi estimate (in decimal hours) for the 12 hour period in full model'),
            
            p_amp_day24 = cog(.x$p_amp_day24, desc = "p value for test of zero amplitude (day shift, 24hr)", group = "p values"),
            p_amp_night24 = cog(.x$p_amp_night24, desc = "p value for test of zero amplitude (night shift, 24hr)", group = "p values"),
            p_amp_diff24 = cog(.x$p_amp_diff24, desc = "p value for test of no difference in amplitudes (night vs day shifts, 24hr)", group = "p values"),
            p_phi_diff24 = cog(.x$p_phi_diff24, desc = "p value for test of no difference in peak time (night vs day shifts, 24hr)", group = "p values"),
            
            p_amp_day12 = cog(.x$p_amp_day12 , desc = "p value for test of zero amplitude (day shift, 12hr)", group = "p values"),
            p_amp_night12 = cog(.x$p_amp_night12, desc = "p value for test of zero amplitude (night shift, 12hr)", group = "p values"),
            p_amp_diff12 = cog(.x$p_amp_diff12, desc = "p value for test of no difference in amplitudes (night vs day shifts, 12hr)", group = "p values"),
            p_phi_diff12 = cog(.x$p_phi_diff12, desc = "p value for test of no difference in peak time (night vs day shifts, 12hr)", group = "p values")
            )
          ),
        adj_p_amp_day24 = cog(adj_p_amp_day24, desc = "Benjamini-Hochberg adjusted p value for test of zero amplitude (day shift, 24hr)", group = "p values"),
        adj_p_amp_night24 = cog(adj_p_amp_night24, desc = "Benjamini-Hochberg adjusted p value for test of zero amplitude (night shift, 24hr)", group = "p values"),
        adj_p_amp_diff24 = cog(adj_p_amp_diff24, desc = "Benjamini-Hochberg adjusted p value for test of no difference in amplitudes (night vs day shifts, 24hr)", group = "p values"),
        adj_p_phi_diff24 = cog(adj_p_phi_diff24, desc = "Benjamini-Hochberg adjusted p value for test of no difference in peak time (night vs day shifts, 24hr)", group = "p values"),
        
        adj_p_amp_day12 = cog(adj_p_amp_day12 , desc = "Benjamini-Hochberg adjusted p value for test of zero amplitude (day shift, 12hr)", group = "p values"),
        adj_p_amp_night12 = cog(adj_p_amp_night12, desc = "Benjamini-Hochberg adjusted p value for test of zero amplitude (night shift, 12hr)", group = "p values"),
        adj_p_amp_diff12 = cog(adj_p_amp_diff12, desc = "Benjamini-Hochberg adjusted p value for test of no difference in amplitudes (night vs day shifts, 12hr)", group = "p values"),
        adj_p_phi_diff12 = cog(adj_p_phi_diff12, desc = "Benjamini-Hochberg adjusted p value for test of no difference in peak time (night vs day shifts, 12hr)", group = "p values")
        ) 
  }
  
  return(df)
}