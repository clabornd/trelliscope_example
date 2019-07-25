# Model function which returns:  
# 1) Point estimates for parameters and day/night differences of those parameters.  Standard errors and p-values for each estimate
# 2) p-value from test of equality of means between vectors of amplitude estimates for Day/Night groups for individual
#' @param df dataframe for a single lipid containing the columns "hours", "SubjectID", "Value", and "DayNight", "cost", and "sint"
model_fs = function(df){
  
  mixedmod1 <- lmer(Value ~ DayNight + DayNight:(cost24 + sint24) + DayNight:(cost12 + sint12) + DayNight:(cost8 + sint8) + (1|SubjectID) - 1 , data = df)
  
  # Fit Parameters
  summary_obj = summary(mixedmod1)

  var_cov <- vcov(mixedmod1) # in order: M_day, M_night, beta_day, beta_night, gamma_day, gamma_night
  pt_est <- summary_obj$coefficients[1:14] 
  stderrs <- summary_obj$coefficients[15:28]
  
  ### Mesor estimates and SE's for day, night and Day/Night difference
  mesor_day <- pt_est[1]
  mesor_night <- pt_est[2]
  mesor_diff <- mesor_night-mesor_day
  SE_mesor_day <- stderrs[1]
  SE_mesor_night <- stderrs[2]
  SE_mesordiff <- sqrt(SE_mesor_day^2+SE_mesor_night^2)
  t_mesordiff <- mesor_diff/SE_mesordiff
  p_mesordiff <- pt(t_mesordiff, 12, lower.tail = t_mesordiff < 0)*2
  #
  
  # Amplitude
  amp_day24 <- sqrt(pt_est[3]^2+pt_est[5]^2)
  amp_day12 <- sqrt(pt_est[7]^2+pt_est[9]^2)
  amp_day8 <- sqrt(pt_est[11]^2+pt_est[13]^2)
  
  amp_night24 <- sqrt(pt_est[4]^2+pt_est[6]^2)
  amp_night12 <- sqrt(pt_est[8]^2+pt_est[10]^2)
  amp_night8 <- sqrt(pt_est[12]^2+pt_est[14]^2)
  
  # Phi
  phi_day24 <- atan2(-pt_est[5], pt_est[3])
  phi_day12 <- atan2(-pt_est[9], pt_est[7])
  phi_day8 <- atan2(-pt_est[13], pt_est[11])
  
  phi_night24 <- atan2(-pt_est[6], pt_est[4])
  phi_night12 <- atan2(-pt_est[10], pt_est[8])
  phi_night8 <- atan2(-pt_est[14], pt_est[12])
  
  # Phi
  phi_day24 <- atan2(-pt_est[5], pt_est[3])
  phi_day12 <- atan2(-pt_est[9], pt_est[7])
  phi_day8 <- atan2(-pt_est[13], pt_est[11])
  
  phi_day24_HM <- to_hours(phi_day24)
  phi_day12_HM <- to_hours(phi_day12)
  phi_day8_HM <- to_hours(phi_day8)
  
  phi_night24 <- atan2(-pt_est[6], pt_est[4])
  phi_night12 <- atan2(-pt_est[10], pt_est[8])
  phi_night8 <- atan2(-pt_est[14], pt_est[12])
  
  phi_night24_HM <- to_hours(phi_night24)
  phi_night12_HM <- to_hours(phi_night12)
  phi_night8_HM <- to_hours(phi_night8)
  
  # create 2x2 variance covariance matrix for the amplitude_day and phi_day
  sigma_day24 <- matrix(c(var_cov[3,3], var_cov[5,3], var_cov[3,5], var_cov[5,5]), nrow = 2)
  sigma_day12 <- matrix(c(var_cov[7,7], var_cov[9,7], var_cov[7,9], var_cov[9,9]), nrow = 2)
  sigma_day8 <- matrix(c(var_cov[11,11], var_cov[13,11], var_cov[11,13], var_cov[13,13]), nrow = 2)
  
  sigma_night24 <- matrix(c(var_cov[4,4], var_cov[4,6], var_cov[6,4], var_cov[6,6]), nrow = 2)
  sigma_night12 <- matrix(c(var_cov[8,8], var_cov[10,8], var_cov[8,10], var_cov[10,10]), nrow = 2)
  sigma_night8 <- matrix(c(var_cov[12,12], var_cov[14,12], var_cov[12,14], var_cov[14,14]), nrow = 2)
  
  # pass point estimates and var-covariance matrix to delta method functions #
  
  # day
  SE_amp_day24 <- (deltamethod_amp(pt_est[3], pt_est[5], sigma_day24)^0.5)[1]
  SE_phi_day24 <- (deltamethod_phi(pt_est[3], pt_est[5], sigma_day24)^0.5)[1]
  SE_phi_day24_HM <- to_hours(SE_phi_day24, FALSE)
  
  SE_amp_day12 <- (deltamethod_amp(pt_est[7], pt_est[9], sigma_day12)^0.5)[1]
  SE_phi_day12 <- (deltamethod_phi(pt_est[7], pt_est[9], sigma_day12)^0.5)[1]
  SE_phi_day12_HM <- to_hours(SE_phi_day12, FALSE)
  
  SE_amp_day8 <- (deltamethod_amp(pt_est[11], pt_est[13], sigma_day8)^0.5)[1]
  SE_phi_day8 <- (deltamethod_phi(pt_est[11], pt_est[13], sigma_day8)^0.5)[1]
  SE_phi_day8_HM <- to_hours(SE_phi_day8, FALSE)
  
  # night
  SE_amp_night24 <- (deltamethod_amp(pt_est[3], pt_est[5], sigma_night24)^0.5)[1]
  SE_phi_night24 <- (deltamethod_phi(pt_est[3], pt_est[5], sigma_night24)^0.5)[1]
  SE_phi_night24_HM <- to_hours(SE_phi_night24, FALSE)
  
  SE_amp_night12 <- (deltamethod_amp(pt_est[7], pt_est[9], sigma_night12)^0.5)[1]
  SE_phi_night12 <- (deltamethod_phi(pt_est[7], pt_est[9], sigma_night12)^0.5)[1]
  SE_phi_night12_HM <- to_hours(SE_phi_night12, FALSE)
  
  SE_amp_night8 <- (deltamethod_amp(pt_est[11], pt_est[13], sigma_night8)^0.5)[1]
  SE_phi_night8 <- (deltamethod_phi(pt_est[11], pt_est[13], sigma_night8)^0.5)[1]
  SE_phi_night8_HM <- to_hours(SE_phi_night8, FALSE)
  
  ## t statistics and p values for amplitude:
  
  #...day
  t_amp_day24 <- amp_day24/SE_amp_day24
  p_amp_day24 <- pt(t_amp_day24, 6, lower.tail = FALSE)
  
  t_amp_day12 <- amp_day12/SE_amp_day12
  p_amp_day12 <- pt(t_amp_day12, 6, lower.tail = FALSE)
  
  t_amp_day8 <- amp_day8/SE_amp_day8
  p_amp_day8 <- pt(t_amp_day8, 6, lower.tail = FALSE)
  
  #...night
  t_amp_night24 <- amp_night24/SE_amp_night24
  p_amp_night24 <- pt(t_amp_night24, 6, lower.tail = FALSE)
  
  t_amp_night12 <- amp_night12/SE_amp_night12
  p_amp_night12 <- pt(t_amp_night12, 6, lower.tail = FALSE)
  
  t_amp_night8 <- amp_night8/SE_amp_night8
  p_amp_night8 <- pt(t_amp_night8, 6, lower.tail = FALSE)
  ##
  
  # point estimates for differences in amplitude/phase between day and night groups
  amp_diff24 <- amp_day24 - amp_night24
  amp_diff12 <- amp_day12 - amp_night12
  amp_diff8 <- amp_day8 - amp_night8
  
  phi_diff24 <- pi - abs(abs(phi_night24 - phi_day24) - pi)
  phi_diff24_HM <- to_hours(phi_diff24, signed = FALSE)
  phi_diff12 <- pi - abs(abs(phi_night12 - phi_day12) - pi)
  phi_diff12_HM <- to_hours(phi_diff12, signed = FALSE)
  phi_diff8 <- pi - abs(abs(phi_night8 - phi_day8) - pi)
  phi_diff8_HM <- to_hours(phi_diff8, signed = FALSE)
  
  # standard error for differences between day and night group
  SE_ampdiff24 <- sqrt(SE_amp_day24^2+SE_amp_night24^2)
  SE_phidiff24 <- sqrt(SE_phi_day24^2+SE_phi_night24^2)
  SE_ampdiff12 <- sqrt(SE_amp_day12^2+SE_amp_night12^2)
  SE_phidiff12 <- sqrt(SE_phi_day12^2+SE_phi_night12^2)
  SE_ampdiff8 <- sqrt(SE_amp_day8^2+SE_amp_night8^2)
  SE_phidiff8 <- sqrt(SE_phi_day8^2+SE_phi_night8^2)
  
  ## t-statistic and p-value for differences
  
  # amplitude
  t_amp_diff24 <- amp_diff24/SE_ampdiff24
  p_amp_diff24 <- 2*pt(t_amp_diff24, 12, lower.tail = t_amp_diff24 < 0)
  t_amp_diff12 <- amp_diff12/SE_ampdiff12
  p_amp_diff12 <- 2*pt(t_amp_diff12, 12, lower.tail = t_amp_diff12 < 0)
  t_amp_diff8 <- amp_diff8/SE_ampdiff8
  p_amp_diff8 <- 2*pt(t_amp_diff8, 12, lower.tail = t_amp_diff8 < 0)
  
  # phi
  t_phi_diff24 <- phi_diff24/SE_phidiff24
  p_phi_diff24 <- 2*pt(t_phi_diff24, 12, lower.tail = t_phi_diff24 < 0)
  t_phi_diff12 <- phi_diff12/SE_phidiff12
  p_phi_diff12 <- 2*pt(t_phi_diff12, 12, lower.tail = t_phi_diff12 < 0)
  t_phi_diff8 <- phi_diff8/SE_phidiff8
  p_phi_diff8 <- 2*pt(t_phi_diff8, 12, lower.tail = t_phi_diff8 < 0)
  
  power_1_day <- 8*amp_day24^2/4
  power_1_night <- 8*amp_night24^2/4
  power_2_day <- 8*amp_day12^2/4
  power_2_night <- 8*amp_night12^2/4
  power_3_day <- 8*amp_day8^2/4
  power_3_night <- 8*amp_night8^2/4
  
  P_tot <- power_1_day + power_1_night + power_2_day + power_2_night + power_3_day + power_3_night
  P_tot_day <- power_1_day + power_2_day + power_3_day 
  P_tot_night <- power_1_night + power_2_night + power_3_night
  log_SNR <- 10*log10(P_tot/sigma(mixedmod1)^2)
  log_SNR_day <- 10*log10(P_tot_day/sigma(mixedmod1)^2)
  log_SNR_night <- 10*log10(P_tot_night/sigma(mixedmod1)^2)
  
  # store values in a list
  params <- list(mesor_day = mesor_day,
                 mesor_night = mesor_night,
                 mesor_diff = mesor_diff,
                 
                 log_SNR = log_SNR,
                 log_SNR_day = log_SNR_day,
                 log_SNR_night = log_SNR_night,
                 
                 # 24 HOUR ESTIMATES
                 amp_day24 = amp_day24,
                 amp_night24 = amp_night24, 
                 amp_diff24 = amp_diff24,
                 phi_day24 = phi_day24, 
                 phi_day24_HM = phi_day24_HM, 
                 phi_night24 = phi_night24, 
                 phi_night24_HM = phi_night24_HM,
                 phi_diff24 = phi_diff24, 
                 phi_diff24_HM = phi_diff24_HM, 
                 
                 SE_mesor_day = SE_mesor_day, 
                 SE_mesor_night = SE_mesor_night, 
                 SE_mesordiff = SE_mesordiff, 
                 SE_amp_day24 = SE_amp_day24, 
                 SE_amp_night24 = SE_amp_night24, 
                 SE_ampdiff24 = SE_ampdiff24, 
                 SE_phi_day24 = SE_phi_day24, 
                 SE_phi_day24_HM = SE_phi_day24_HM,
                 SE_phi_night24 = SE_phi_night24, 
                 SE_phi_night24_HM = SE_phi_night24_HM,
                 SE_phidiff24 = SE_phidiff24,
                 
                 p_mesordiff = p_mesordiff,   
                 p_amp_day24 = p_amp_day24,
                 p_amp_night24 = p_amp_night24,
                 p_amp_diff24 = p_amp_diff24,
                 p_phi_diff24 = p_phi_diff24,
                 
                 # 12 HOUR ESTIMATES
                 amp_day12 = amp_day12, 
                 amp_night12 = amp_night12, 
                 amp_diff12 = amp_diff12,
                 phi_day12 = phi_day12, 
                 phi_day12_HM = phi_day12_HM, 
                 phi_night12 = phi_night12, 
                 phi_night12_HM = phi_night12_HM,
                 phi_diff12 = phi_diff12, 
                 phi_diff12_HM = phi_diff12_HM, 
                 
                 SE_amp_day12 = SE_amp_day12, 
                 SE_amp_night12 = SE_amp_night12, 
                 SE_ampdiff12 = SE_ampdiff12, 
                 SE_phi_day12 = SE_phi_day12, 
                 SE_phi_day12_HM = SE_phi_day12_HM,
                 SE_phi_night12 = SE_phi_night12, 
                 SE_phi_night12_HM = SE_phi_night12_HM,
                 SE_phidiff12 = SE_phidiff12,
                 
                 p_amp_day12 = p_amp_day12,
                 p_amp_night12 = p_amp_night12,
                 p_amp_diff12 = p_amp_diff12,
                 p_phi_diff12 = p_phi_diff12,
                 
                 # 8 hour estimates
                 amp_day8 = amp_day8, 
                 amp_night8 = amp_night8, 
                 amp_diff8 = amp_diff8,
                 phi_day8 = phi_day8, 
                 phi_day8_HM = phi_day8_HM, 
                 phi_night8 = phi_night8, 
                 phi_night8_HM = phi_night8_HM,
                 phi_diff8 = phi_diff8, 
                 phi_diff8_HM = phi_diff8_HM, 
                 
                 SE_amp_day8 = SE_amp_day8, 
                 SE_amp_night8 = SE_amp_night8, 
                 SE_ampdiff8 = SE_ampdiff8, 
                 SE_phi_day8 = SE_phi_day8, 
                 SE_phi_day8_HM = SE_phi_day8_HM,
                 SE_phi_night8 = SE_phi_night8, 
                 SE_phi_night8_HM = SE_phi_night8_HM,
                 SE_phidiff8 = SE_phidiff8,
                 
                 p_amp_day8 = p_amp_day8,
                 p_amp_night8 = p_amp_night8,
                 p_amp_diff8 = p_amp_diff8,
                 p_phi_diff8 = p_phi_diff8
  )
  
  return(params)
  
}