#######
library(dplyr)
library(tidyr)
library(readr)

n_mols <- 100
n_subjects <- 20
subject_effects <- runif(n_subjects, -1, 1)

amps_g1_24 <- runif(n_mols, 1, 10)
amps_g1_12 <- runif(n_mols, 1, 5)
amps_g1_8 <- runif(n_mols, 1, 2.5)
amps_g2_24 <- runif(n_mols, 1, 10)
amps_g2_12 <- runif(n_mols, 1, 5)
amps_g2_8 <- runif(n_mols, 1, 2.5)
phs_g1_24 <- runif(n_mols, -pi, pi)
phs_g1_12 <- runif(n_mols, -pi, pi)
phs_g1_8 <- runif(n_mols, -pi, pi)
phs_g2_24 <- runif(n_mols, -pi, pi)
phs_g2_12 <- runif(n_mols, -pi, pi)
phs_g2_8 <- runif(n_mols, -pi, pi)
mesors_g1 <- runif(n_mols, -4, 4)
mesors_g2 <- runif(n_mols, -4, 4)

times = c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5, 22.5)

# a function which generates data based on a lipid index, and a subject effect
make_subject_data <- function(i, subject_effect, group = 'day'){
  if(group == 'day'){
    mesors_g1[i] + amps_g1_24[i]*cos(2*pi*times/24 + phs_g1_24[i]) + amps_g1_12[i]*cos(2*pi*times/12 + phs_g1_12[i]) + amps_g1_8[i]*cos(2*pi*times/8 + phs_g1_8[i]) + rnorm(length(times), rnorm(1, subject_effect, 0.1), 0.1)
  }
  else{
    mesors_g2[i] + amps_g2_24[i]*cos(2*pi*times/24 + phs_g2_24[i]) + amps_g2_12[i]*cos(2*pi*times/12 + phs_g2_12[i]) + amps_g2_8[i]*cos(2*pi*times/8 + phs_g2_8[i]) + rnorm(length(times), rnorm(1, subject_effect, 0.1), 0.1)
  }
}

df <- data.frame()

for(i in 1:n_mols){
  row <- NULL
  for(j in 1:n_subjects){
    grp <- if(j > 10) 'night' else 'day'
    row <- c(row, make_subject_data(i, subject_effects[j], group = grp))
  }  
  df <- df %>% rbind(row)
}

subject_names <- paste0('Subject_', 1:n_subjects)
colnames <- lapply(subject_names, function(name){
  paste0(name, '_time_', 1:length(times))
}) %>% unlist()

colnames(df) <- colnames
lipids <- readRDS('data/lipids.RDS')
df <- df %>% mutate(Lipid = sample(lipids, 100))
plot_df <- df %>% gather(Sample, Value, -'Lipid')

groupdata <- data.frame(Sample = colnames, DayNight = rep(c('day', 'night'), each = 80)) %>%
  mutate(hours = rep(times, 20), SubjectID = rep(subject_names, each = 8))

plot_df <- plot_df %>% left_join(groupdata)

write_csv(plot_df, 'data/plotting_df.csv')
