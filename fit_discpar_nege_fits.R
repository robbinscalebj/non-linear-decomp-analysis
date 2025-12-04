
litterfitter_pkg<-devtools::as.package("D:/Github_Repos/litterfitter") 
devtools::install(litterfitter_pkg); library(litterfitter) #need to build/install the package into the library to work with foreach
library(tidyverse)
library(broom)
library(foreach); library(doSNOW); library(parallel)
library(here)

# Load and Tidy Data
litter_df <- read_csv("https://github.com/robbinscalebj/revisiting-k/raw/refs/heads/main/data/derived_data/tidied_detrital_nutrients_data.csv")|> # tidied litter time series
  mutate(mass_prop = Mass_per_remaining/100)|>
  group_by(cohort_id)|>
  mutate(series_length = n())|>
  ungroup()|>
  filter(series_length > 5)|>
  rename(mass = "mass_prop", time = "Meas_Day")

nested_data <- litter_df|>
  select(cohort_id, time, mass)|>
  group_by(cohort_id)|> # needs group_by()|>nest() structure... nest_by() doesn't work!
  nest()|>
  ungroup()|>
  mutate(batch_number = row_number())#|>
  #slice_head(n = 4) # for testing

# Definitions for use in parallel loops
n_data_batches <- nrow(nested_data)
cohort_names <- nested_data|>pull(cohort_id)
n_iters <- 9999

pb <- txtProgressBar(max = n_data_batches, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Set parallel processing
cl <- makePSOCKcluster(16) 
registerDoSNOW(cl)


# Loop

nege_fits <- foreach(i = 1:nrow(nested_data), .combine = 'rbind', .options.snow = opts, 
                       .packages = c("tidyverse", "litterfitter", "future"), .inorder = TRUE#, .final = function(x) setNames(x, cohort_names)
) %dopar% {
  
  # filter to cohort data and unnest into normal data frame
  df <- nested_data|>filter(batch_number == i)|>unnest()
  
  cohort_id <- df|>first()|>pull(cohort_id)
  
  nege_fit <- fit_litter(time = df$time,  
                           mass.remaining = df$mass,
                           upper = c(0.65, 20,35,5000), # k with upper at 2 and lower at 0.1 lots of boundary fits
                           lower = c(0,0.005,0.1,10), #r, k, alpha, beta
                           model = "discpar.negexp.genexp",
                           iters = 9999) # need far fewer iterations for only one parameter
  
 
  
  # Package up results for export
  tibble(cohort_id = cohort_id,
         nege_fit = list(nege_fit)) 
  
  
}     # ~9 minutes on 16 cores; all converged

stopCluster(cl)


# save to file
saveRDS(nege_fits, file = here("DetNut_nege_fits.rds")) 
