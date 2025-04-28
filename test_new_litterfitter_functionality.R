# Test ground for adding Hessian calculation to litterfitter package fork
devtools::load_all("D:/Github_Repos/litterfitter") #loads local litterfitter without installing into library
library(tidyverse)
library(here)
library(devtools)



# load test time series

test_series <- read_csv("https://github.com/robbinscalebj/revisiting-k/raw/refs/heads/main/data/derived_data/tidied_detrital_nutrients_data.csv")|> # tidied litter time series
  mutate(mass_prop = Mass_per_remaining/100)|>
  group_by(cohort_id)|>
  mutate(series_length = n())|>
  ungroup()|>
  filter(series_length > 5)|>
  rename(mass = "mass_prop", time = "Meas_Day")|>
  filter(cohort_id == 118)

# fit

fit_16 <-fit_litter(mass.remaining = test_series$mass, time = test_series$time, 
                    iters = 1000, model = "weibull",
                    lower = c(0.1,0.001), #beta, then alpha 
                    upper = c(10000,100)) 

fit_16_dp <-fit_litter(mass.remaining = test_series$mass, time = test_series$time, 
                    iters = 1000, model = "discrete.parallel") 
plot(fit_16)

fit_16_cov <- solve(fit_16$optimFit$hessian/fit_16$optimFit$value)

#fit_16_cov2 <- solve(hess_weibull(fit_16)/fit_16$optimFit$value)

# we are minimizing (sum of squares) here, so don't negate Hessian
fit_16_coefs <- fit_16$optimFit$par

beta_hats <- MASS::mvrnorm(9999,fit_16_coefs,fit_16_cov)|>
  as_tibble()|>rename(beta = "V1", alpha = "V2")



#predict

predict_weibull <- function(time, beta, alpha){
  exp(-(time/beta)^alpha)
}

blah <- beta_hats|>
  slice_sample(n = 1000)|>
  rowwise()|>
  mutate(time = list(seq(0,max(test_series$time), by = 0.2)),
         predictions = list(predict_weibull(time = unlist(time), beta = beta, alpha = alpha)))|>
  unnest(cols = c(time, predictions))

ggplot(blah, aes(x = time, y = predictions))+
  geom_line(aes(group = beta))+
  geom_point(data = test_series, aes(x = time, y = mass), color = "red")
