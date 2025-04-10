# Fit Bayesian litter fits

# Load data
litter_df <- read_csv("https://github.com/robbinscalebj/revisiting-k/raw/refs/heads/main/data/derived_data/tidied_detrital_nutrients_data.csv")|> # tidied litter time series
  mutate(mass_prop = Mass_per_remaining/100)|>
  group_by(cohort_id)|>
  mutate(series_length = n())|>
  ungroup()|>
  filter(series_length >6)|>
  rename(mass = "mass_prop", time = "Meas_Day")

# Set test data
series_ids <- unique(litter_df$cohort_id)
series_ids <- series_ids[1:4]


# Set priors and formulae

## K

k_bf <- bf(mass ~ exp(-k*time), k ~1, nl = TRUE)
k_priors <- prior(lognormal(-6,2), nlpar = "k", lb = 0)

## Weibull
weibull_bf <- bf(mass ~ exp(-(time/beta)^alpha), beta + alpha ~1, nl = TRUE)
weibull_priors <- prior(gamma(1,7), nlpar = "alpha", lb = 0)+
  prior(lognormal(5.5,1), nlpar = "beta", lb = 0)

## Discrete Series
discser_bf <- bf(mass ~ a*exp(-k1 * time) + (1 - a) * exp(-k2 * time),
                 a + k1 + k2 ~ 1, nl = TRUE)



# Fitting loops
nested_data <- litter_df|>
  filter(cohort_id %in% series_ids)|>
  select(cohort_id, time, mass)|>
  #filter(cohort_id %in% seq(1:38))|> #for testing
  group_by(cohort_id)|> # needs group_by()|>nest() structure... nest_by() doesn't work!
  nest()|>
  ungroup()|>
  mutate(batch_number = row_number())

# Definitions for use in parallel loops
n_data_batches <- nrow(nested_data)
cohort_names <- nested_data|>pull(cohort_id)

pb <- txtProgressBar(max = n_data_batches, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Set parallel processing
cl <- makePSOCKcluster(4)
registerDoSNOW(cl)
tic()
# Loop
weibull_summaries <- foreach(i = 1:nrow(nested_data), .combine = 'rbind', .options.snow = opts,
                         .packages = c("tidyverse", "brms", "here", "tidybayes"), .inorder = FALSE
) %dopar% {

  # filter to cohort data and unnest into normal data frame
  df <- nested_data|>filter(batch_number == i)|>unnest()
  cohort_id <- df|>first()|>pull(cohort_id)
  file_name <- paste0("cohort_", cohort_id)

  fit <- brm(weibull_bf,
             prior = weibull_priors,
             data = df,
             chains = 4, cores = 1, iter = 4000,
             warmup = 500, control= list(adapt_delta = 0.99),
             file = paste0("analysis/data/derived_data/brms_fits/weibull/",
                           file_name)
  )

  summarise_draws(fit)|>
    tibble(cohort_id = cohort_id,
           model_type = "Weibull")|>
      relocate(cohort_id, model_type)


}
stopCluster(cl)
toc()
# save model objects, summarize parameter HDIs, diagnostic values, to table
fit <- readRDS(here("analysis/data/derived_data/brms_fits/weibull/cohort_16.Rds"))
summary(fit)
summarise_draws(fit)


# use this to compare with litterfitter output. Need to adjust litterfitter to produce all
# desired posterior estimands
