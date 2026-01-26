#functions that work with posterior predictions

#generate predictions
# This function uses marginaleffects::predictions() to generate predictions for each timestep (default 1 day) covering the time series in the model's data
# get_draws controls whether posterior draws of the predictions are provided. If FALSE, predictions() summarizes all draws to upper and lower 95% credible 
# intervals
# if get_draws == TRUE, trim_every chooses the number of draws to skip to reduce data size and processing time for for downstream analysis - defaults to all draws
# timestep controls number of predictions made evenly from day 0 to the last observation in the day. Units are days, default is 1 (generate prediction for each day) 


generate_predictions <- function(model, trim_every = 1, get_draws = FALSE, timestep=1){
  
  if(get_draws == TRUE){
  predictions(model,
              newdata = datagrid(time = \(x) seq(0,max(x), by = timestep)))|>
    get_draws()|>
    mutate(drawid = as.numeric(drawid))|>
    filter(drawid %% trim_every == 0)
  }else{
    predictions(model,
                newdata = datagrid(time = \(x) seq(0,max(x), by = timestep)))
  }
}





# function to numerically estimate derived quantity time-to-mass loss (e.g., t0.5 - half-life)
# takes existing E[t] predictions from generate_predictions, assuming get_draws is true
# M_x is proportion of mass remaining to estimate time function crosses it - e.g., 0.5 would be time to hit 0.5 
generate_tx <- function(predictions, M_x){
  

  if(!("draw" %in% names(predictions))){
    stop("Draws are missing - use generate_predictions() with get_draws = TRUE and a trim value less than the total number of posterior draws")
  }

  
  n_draws <- length(unique(predictions$drawid))
  est_timestep <- (max(predictions$time)-min(predictions$time))/length(unique(predictions$time))
  

  if(n_draws < 1000){
    message("Number of draws is less than 1000 - consider whether credible interval tails are stable")
  }
  
  if(est_timestep >1.09){
    message(paste("Timestep is ~", round(est_timestep, digits = 1), " days. Timestep may be too large for accurately identifying bisection point"))
  }

t_xf0_name <- paste0("t_",M_x)
    predictions|>
      group_by(rowid)|>
      mutate(fx_diff = draw - M_x)|>
      arrange(drawid, time)|>
      group_by(drawid)|>
      mutate(sign_change = fx_diff * lag(fx_diff) < 0)|>
      filter(sign_change | lead(sign_change, default = FALSE))|>
      mutate("{t_xf0_name}" := time + (time - lag(time)) * (-fx_diff) / (fx_diff - lag(fx_diff)))|>
      filter(sign_change == TRUE)
  }


#generate_tx(summed_preds)
blah <- generate_tx(draws_preds, M_x = 0.8)


ggplot(blah, aes(x = t_0.8))+geom_density()
    
  
  


# function to generate kapp along function
generate_kapp_preds <- function(predictions, est_lambda = FALSE){
  
  #if est_lambda is TRUE, evaluate dataframe for a coarse/fine mesh variable, to be ID'd in args?
  
  if(!("draw" %in% names(predictions))){
    stop("Draws are missing - use generate_predictions() with get_draws = TRUE and a trim value less than the total number of posterior draws")
  }
  
  n_draws <- length(unique(predictions$drawid))
  est_timestep <- (max(predictions$time)-min(predictions$time))/length(unique(predictions$time))
  
  
  if(n_draws < 1000){
    message("Number of draws is less than 1000 - consider whether credible interval tails are stable")
  }
  
  if(est_timestep >1.09){
    message(paste("Timestep is ~", round(est_timestep, digits = 1), " days. Timestep may be too large for accurately identifying instantaneous slopes"))
  }
  
 predictions|>ungroup()|>
   group_by(drawid)|>
   arrange(draw, time)|>
   mutate(logmass = log(draw),
          timestep = time-lead(time),
          k_app = (logmass-lead(logmass))/timestep)
  
}


draws_preds <- generate_predictions(pnge_436_brm, trim_every = 10, get_draws = TRUE, timestep = 1)

blah2 <- generate_kapp_preds(draws_preds)

# function to generate lambda microbe and lambda invert along function