# Prediction functions
### modified from litterfitter functions of Cornwell and Weedon
### https://github.com/traitecoevo/litterfitter/blob/master/R/internal.R


## negative exponential
predict_negexp <- function(time, k){
  exp(-k*time)
}

## discrete series
predict_discser <- function(time, r, k1, k2){
  (((1 - r) * k1 * exp(-k2 * time)) - ((k2 - k1 * r) * exp(-k1 * time)))/(k1 - k2)
}

## discete parallel
predict_discpar <- function(time, a, k1, k2){
  a * exp(-k1 * time) + (1 - a) * exp(-k2 * time)
}
## continuous quality
predict_contqual <- function(time, a, b) {
  1/((1 + b * time)^a)
}
## weibull
predict_weibull <- function(time, beta, alpha){
  exp(-(time/beta)^alpha)
}


newtimes <- seq(0,300, by = 1)
