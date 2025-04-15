fLD <- function(beta,p){ (log(p/(1-p))-beta[1])/beta[2] }  # function to calculate LD(p)

n <- 20                        # number of dose-levels
m <- 20                        # number of individuals per dose
p <- 0.5                       # type of LD to estimate (e.g. LD(0.5) = LD50)
b0 <- -3                       # parameters in the true model
b1 <- 6
nsim <- 10^2                   # number of simulations
steps <- 10^3                  # how often to print progress
B <- 10^4                      # number of bootstrap samples
conf <- 0.95                   # confidence level

### derived quantities
alpha <- (1-conf)/2
z <- qnorm(1-alpha)
LD <- (log(p/(1-p))-b0)/b1     # true LD(p)
x <- seq(0,1,length.out=n)     # set dose levels
pi <- 1/(1+exp(-(b0+b1*x)))    # true probabilities of death

### used to store results
rd.low <- rd.upp <- array(NA,c(nsim,2))
wald.low <- wald.upp <- vector()
sfb.low <- sfb.upp <- vector()
pb.low <- pb.upp <- vector()
LDs <- LD.sim.pb <- vector()

### simulation loop
for(j in 1:nsim){
  
  # generate the data
  died <- rbinom(n,m,pi)
  y <- cbind(died,m-died)
  
  
  # fit the model
  a <- glm(y~x,family=binomial)
  
  # calculate wald interval for LD(p)
  beta.hat <- coef(a)
  cov.hat <- vcov(a)
  LD.hat <- fLD(beta.hat,p)
  dg <- grad(fLD,beta.hat,p=p)
  se.hat <- sqrt(t(dg)%*%cov.hat%*%dg)
  wald.low[j] <- LD.hat-z*se.hat           # wald confidence limits for LD(p)
  wald.upp[j] <- LD.hat+z*se.hat
  
  # compare wald and profile-likelihood limits for each element of beta  
  se.beta.hat <- sqrt(diag(cov.hat))
  plik <- confint(profile(a),level=conf)
  w.low <- beta.hat-z*se.beta.hat          # wald limits for each element of beta
  w.upp <- beta.hat+z*se.beta.hat 
  p.low <- plik[,1]                        # profile-likelihood limits for each element of beta
  p.upp <- plik[,2]
  w.wid <- w.upp-w.low
  p.wid <- p.upp-p.low
  rd.low[j,] <- abs((w.low-p.low)/p.wid)   # compare wald and profile-likelihood limits
  rd.upp[j,] <- abs((w.upp-p.upp)/p.wid)
  
  # calculate SFB interval (if wald and profile-likelihood limits for beta similar)
  beta.hats <- MASS::mvrnorm(B,beta.hat,cov.hat)            # generate B estimates of beta
  for (i in 1:B){ LDs[i] <- fLD(beta.hats[i,],p) }    # corresponding estimates of LD(P)
  sfb.low[j] <- quantile(LDs,probs=alpha)             # SFB limits for LD(p)
  sfb.upp[j] <- quantile(LDs,probs=1-alpha)
  
  # calculate PB interval
  pi.hat <- fitted(a)                                 # estimates of probabilities of death
  for(i in 1:B){                                    
    died.sim.pb <- rbinom(n,m,pi.hat)                 # parametric simulation of new data
    y.sim.pb <- cbind(died.sim.pb,m-died.sim.pb)
    a.sim.pb <- glm(y.sim.pb~x,family=binomial)       # fit model to new data
    beta.hat.sim.pb <- coef(a.sim.pb)                 
    
    LD.sim.pb[i] <- fLD(beta.hat.sim.pb,p)            # calculate estimate of LD(p)
  }  
  pb.low[j] <- quantile(LD.sim.pb,probs=alpha)        # PB limits for LD(p)
  pb.upp[j] <- quantile(LD.sim.pb,probs=1-alpha)
  
  # print progress
  if(j%%steps==0) { print(j) }
  
}
