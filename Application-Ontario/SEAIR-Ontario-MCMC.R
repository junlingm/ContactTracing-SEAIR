library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)


source("Application-Ontario/SEAIR-3beta-3tau.R")

likelihood <- function(parameters, fixed, data, times=c(0, data$t)) {
  if (any(parameters<0))
    return(-Inf)
  model_output <- seair_model( times = times,  params = parameters, fixed = fixed)  
  observations1<- data$newcase
  predictions1<- diff(model_output[,"new"])
  observations2<- data$symptom
  predictions2<- diff(model_output[,"sym"])
  observations3<- data$CC
  predictions3<- diff(model_output[,"con"])
  
  
  d1<-dpois(x=observations1,lambda=predictions1, log = TRUE)
  d2<-dpois(x=observations2,lambda=predictions2, log = TRUE)
  d3<-dpois(x=observations3,lambda=predictions3, log = TRUE)
  log_likelihood <- sum(d1)+ sum(d2)+ sum(d3)
  return(log_likelihood)
}


prior <- function(parameters) {
  if(any(parameters < 0))
    return(-Inf)
  with(as.list(parameters), {
    betaI1_prior = dnorm(betaI1, mean=0.3, sd=0.2, log = T)
    betaI2_prior = dnorm(betaI2, mean=0.25, sd=0.2, log = T)
    betaI3_prior = dnorm(betaI3, mean=0.15, sd=0.1, log = T)
    theta_prior = dnorm(theta, mean=0.5, sd=0.2, log = T)
    p_prior = dunif(p, 0,1 , log = T)
    sigma_prior = dnorm(sigma, mean=0.2, sd=0.1, log = T)
    tauI1_prior = dnorm(tauI1, mean=0.2, sd=0.1, log = T)
    tauI2_prior = dnorm(tauI2, mean=0.2, sd=0.1, log = T)
    tauI3_prior = dnorm(tauI3, mean=0.2, sd=0.1, log = T)
    I0_prior = dnorm(I0, mean=1, sd=0.4, log = T)
    return(betaI1_prior+betaI2_prior+betaI3_prior+theta_prior+p_prior+sigma_prior+tauI1_prior+tauI2_prior+tauI3_prior+I0_prior)
  })
}




parameters_initial <- c(betaI1= 0.5, betaI2= 0.35,betaI3=0.3, theta= 1, p = 0.2, sigma= 0.3,tauI1=0.1,tauI2= 0.15,tauI3=0.2,I0=1)

############

joint_log_prob <- function(parameters, fixed, data, times=c(0, data$t)) {
  names(parameters) <- names(parameters_initial)
  log_likelihood <- likelihood(parameters, fixed=fixed, data=data, times=times)
  log_prior <- prior(parameters)
  return(log_likelihood + log_prior)
}


mcmc_samples <- metrop(obj=joint_log_prob,
                       parameters_initial,
                       nbatch = 100000,
                       blen = 10,
                       scale = 0.002,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))

mcmc_samples <- metrop(obj=mcmc_samples,
                       nbatch = 200000,
                       blen = 10,
                       scale = 0.002,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))



# sample1 = mcmc_samples
# save(sample1, file="sample1_SEAIR_ON-3.16-5.1.RData")
# 


