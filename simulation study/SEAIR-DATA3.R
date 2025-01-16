library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)

source("simulation study/SEAIR-3beta.R")


likelihood <- function(parameters, fixed, data, times=c(0, data$t)) {
  if (any(parameters<0))
    return(-Inf)
  model_output <- seair_model( times = times,  params = parameters, fixed = fixed)  
  observations1<- data$newcase
  predictions1<- diff(model_output[,"new"])
  observations3<- data$conQT
  predictions3<- diff(model_output[,"con"])
  observations4<- data$voluntary
  predictions4<- diff(model_output[,"vol"])
  observations5<- data$symptom
  predictions5<- diff(model_output[,"sym"])
  
  d1<-dpois(x=observations1,lambda=predictions1, log = TRUE)
  d3<-dpois(x=observations3,lambda=predictions3, log = TRUE)
  d4<-dpois(x=observations4,lambda=predictions4, log = TRUE)
  d5<-dpois(x=observations5,lambda=predictions5, log = TRUE)
  log_likelihood <- sum(d1)+ sum(d3)+ sum(d4)+ sum(d5)
  return(log_likelihood)
}

parameters_initial <- c(betaI1 =  0.6036950,betaI2 = 0.3255494, betaI3 = 0.4452518,  theta= 2.2819232, p= 0.1924505, tauI=0.1610236,I0=29)

mcmc_samples <- metrop(obj=joint_log_prob,
                       parameters_initial,
                       nbatch = 200000,
                       blen = 1,
                       scale = 0.0012,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))

mcmc_samples <- metrop(obj=mcmc_samples,
                       nbatch = 250000,
                       blen = 1,
                       scale = 0.0012,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))