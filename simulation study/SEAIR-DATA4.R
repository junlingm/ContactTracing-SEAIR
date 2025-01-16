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
  observations3<- data$contact
  predictions3<- diff(model_output[,"con"])
  observations4<- data$voluntary
  predictions4<- diff(model_output[,"vol"])
  observations5<- data$QT
  predictions5<- diff(model_output[,"QT"])
  observations6<- data$symptom
  predictions6<- diff(model_output[,"sym"])
  
  
  d1<-dpois(x=observations1,lambda=predictions1, log = TRUE)
  d3<-dpois(x=observations3,lambda=predictions3, log = TRUE)
  d4<-dpois(x=observations4,lambda=predictions4, log = TRUE)
  d5<-dpois(x=observations5,lambda=predictions5, log = TRUE)
  d6<-dpois(x=observations6,lambda=predictions6, log = TRUE)
  log_likelihood <- sum(d1)+ sum(d3)+ sum(d4)+ sum(d5)+ sum(d6)
  return(log_likelihood)
}

parameters_initial <- c(betaI1 =  0.6050257,betaI2 = 0.3333691, betaI3 = 0.4489994,  theta=2.1198038, p=0.2066822,tauI=0.1628169,I0=30)


mcmc_samples <- metrop(obj=joint_log_prob,
                       parameters_initial,
                       nbatch = 200000,
                       blen = 1,
                       scale = 0.0011,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))
mcmc_samples <- metrop(obj=mcmc_samples,
                       nbatch = 250000,
                       blen = 1,
                       scale = 0.0011,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))

