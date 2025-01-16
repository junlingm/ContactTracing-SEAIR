library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)

source("simulation study/SEIR-3beta.R")

likelihood <- function(parameters, fixed, data, times=c(0, data$t)) {
  if (any(parameters<0))
    return(-Inf)
  model_output <- seir_model( times = times,  params = parameters, fixed = fixed)  
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


parameters_initial <- c(betaI1 = 0.4729350,betaI2 =  0.2635146, betaI3 = 0.3439327,  theta=1.8214296, p=0.1688085,tau=0.1553022, I0=32.910579)

mcmc_samples <- metrop(obj=joint_log_prob,
                       parameters_initial,
                       nbatch =300000,
                       blen = 1,
                       scale = 0.001,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))

mcmc_samples <- metrop(obj=mcmc_samples,
                       nbatch = 500000,
                       blen = 1,
                       scale = 0.001,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))

