library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)

source("SEAIR.R")

seair_equations_3beta <- function(time, variables, parameters) {
  if(time >60) betaI<-parameters["betaI3"]
  else if (time>=40) betaI<-parameters["betaI2"]
  else betaI<-parameters["betaI1"]
  parameters[["betaI"]] <- betaI
  parameters[["betaA"]] <- betaI / 3
  seair_equations(time, variables, parameters)
}

seair_model <- function(params, times, fixed) {
  I0 = params[["I0"]]
  initial_values <- c( S=fixed[["N"]] - 3*I0, E=I0, A=I0, I=I0, QE=0, QA = 0,
                       EI=0,EA=0,ET=0,AI = 0,IA =0,II =0,
                       AT = 0,TA = 0,
                       IT=0,TI=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0,sym=0 )
  # the parameters values:
  parameters_values <- c(c(betaI1  = params[["betaI1"]],
                           betaI2  = params[["betaI2"]],
                           betaI3  = params[["betaI3"]],
                           theta = params[["theta"]],
                           p = params[["p"]]), fixed)
  
  # solving
  ode(initial_values, times, seair_equations_3beta, parameters_values)
}


fixed <- c(sigma = 0.27,
           gammaA= 0.2,
           gammaI= 0.1, 
           q= 0.3,
           tauI=0.15,
           N=300000)



# 观测数据


data <- read.csv("simulation study/T80-Nlarge-3beta-2.csv", header = TRUE)


likelihood <- function(parameters, fixed, data, times=c(0, data$t)) {
  if (any(parameters<0))
    return(-Inf)
  names(parameters) <- names(parameters_initial)
  model_output <- seair_model( times = times,  params = parameters, fixed = fixed)  
  observations1<- data$newcase
  predictions1<- diff(model_output[,"new"])
  observations3<- data$conQT
  predictions3<- diff(model_output[,"con"])
  observations4<- data$voluntary
  predictions4<- diff(model_output[,"vol"])
  
  d1<-dpois(x=observations1,lambda=predictions1, log = TRUE)
  d3<-dpois(x=observations3,lambda=predictions3, log = TRUE)
  d4<-dpois(x=observations4,lambda=predictions4, log = TRUE)
  log_likelihood <- sum(d1)+ sum(d3)+ sum(d4)
  return(log_likelihood)
}


prior <- function(parameters) {
  if(any(parameters < 0))
    return(-Inf)
  betaI1_prior = dnorm(parameters[[1]], mean=0.55, sd=0.2, log = T)
  betaI2_prior = dnorm(parameters[[2]], mean=0.25, sd=0.2, log = T)
  betaI3_prior = dnorm(parameters[[3]], mean=0.4, sd=0.2, log = T)
  theta_prior = dnorm(parameters[[4]], mean=1.5, sd=1, log = T)
  p_prior = dunif(parameters[[5]], 0, 1, log = T)
  I0_prior = dnorm(parameters[[6]], mean=1, sd=1, log = T)
  return(betaI1_prior+betaI2_prior+betaI3_prior+theta_prior+p_prior+I0_prior)
}



parameters_initial <- c(betaI1 =  0.6036950,betaI2 = 0.3255494, betaI3 = 0.4452518,  theta= 2.2819232, p= 0.1924505, I0=29)



joint_log_prob <- function(parameters, fixed, data, times=c(0, data$t)) {
  log_likelihood <- likelihood(parameters, fixed=fixed, data=data, times=times)
  log_prior <- prior(parameters)
  return(log_likelihood + log_prior)
}


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