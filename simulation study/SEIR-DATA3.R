
rm(list = ls())

library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)

seir_model <- function(params, times,fixed) {
  # the differential equations:
  seir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      
      if(time >60) beta<-betaI3
      else if (time>=40) beta<-betaI2
      else beta<-betaI1
      
      
      
      dS <- -beta*S*I/N
      dE <- beta*S*I/N-sigma*E-theta*p*ET
      dI<- sigma*E-(gamma+tau)*I-theta*p*IT-theta*p*TI
      dEI<- beta*S*I/N-(sigma+tau+gamma)*EI-theta*p*EI*(IT/I+TI/I)
      dET<- tau*EI+theta*p*EI*(IT/I+TI/I)-(sigma+theta)*ET
      dII<- sigma*EI-2*(gamma+tau)*II-theta*p*II*(IT/I+TI/I+TI/I)
      dIT<-sigma*ET+theta*p*II*(IT/I+TI/I)+tau*II-(tau+theta+gamma)*IT-theta*p*TI*IT/I
      dTI<- theta*p*TI*II/I+tau*II-(tau+theta+gamma)*TI-theta*p*(TI*IT/I+TI*TI/I)
      dQ<-theta*p*ET-sigma*Q
      dT<-sigma*Q+tau*I+theta*p*(IT+TI)-theta*T
      dX<-theta*T
      dR<-gamma*I
      dnew <-tau*I+theta*p*(IT+TI)+sigma*Q
      dcon<-theta*p*(IT+TI)+sigma*Q
      dvol<- tau*I
      dsym <- sigma*E
      
      
      return(list(c(dS, dE,dI,dEI,dET,dII,dIT,dTI,dQ,dT,dX,dR,dnew,dcon,dvol,dsym)))
    })
  }
  
  
  # the initial values of variables
  
  I0 = params[[7]]*10
  
  initial_values <- c( S=fixed[["N"]] - 2*I0, E=I0, I=I0, 
                       EI=0,ET=0,II =0,
                       IT=0,TI=0,Q=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0,sym=0 )
  # the parameters values:
  parameters_values <- c(c(betaI1  = params[[1]],
                           betaI2  = params[[2]],
                           betaI3  = params[[3]],
                           theta = params[[4]],
                           p = params[[5]],
                           tau = params[[6]]), fixed)
  
  
  # solving
  ode(initial_values, times, seir_equations, parameters_values)
}


fixed <- c(sigma = 0.27,
           gamma= 0.1, 
           q= 0.3,
           N=14726022
)



data <- read.csv("T80-Nlarge-3beta-2.csv", header = TRUE)


likelihood <- function(parameters, fixed, data, times=c(0, data$t)) {
  if (any(parameters<0))
    return(-Inf)
  model_output <- seir_model( times = times,  params = parameters, fixed = fixed) 
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
  
  log_likelihood <- sum(d1)+sum(d3)+ sum(d4)+ sum(d5)
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
  tau_prior = dnorm(parameters[[6]], mean=0.15, sd=0.1, log = T)
  I0_prior = dnorm(parameters[[7]], mean=1, sd=1, log = T)
  return(betaI1_prior+betaI2_prior+betaI3_prior+theta_prior+p_prior+tau_prior+I0_prior)
}



parameters_initial <- c(betaI1 = 0.4894365,betaI2 = 0.2606752, betaI3 = 0.3484864,  theta=1.9333680, p=0.1693686,tau=0.1578712,I0=2.8164169)


joint_log_prob <- function(parameters, fixed, data, times=c(0, data$t)) {
  log_likelihood <- likelihood(parameters, fixed=fixed, data=data, times=times)
  log_prior <- prior(parameters)
  return(log_likelihood + log_prior)
}


mcmc_samples <- metrop(obj=joint_log_prob,
                       parameters_initial,
                       nbatch = 200000,
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


