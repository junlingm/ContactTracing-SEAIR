rm(list = ls())

library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)


seair_model <- function(params, times,fixed) {
  # the differential equations:
  seair_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      if(time >60) betaI<-betaI3
      else if (time>=40) betaI<-betaI2
      else betaI<-betaI1
      
      betaA=betaI/3
      
      dS <- -betaI*S*I/N-betaA*S*A/N
      dE <-  betaA*S*A/N+betaI*S*I/N-sigma*E-theta*p*ET
      dA <-  sigma*q*E-theta*p*AT-theta*p*TA-gammaA*A
      dI <-  sigma*(1-q)*E-(gammaI+tauI)*I-theta*p*IT-theta*p*TI
      dQE <- theta*p*ET-sigma*q*QE-sigma*(1-q)*QE
      dQA <- theta*p*AT+theta*p*TA+sigma*q*QE-gammaA*QA
      dEI <- betaI*S*I/N-(sigma+tauI+gammaI)*EI-theta*p*EI*(IT/I+TI/I)
      dEA <- betaA*S*A/N-(sigma+gammaA)*EA-theta*p*EA*(AT/A+TA/A)
      dET <- tauI*EI+theta*p*EI*(IT/I+TI/I)-(sigma+theta)*ET
      dAI <- sigma*q*EI-(gammaA+gammaI+tauI)*AI-theta*p*(AI*TI/I+AI*IT/I+AI*TA/A)
      dIA <- sigma*(1-q)*EA-(gammaA+gammaI+tauI)*IA-theta*p*(IA*TA/A+IA*AT/A+IA*TI/I)
      dII <- sigma*(1-q)*EI-2*(gammaI+tauI)*II-theta*p*II*(IT/I+TI/I+TI/I)
      dAT <- theta*p*AI*(IT/I+TI/I)+tauI*AI+sigma*q*ET-(theta+gammaA)*AT-theta*p*TA*AT/A
      dTA <- theta*p*TI*IA/I+tauI*IA-(theta+gammaA)*TA-theta*p*(TA*AT/A+TA*TA/A)
      dIT <- sigma*(1-q)*ET+theta*p*II*(IT/I+TI/I)+tauI*II-(tauI+theta+gammaI)*IT-theta*p*TI*IT/I
      dTI <- theta*p*TI*II/I+tauI*II-(tauI+theta+gammaI)*TI-theta*p*(TI*IT/I+TI*TI/I)
      dT <-  tauI*I+theta*p*(IT+TI)+sigma*(1-q)*QE-theta*T;
      dX <-  theta*T;
      dR <- gammaA*A+gammaI*I+gammaA*QA
      dnew <-  tauI*I+theta*p*(IT+TI)+sigma*(1-q)*QE
      dcon <-  theta*p*(IT+TI)
      dQT <- sigma*(1-q)*QE
      dvol <-  tauI*I
      dsym <- sigma*(1-q)*E
      return(list(c(dS, dE,dA,dI,dQE, dQA,dEI,dEA,dET,dAI,dIA,dII,dAT,dTA,dIT,dTI,dT,dX,dR,dnew,dcon,dQT,dvol,dsym)))
    })
  }
  
  
  # the initial values of variables
  
  I0 = params[[6]]*10
  initial_values <- c( S=fixed[["N"]] - 3*I0, E=I0, A=I0, I=I0, QE=0, QA = 0,
                       EI=0,EA=0,ET=0,AI = 0,IA =0,II =0,
                       AT = 0,TA = 0,
                       IT=0,TI=0,T=0,X=0,R=0,
                       new=0,con=0,QT=0,vol=0,sym=0 )
  # the parameters values:
  parameters_values <- c(c(betaI1  = params[[1]],
                           betaI2  = params[[2]],
                           betaI3  = params[[3]],
                           theta = params[[4]],
                           p = params[[5]]), fixed)
  
  
  # solving
  ode(initial_values, times, seair_equations, parameters_values)
}


fixed <- c(sigma = 0.27,
           gammaA= 0.2,
           gammaI= 0.1, 
           q= 0.3,
           tauI=0.15,
           N=300000)



data <- read.csv("T80-Nlarge-3beta-2.csv", header = TRUE)

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
  
  
  d1<-dpois(x=observations1,lambda=predictions1, log = TRUE)
  d3<-dpois(x=observations3,lambda=predictions3, log = TRUE)
  d4<-dpois(x=observations4,lambda=predictions4, log = TRUE)
  d5<-dpois(x=observations5,lambda=predictions5, log = TRUE)
  log_likelihood <- sum(d1)+ sum(d3)+ sum(d4)+ sum(d5)
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


parameters_initial <- c(betaI1 =  0.6050257,betaI2 = 0.3333691, betaI3 = 0.4489994,  theta=2.1198038, p=0.2066822,I0=2.9976894)

joint_log_prob <- function(parameters, fixed, data, times=c(0, data$t)) {
  log_likelihood <- likelihood(parameters, fixed=fixed, data=data, times=times)
  log_prior <- prior(parameters)
  return(log_likelihood + log_prior)
}



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

