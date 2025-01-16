library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)


seir_model <- function(params, times,fixed) {
  # the differential equations:
  seir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      beta <- ifelse(time > 32, betaI3, ifelse(time >= 15, betaI2, betaI1))
      tau <- ifelse(time > 32, tauI3, ifelse(time >= 15, tauI2, tauI1))

      
      
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
  
  I0 = params[["I0"]]*1000
  parameters_values = c(params, fixed)
  
  initial_values <- c( S=fixed[["N"]] - 2*I0, E=I0, I=I0, 
                       EI=0,ET=0,II =0,
                       IT=0,TI=0,Q=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0,sym=0 )
  
  
  # solving
  ode(initial_values, times, seir_equations, parameters_values)
}


fixed <- c(gamma= 0, 
           q= 0.3,
           N=14726022
)





data <- read.csv("3.16-5.1-3beta.csv", header = TRUE)



likelihood <- function(parameters, fixed, data, times=c(0, data$t)) {
  if (any(parameters<0))
    return(-Inf)
  model_output <- seir_model( times = times,  params = parameters, fixed = fixed)  # 求解微分方程
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
    theta_prior = dnorm(theta1, mean=0.5, sd=0.2, log = T)
    p_prior = dunif(p, 0,1 , log = T)
    sigma_prior = dnorm(sigma, mean=0.2, sd=0.1, log = T)
    tauI1_prior = dnorm(tauI1, mean=0.2, sd=0.1, log = T)
    tauI2_prior = dnorm(tauI2, mean=0.2, sd=0.1, log = T)
    tauI3_prior = dnorm(tauI3, mean=0.2, sd=0.1, log = T)
    I0_prior = dnorm(I0, mean=1, sd=0.4, log = T)
    return(betaI1_prior+betaI2_prior+betaI3_prior+theta_prior+p_prior+sigma_prior+tauI1_prior+tauI2_prior+tauI3_prior+I0_prior)
  })
}



parameters_initial <- c(betaI1= 0.5, betaI2= 0.35,betaI3=0.3, theta1 = 1,p=0.3, sigma= 0.3,tauI1=0.1,tauI2= 0.15,tauI3=0.2,I0=1)



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
                       scale = 0.003,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))

mcmc_samples <- metrop(obj=mcmc_samples,
                       nbatch = 200000,
                       blen = 10,
                       scale = 0.003,
                       debug = FALSE,
                       data = data,
                       fixed=fixed,
                       times = c(0, data$t))

# 
# sample1 = mcmc_samples
# save(sample1, file="sample1_SEIR_On-3.16-5.1.RData")

