library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)


seair_model <- function(params, times,fixed) {
  # the differential equations:
  seair_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      betaI <- ifelse(time > 32, betaI3, ifelse(time >= 15, betaI2, betaI1))
      tauI <- ifelse(time > 32, tauI3, ifelse(time >= 15, tauI2, tauI1))
      
      betaA <- betaI/3
      
 
      dS <- -betaI*S*I/N-betaA*S*A/N
      dE <-  betaA*S*A/N+betaI*S*I/N-sigma*E-theta*p*ET
      dA <-  sigma*q*E-theta*p*AT-theta*p*TA-gammaA*A
      dI <-  sigma*(1-q)*E-(gammaI+tauI)*I-theta*p*IT-theta*p*TI
      dQE <- theta*p*ET-sigma*q*QE-sigma*(1-q)*QE
      dQA <- theta*p*AT+theta*p*TA+sigma*q*QE-gammaA*QA
      dEI <- betaI*S*I/N-(sigma+tauI+gammaI)*EI-theta*p*EI*(IT/I+TI/I)
      dEA <- betaA*S*A/N-(sigma+gammaA)*EA
      dET <- tauI*EI+theta*p*EI*(IT/I+TI/I)-(sigma+theta)*ET
      dAI <- sigma*q*EI-(gammaA+gammaI+tauI)*AI-theta*p*(AI*TI/I+AI*IT/I+AI*TA/A)
      dIA <- sigma*(1-q)*EA-(gammaA+gammaI+tauI)*IA-theta*p*(IA*TA/A+IA*AT/A+IA*TI/I)
      dII <- sigma*(1-q)*EI-2*(gammaI+tauI)*II-theta*p*II*(IT/I+TI/I+TI/I)
      dAT <- theta*p*AI*(IT/I+TI/I)+tauI*AI+sigma*q*ET-(theta+gammaA)*AT-theta*p*TA*AT/A
      dTA <- theta*p*TI*IA/I+tauI*IA-(theta+gammaA)*TA-theta*p*(TA*AT/A+TA*TA/A)
      dIT <- sigma*(1-q)*ET+theta*p*II*(IT/I+TI/I)+tauI*II-(tauI+theta+gammaI)*IT-theta*p*TI*IT/I
      dTI <- theta*p*TI*II/I+tauI*II-(tauI+theta+gammaI)*TI-theta*p*(TI*IT/I+TI*TI/I)
      dT <-  tauI*I+theta*p*(IT+TI)+sigma*(1-q)*QE-theta*T
      dX <-  theta*T
      dR <- gammaA*A+gammaI*I+gammaA*QA
      dnew <-  tauI*I+theta*p*(IT+TI)+sigma*(1-q)*QE
      dcon <-  theta*p*(IT+TI)+sigma*(1-q)*QE
      dvol <-  tauI*I
      dsym <- sigma*(1-q)*E
      
      return(list(c(dS, dE,dA,dI,dQE, dQA,dEI,dEA,dET,dAI,dIA,dII,dAT,dTA,dIT,dTI,dT,dX,dR,dnew,dcon,dvol,dsym)))
    })
  }
  
  
  # the initial values of variables
  
  I0 = params[["I0"]]*1000
  parameters_values = c(params, fixed)
  initial_values <- c( S=fixed[["N"]] - 3*I0, E=I0, A=I0, I=I0, 
                       QE=0, QA = 0,
                       EI=0,EA=0,ET=0,AI = 0,IA =0,II =0,
                       AT = 0,TA = 0,
                       IT=0,TI=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0,sym=0 )
  
  
  # solving
  ode(initial_values, times, seair_equations, parameters_values)
}


fixed <- c( gammaA= 0.2,
            gammaI= 0, 
            q= 0.3,
            N=14726022
)




data <- read.csv("3.16-5.1-3beta.csv", header = TRUE)



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


