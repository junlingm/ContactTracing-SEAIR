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
                       new=0,con=0,vol=0,sym=0,QT=0)
  # the parameters values:
  parameters_values <- c(c(betaI1  = params[["betaI1"]],
                           betaI2  = params[["betaI2"]],
                           betaI3  = params[["betaI3"]],
                           theta = params[["theta"]],
                           p = params[["p"]]), fixed)
  
  # solving
  ode(initial_values, times, seair_equations_3beta, parameters_values)
}

# read simulated epidemic data
data <- read.csv("simulation study/T80-Nlarge-3beta-2.csv", header = TRUE)

fixed <- c(sigma = 0.27,
           gammaA= 0.2,
           gammaI= 0.1, 
           q= 0.3,
           tauI=0.15,
           N=300000)

prior <- function(parameters) {
  if(any(parameters < 0))
    return(-Inf)
  betaI1_prior = dnorm(parameters[["betaI1"]], mean=0.55, sd=0.2, log = T)
  betaI2_prior = dnorm(parameters[["betaI2"]], mean=0.25, sd=0.2, log = T)
  betaI3_prior = dnorm(parameters[["betaI3"]], mean=0.4, sd=0.2, log = T)
  theta_prior = dnorm(parameters[["theta"]], mean=1.5, sd=1, log = T)
  p_prior = dunif(parameters[["p"]], 0, 1, log = T)
  tauI_prior = if (is.na(parameters["tau"])) 0 else
    dnorm(parameters[["tauI"]], mean=0.15, sd=0.1, log = T)
  I0_prior = dnorm(parameters[["I0"]], mean=1, sd=1, log = T)
  return(betaI1_prior+betaI2_prior+betaI3_prior+theta_prior+p_prior+tauI_prior+I0_prior)
}

joint_log_prob <- function(parameters, fixed, data, times=c(0, data$t)) {
  names(parameters) <- names(parameters_initial)
  log_likelihood <- likelihood(parameters, fixed=fixed, data=data, times=times)
  log_prior <- prior(parameters)
  return(log_likelihood + log_prior)
}

