source("SEAIR.R")

seair_equations_3beta_3tau <- function(time, variables, parameters) {
  betaI <- ifelse(time > 32, parameters[["betaI3"]], 
                  ifelse(time >= 15, parameters[["betaI2"]], parameters[["betaI1"]]))
  tauI <- ifelse(time > 32, parameters[["tauI3"]], 
                 ifelse(time >= 15, parameters[["tauI2"]], parameters[["tauI1"]]))
  parameters[["betaI"]] <- betaI
  parameters[["betaA"]] <- betaI / 3
  parameters[["tauI"]] <- tauI
  seair_equations(time, variables, parameters)
}


seair_model <- function(params, times,fixed) {
# the initial values of variables
  I0 = params[["I0"]]*1000
  parameters_values = c(params, fixed)
  initial_values <- c( S=fixed[["N"]] - 3*I0, E=I0, A=I0, I=I0, 
                       QE=0, QA = 0,
                       EI=0,EA=0,ET=0,AI = 0,IA =0,II =0,
                       AT = 0,TA = 0,
                       IT=0,TI=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0,sym=0,QT=0)
  # solving
  ode(initial_values, times, seair_equations_3beta_3tau, parameters_values)
}

fixed <- c( gammaA= 0.2,
            gammaI= 0, 
            q= 0.3,
            N=14726022
)

data <- read.csv("Application-Ontario/3.16-5.1-3beta.csv", header = TRUE)
