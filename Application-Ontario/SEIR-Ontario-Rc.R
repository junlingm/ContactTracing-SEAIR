library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)

source("Application-Ontario/SEIR-3beta-3tau.R")

parameters <- c(betaI1=0.22470651, 
           betaI2=0.20214211,
           betaI3=  0.14361733, 
           theta =  0.89553045, 
           p = 0.26991923, 
           sigma= 0.57844540, 
           tauI1= 0.08081608, 
           tauI2= 0.12328228, 
           tauI3= 0.13013146, 
           I0=  362.098110
)

predictions <- seir_model(params=parameters, times = c(0,data$t),fixed = fixed)

# Extract predictions
I <- predictions[,"I"]
E <- predictions[,"E"]
IT <- predictions[,"IT"]
TI <- predictions[,"TI"]
ET <- predictions[,"ET"]

# Ensure predictions are numeric and handle divisions safely
vIT <- ifelse(I != 0, IT / I, 0)
vTI <- ifelse(I != 0, TI / I, 0)
uET <- ifelse(E != 0, ET / E, 0)

# Initialize Rc vector
Rc <- numeric(length(I))

# Iterate over each time step to calculate Rc
for (t in seq_along(I)) {
  # Determine betaI and tauI based on the time step
  beta <- ifelse(t > 32, as.numeric(parameters[["betaI3"]]), 
                  ifelse(t >= 15, as.numeric(parameters[["betaI2"]]), as.numeric(parameters[["betaI1"]])))
  tau <- ifelse(t > 32, as.numeric(parameters[["tauI3"]]), 
                 ifelse(t >= 15, as.numeric(parameters[["tauI2"]]), as.numeric(parameters[["tauI1"]])))
  
  # Ensure parameters are numeric
  sigma <- as.numeric(parameters[["sigma"]])
  theta <- as.numeric(parameters[["theta"]])
  p <- as.numeric(parameters[["p"]])
  gamma <- as.numeric(fixed[["gamma"]])
  

  b <- beta * sigma / (theta * p * vIT[t] + theta * p * vTI[t] + gamma + tau)
  c <- 1 / (sigma + theta * p * uET[t])
  
  # Calculate Rc for current time step
  Rc[t] <- b * c
}

# Plot Rc over time
time_points <- 1:length(Rc)
plot(time_points, Rc, type = "l", xlab = "Time", ylab = "Rc", main = "Time-Dependent Rc")


