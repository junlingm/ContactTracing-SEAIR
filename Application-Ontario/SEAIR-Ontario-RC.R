library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)

source("Application-Ontario/SEAIR-3beta-3tau.R")

parameters <- c( betaI1=0.29751292, 
            betaI2=0.27205937,
            betaI3= 0.20167404 , 
            theta = 1.15275406, 
            p = 0.29329760, 
            sigma = 0.75168476, 
            tauI1= 0.08540103, 
            tauI2= 0.13119092, 
            tauI3= 0.14093923, 
            I0= 339.222464
)


predictions <- seair_model(params = parameters, times = c(0,data$t), fixed=fixed)

# Extract predictions
I <- predictions[,"I"]
A <- predictions[,"A"]
E <- predictions[,"E"]
IT <- predictions[,"IT"]
TI <- predictions[,"TI"]
AT <- predictions[,"AT"]
TA <- predictions[,"TA"]
ET <- predictions[,"ET"]

# Ensure predictions are numeric and handle divisions safely
vAT <- ifelse(A != 0, AT / A, 0)
vTA <- ifelse(A != 0, TA / A, 0)
wIT <- ifelse(I != 0, IT / I, 0)
wTI <- ifelse(I != 0, TI / I, 0)
uET <- ifelse(E != 0, ET / E, 0)

# Initialize Rc vector
Rc <- numeric(length(I))

# Iterate over each time step to calculate Rc
for (t in seq_along(I)) {
  # Determine betaI and tauI based on the time step
  betaI <- ifelse(t > 32, as.numeric(parameters[["betaI3"]]), 
                  ifelse(t >= 15, as.numeric(parameters[["betaI2"]]), as.numeric(parameters[["betaI1"]])))
  tauI <- ifelse(t > 32, as.numeric(parameters[["tauI3"]]), 
                 ifelse(t >= 15, as.numeric(parameters[["tauI2"]]), as.numeric(parameters[["tauI1"]])))
  
  # Ensure parameters are numeric
  sigma <- as.numeric(parameters[["sigma"]])
  q <- as.numeric(fixed[["q"]])
  theta <- as.numeric(parameters[["theta"]])
  p <- as.numeric(parameters[["p"]])
  gammaA <- as.numeric(fixed[["gammaA"]])
  gammaI <- as.numeric(fixed[["gammaI"]])
  
  # Calculate a, b, and c for current time step
  a <- (betaI / 3) * sigma * q / (theta * p * vAT[t] + theta * p * vTA[t] + gammaA)
  b <- betaI * sigma * (1 - q) / (theta * p * wIT[t] + theta * p * wTI[t] + gammaI + tauI)
  c <- (sigma + theta * p * uET[t])
  
  # Calculate Rc for current time step
  Rc[t] <- (a + b) / c
}

# Plot Rc over time
time_points <- 1:length(Rc)
plot(time_points, Rc, type = "l", xlab = "Time", ylab = "Rc", main = "Time-Dependent Rc")


