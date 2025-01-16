#SEAIR three BETA
rm(list = ls())

library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)


seair_model <- function(times,fixed) {
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
  
  I0 <- fixed[["I0"]] 
  parameters_values <- fixed
  initial_values <- c( S=fixed[["N"]] - 3*I0, E=I0, A=I0, I=I0, 
                       QE=0, QA = 0,
                       EI=0,EA=0,ET=0,AI = 0,IA =0,II =0,
                       AT = 0,TA = 0,
                       IT=0,TI=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0,sym=0 )
  
  
  # solving
  ode(initial_values, times, seair_equations, parameters_values)
}


fixed <- c( betaI1=0.29751292, 
            betaI2=0.27205937,
            betaI3= 0.20167404 , 
            theta = 1.15275406, 
            p = 0.29329760, 
            sigma = 0.75168476, 
            tauI1= 0.08540103, 
            tauI2= 0.13119092, 
            tauI3= 0.14093923, 
            I0= 339.222464,
            gammaA= 0.2,
            gammaI= 0, 
            q= 0.3,
            N=14726022
)




data <- read.csv("3.16-5.1-3beta.csv", header = TRUE)

predictions <- seair_model( times = c(0,data$t),fixed = fixed)


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
  betaI <- ifelse(t > 32, as.numeric(fixed[["betaI3"]]), 
                  ifelse(t >= 15, as.numeric(fixed[["betaI2"]]), as.numeric(fixed[["betaI1"]])))
  tauI <- ifelse(t > 32, as.numeric(fixed[["tauI3"]]), 
                 ifelse(t >= 15, as.numeric(fixed[["tauI2"]]), as.numeric(fixed[["tauI1"]])))
  
  # Ensure parameters are numeric
  sigma <- as.numeric(fixed[["sigma"]])
  q <- as.numeric(fixed[["q"]])
  theta <- as.numeric(fixed[["theta"]])
  p <- as.numeric(fixed[["p"]])
  gammaA <- as.numeric(fixed[["gammaA"]])
  gammaI <- as.numeric(fixed[["gammaI"]])
  
  # Calculate a, b, and c for current time step
  a <- (betaI / 3) * sigma * q / (theta * p * vAT[t] + theta * p * vTA[t] + gammaA)
  b <- betaI * sigma * (1 - q) / (theta * p * wIT[t] + theta * p * wTI[t] + gammaI + tauI)
  c <- 1 / (sigma + theta * p * uET[t])
  
  # Calculate Rc for current time step
  Rc[t] <- (a + b) * c
}

# Plot Rc over time
time_points <- 1:length(Rc)
plot(time_points, Rc, type = "l", xlab = "Time", ylab = "Rc", main = "Time-Dependent Rc")


