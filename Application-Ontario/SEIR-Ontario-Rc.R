library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)


seir_model <- function(times,fixed) {
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
  
  I0 <- fixed[["I0"]] 
  parameters_values <- fixed
  
  initial_values <- c( S=fixed[["N"]] - 2*I0, E=I0, I=I0, 
                       EI=0,ET=0,II =0,
                       IT=0,TI=0,Q=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0,sym=0 )
  
  
  # solving
  ode(initial_values, times, seir_equations, parameters_values)
}



fixed <- c(betaI1=0.22470651, 
           betaI2=0.20214211,
           betaI3=  0.14361733, 
           theta =  0.89553045, 
           p = 0.26991923, 
           sigma= 0.57844540, 
           tauI1= 0.08081608, 
           tauI2= 0.12328228, 
           tauI3= 0.13013146, 
           I0=  362.098110,
           gamma= 0, 
           q= 0.3,
           N=14726022
)




data <- read.csv("3.16-5.1-3beta.csv", header = TRUE)

predictions <- seir_model( times = c(0,data$t),fixed = fixed)

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
  beta <- ifelse(t > 32, as.numeric(fixed[["betaI3"]]), 
                  ifelse(t >= 15, as.numeric(fixed[["betaI2"]]), as.numeric(fixed[["betaI1"]])))
  tau <- ifelse(t > 32, as.numeric(fixed[["tauI3"]]), 
                 ifelse(t >= 15, as.numeric(fixed[["tauI2"]]), as.numeric(fixed[["tauI1"]])))
  
  # Ensure parameters are numeric
  sigma <- as.numeric(fixed[["sigma"]])
  theta <- as.numeric(fixed[["theta"]])
  p <- as.numeric(fixed[["p"]])
  gamma <- as.numeric(fixed[["gamma"]])
  

  b <- beta * sigma / (theta * p * vIT[t] + theta * p * vTI[t] + gamma + tau)
  c <- 1 / (sigma + theta * p * uET[t])
  
  # Calculate Rc for current time step
  Rc[t] <- b * c
}

# Plot Rc over time
time_points <- 1:length(Rc)
plot(time_points, Rc, type = "l", xlab = "Time", ylab = "Rc", main = "Time-Dependent Rc")


