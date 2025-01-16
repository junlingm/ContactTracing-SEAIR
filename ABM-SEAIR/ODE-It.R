#SEAIR three BETA
library(tikzDevice)
library(ggplot2)
library(deSolve)

source("SEAIR.R")
source("ABM-SEAIR/parameters.R")

seair_model <- function(params, times,fixed) {
  # the initial values of variables
  I0 = params[[6]]
  initial_values <- c( S=fixed[["N"]] - I0*3, E=I0, A=I0, I=I0, QE=0, QA = 0,
                       EI=0,EA=0,ET=0,AI = 0,IA =0,II =0,
                       AT = 0,TA = 0,
                       IT=0,TI=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0, sym=0)
  # the parameters values:
  parameters_values <- c(c(betaI  = params[["betaI"]],
                           betaA  = params[["betaA"]],
                           theta = params[["theta"]],
                           p = params[["p"]],
                           tauI = params[["tauI"]]), fixed)
  
  
  # solving
  sol <- ode(initial_values, times, seair_equations, parameters_values)
  sol[,"I"]
}

fixed <- c(sigma = sigma,
           gammaA= gammaA,
           gammaI= gammaI, 
           q= q,
           N=N)


t = 0:250
I <- seair_model(c(betaI=beta,betaA = betaA, theta =theta, p = p, tauI = tauI, I0=I0),
                           times = t,fixed = fixed)

sim.data <- read.csv("mean_I_quantile_results.csv")
plot.data <- rbind(
  sim.data,
  data.frame(Time = 0:250, value = I, curve = "ODE")
)

fig <- ggplot(plot.data, aes(x=Time,y=value, color=curve)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Time",
       y = "I(t)")

print(fig)
