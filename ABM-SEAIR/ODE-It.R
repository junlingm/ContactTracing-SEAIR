#SEAIR three BETA
library(tikzDevice)
library(ggplot2)
library(deSolve)

seair_model <- function(params, times,fixed) {
  # the differential equations:
  seair_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
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
      dcon <-  theta*p*(IT+TI)+sigma*(1-q)*QE
      dvol <-  sigma*(1-q)*QE
      
      return(list(c(dS, dE,dA,dI,dQE, dQA,dEI,dEA,dET,dAI,dIA,dII,dAT,dTA,dIT,dTI,dT,dX,dR,dnew,dcon,dvol)))
    })
  }
  
  
  # the initial values of variables
  
  I0 = params[[6]]
  initial_values <- c( S=fixed[["N"]] - I0*3, E=I0, A=I0, I=I0, QE=0, QA = 0,
                       EI=0,EA=0,ET=0,AI = 0,IA =0,II =0,
                       AT = 0,TA = 0,
                       IT=0,TI=0,T=0,X=0,R=0,
                       new=0,con=0,vol=0 )
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
