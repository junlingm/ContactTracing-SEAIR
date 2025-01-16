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
    dsym <- sigma*(1-q)*E
    
    return(list(c(dS, dE,dA,dI,dQE, dQA,dEI,dEA,dET,dAI,dIA,dII,dAT,dTA,dIT,dTI,dT,dX,dR,dnew,dcon,dvol,dsym)))
  })
}
