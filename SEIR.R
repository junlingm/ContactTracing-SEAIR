seir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
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
    dQT<- sigma*Q
    
    return(list(c(dS, dE,dI,dEI,dET,dII,dIT,dTI,dQ,dT,dX,dR,dnew,dcon,dvol,dsym,dQT)))
  })
}
