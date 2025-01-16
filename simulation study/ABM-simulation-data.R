rm(list = ls())
library(ABM)

N = 10000  # the population size
I0 = 20  # the initial number of infections

p = 0.2  # the contact tracing coverage probability
q = 0.3
theta = 2  # the contact tracing rate for diagnosed patients

beta = 0.6  # the contact rate
beta2 = 0.3
beta3 = 0.45
eps2 = beta2 / beta
eps3 = beta3 / beta
T2 = 50
T3 = 70


sigma = 0.27  # the susceptibility, i.e., the probability of transmission per contact

gammaI = 0.1  # the recovery rate
gammaA = 0.2  # the recovery rate
tauI = 0.15  # the diagnose rate
tauA = 0  # the diagnose rate

check_susceptibility = function(time, agent, contact) {
  eps = if (time < T2) 1 else if (time < T3) eps2 else eps3
  runif(1) < eps
}

form_paris = function(time, agent, contact) {
  sa = getState(agent)
  sa$contacts = c(sa$contacts, contact)
  sc = getState(contact)
  sc$contacts = c(sc$contacts, agent)
  setState(agent, sa)
  setState(contact, sc)
}

trace = function(time, agent) {
  for (c in getState(agent)$contacts) {
    if (runif(1) < p) {
      if (matchState(c, "E")) {
        setState(c, "QE")
      } else if (matchState(c, "A")) {
        setState(c, "QA")
      } else if (matchState(c, "I")) {
        setState(c, "Tt")
      }
    }
  }
}

run <- function() {
  sim = Simulation$new(N, function(i) {
    state = if (i <= I0) "I" else if (i <= 2 * I0) "A" else if (i <= 3 * I0) "E" else "S"
    list(state, contacts = list())
  })
  m = newRandomMixing()
  sim$addContact(m)
  
  sim$addTransition("E" -> "I", sigma * (1 - q))
  sim$addTransition("E" -> "A", sigma * q)
  sim$addTransition("A" -> "T", tauA)
  sim$addTransition("A" -> "R", gammaA)
  sim$addTransition("I" -> "R", gammaI)
  sim$addTransition("I" -> "T", tauI)
  sim$addTransition(
    "I" + "S" -> "I" + "E" ~ m, beta,
    to_change_callback = check_susceptibility,
    changed_callback = form_paris
  )
  sim$addTransition(
    "A" + "S" -> "A" + "E" ~ m, beta/3,
    to_change_callback = check_susceptibility,
    changed_callback = form_paris
  )
  sim$addTransition(  # voluntary
    "T" -> "X", theta,
    changed_callback = trace
  )
  sim$addTransition(  # contact tracing
    "Tt" -> "X", theta,
    changed_callback = trace
  )
  sim$addTransition("QA" -> "Tt", tauA)
  sim$addTransition("QA" -> "R", gammaA)
  sim$addTransition("QE" -> "QA", sigma * q)
  sim$addTransition("QE" -> "Tt", sigma * (1 - q))
  
  
  sim$addLogger(newCounter("I", "I"))
  sim$addLogger(newCounter("T", "T"))
  sim$addLogger(newCounter("Tt", "Tt"))
  sim$addLogger(newCounter("VoluntaryI", "I", "T"))
  sim$addLogger(newCounter("QuarantinedE", "QE", "Tt"))
  sim$addLogger(newCounter("onset", "E", "I"))
  sim$addLogger(newCounter("TracedI", "I", "Tt"))
  
 # sim$addLogger(newCounter("VoluntaryA", "A", "T"))
#  sim$addLogger(newCounter("casesTracedA", "QA", "Tt"))
#  sim$addLogger(newCounter("traced", "T", "X"))
  
  return(sim$run(0:250))
}

# Number of runs
n_runs = 100

# Run the simulation 100 times and store results
results = replicate(n_runs, run(), simplify = FALSE)


I_values = matrix(0, nrow = 251, ncol = n_runs)
T_values = matrix(0, nrow = 251, ncol = n_runs)
Tt_values = matrix(0, nrow = 251, ncol = n_runs)
VoluntaryI_values = matrix(0, nrow = 251, ncol = n_runs)
QuarantinedE_values = matrix(0, nrow = 251, ncol = n_runs)
onset_values = matrix(0, nrow = 251, ncol = n_runs)
TracedI_values = matrix(0, nrow = 251, ncol = n_runs)



for (i in 1:n_runs) {
  res = results[[i]]
  I_values[, i] = res$I
  T_values[, i] = res$T
  Tt_values[, i] = res$Tt
  VoluntaryI_values[, i] = res$VoluntaryI
  QuarantinedE_values[, i] = res$QuarantinedE
  onset_values[, i] = res$onset
  TracedI_values[, i] = res$TracedI
}

# Calculate the mean for each time step across all runs
mean_I = rowMeans(I_values)
mean_T = rowMeans(T_values)
mean_Tt = rowMeans(Tt_values)
mean_VoluntaryI = rowMeans(VoluntaryI_values)
mean_QuarantinedE = rowMeans(QuarantinedE_values)
mean_onset = rowMeans(onset_values)
mean_TracedI = rowMeans(TracedI_values)

# Print the mean onset values for each time step
cat("Time\tMean Onset\n")
for (t in 1:length(mean_onset)) {
  cat(t, "\t", mean_I[t], "\n")
  cat(t, "\t", mean_T[t], "\n")
  cat(t, "\t", mean_Tt[t], "\n")
  cat(t, "\t", mean_VoluntaryI[t], "\n")
  cat(t, "\t", mean_QuarantinedE[t], "\n")
  cat(t, "\t", mean_onset[t], "\n")
  cat(t, "\t", mean_TracedI[t], "\n")
}


write.csv(data.frame(
  Time = 0:250,
  MeanI = mean_I,
  MeanT = mean_T,
  MeanTt = mean_Tt,
  MeanVoluntaryI = mean_VoluntaryI,
  MeanQuaranrinedE = mean_QuarantinedE,
  MeanOnset = mean_onset,
  MeanTracedI = mean_TracedI
), "mean_results.csv", row.names = FALSE)
