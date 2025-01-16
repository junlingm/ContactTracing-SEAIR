library(ABM)
if (packageVersion("ABM") < "0.4.2") {
  stop("need ABM version >= 0.4.2")
}

N=3000 # the population size
I0 = 20 # the initial number of infections

p = 0.2 # the contact tracing coverage probability
q = 0.3
theta = 2 # the contact tracing rate for diagnosed patients

beta = 0.6 # the contact rate
betaA = beta/3

sigma = 0.27

gammaI = 0.1 # the recovery rate
gammaA = 0.2 # the recovery rate
tauI = 0.15 
tauA = 0 

t = 0:250

form_pairs = function(time, agent, contact) { 
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
    state = if (i<=I0) "I" else if (i <= 2*I0) "A" else if (i <= 3*I0) "E" else "S"
    list(state, contacts=list())
  })
  m = newRandomMixing()
  sim$addContact(m)
  
  sim$addTransition("E" -> "I", sigma*(1-q))
  sim$addTransition("E" -> "A", sigma*q)
  sim$addTransition("A" -> "T", tauA)
  sim$addTransition("A" -> "R", gammaA)
  sim$addTransition("I"->"R", gammaI)
  sim$addTransition("I"->"T", tauI)
  sim$addTransition(
    "I" + "S" -> "I" + "E" ~ m, beta,
    changed_callback = form_pairs)
  sim$addTransition(
    "A" + "S" -> "A" + "E" ~ m, betaA,
    changed_callback = form_pairs)
  sim$addTransition( # voluntary
    "T"->"X", theta, 
    changed_callback = trace)
  sim$addTransition( # contact tracing
    "Tt"->"X", theta, 
    changed_callback = trace)
  sim$addTransition("QA" -> "Tt", tauA)
  sim$addTransition("QA" -> "R", gammaA)
  sim$addTransition("QE" -> "QA", sigma*q)
  sim$addTransition("QE" -> "Tt", sigma*(1-q))
  
  
  sim$addLogger(newCounter("I", "I"))
  sim$run(t)
}

# Number of runs
n_runs = 10

# Run the simulation 100 times and store results
results = replicate(n_runs, run(), simplify = FALSE)

# Initialize an empty matrix to store "onset" values for each run
I_values = matrix(0, nrow = 251, ncol = n_runs)

# Loop through each run and extract the "onset" values
for (i in 1:n_runs) {
  res = results[[i]]
  I_values[, i] = res$I
}

# Calculate the mean "onset" for each time step across all runs
mean_I = rowMeans(I_values)

# Calculate the quantiles (0.025 and 0.975) for each time step
quantile_I_025 = apply(I_values, 1, function(x) quantile(x, 0.025))
quantile_I_975 = apply(I_values, 1, function(x) quantile(x, 0.975))

# Optionally, save the results to a CSV file
sim.data <- rbind(
  data.frame(
    Time = t,
    value = mean_I,
    curve="mean"
  ),
  data.frame(
    Time = t,
    value = quantile_I_025,
    curve="2.5% quantile"
  ),
  data.frame(
    Time = t,
    value = quantile_I_975,
    curve = "97,.5% quantile"
  )
)
write.csv(sim.data, file="mean_I_quantile_results.csv", row.names = FALSE)




