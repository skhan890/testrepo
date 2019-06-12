
##
## Tutorial 4: Interventions
## Day 2 | Network Models and HIV/STI with EpiModel | Harvard 2017"
##

# Load EpiModel
library(EpiModel)

# Set up the network and estimate a simple random graph model
nw <- network.initialize(n = 100, directed = FALSE)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Parameterizing a simple intervention with two parameters: an intervention
# efficacy of 96% that starts on time step 25
param <- param.net(inf.prob = 0.5, inter.eff = 0.96, inter.start = 25)
init <- init.net(i.num = 5)
control <- control.net(type = "SI", nsteps = 100, nsims = 10)
sim <- netsim(est, param, init, control)
plot(sim)

# Another example for an SIS model
param <- param.net(inf.prob = 0.5, inter.eff = 0.8, inter.start = 100,
                   rec.rate = 0.07)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 250, nsims = 10)
sim <- netsim(est, param, init, control)
plot(sim)

