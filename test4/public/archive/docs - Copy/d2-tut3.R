
##
## Tutorial 3: Time-Varying Biology & Behavior
## Day 2 | Network Models and HIV/STI with EpiModel | Harvard 2017"
##

# Load EpiModel
library(EpiModel)


# Generic Model -----------------------------------------------------------

# Define a time-varying probability of transmission as a function of how many
# time steps the infector has been infected. This code implies a 50% per-act
# transmission probability for the first 10 time steps and then a 5% for the
# remainder of a person's infection.
probs <- c(0.5, 0.05)
durs <- c(10, 1)
inf.probs <- rep(probs, durs)
inf.probs

# Act rates may also be time-varying. Here we chose a Random Poisson function to
# draw variabile number of acts for the first 10 time steps of an infector's
# infection with a mean of 3 acts, and then the deterministic mean for the remainder
# of the infection.
act.rates <- c(rpois(10, lambda = 3), 3)
act.rates

# A two panel plot to examine the results
par(mfrow = c(1, 2))
plot(inf.probs, type = "S", lwd = 2, ylim = 0:1)
plot(act.rates, type = "S", lwd = 2)
abline(h = 3, lty = 2)

# Fit a simple random graph model and simulate an SIS disease with these
# varying parameters
nw <- network.initialize(n = 100, directed = FALSE)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
param <- param.net(inf.prob = inf.probs, act.rate = act.rates, rec.rate = 0.01)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 100, nsims = 1)
sim <- netsim(est, param, init, control)

# The transmission matrix allows for an intuitive demonstration of what has occured
tm <- get_transmat(sim)
head(tm, 10)
mean(tm$infDur <= 10)



# Stage-Specific HIV Model ------------------------------------------------

# Transmission *rates* per month by acute, chronic, early AIDS, and late AIDS
# stages. These stages will last the number of months indicted in durs, at which
# point a person with late-stage AIDS will be assumed to die from disease. In this
# parameterization, it is assumed that there is one "act" per time step (in reality,
# this is a transmission rate per partnership-month, and the number of acts per
# month are/were unknown)
probs <- c(0.2055, 0.0088, 0.0614, 0)
durs <- c(3, 100, 9, 10)
inf.probs <- rep(probs, durs)
inf.probs

# Note in the example above that the late-stage AIDS probability was 0. That is
# was estimated as a behavioral, not biological function: persons in this stage
# effectively had no sexual activity. Here's an alternative parameterization that
# does the same thing, but makes the bio-behavioral interaction clearer.
probs <- c(0.2055, 0.0088, 0.0614, 0.1)
acts <- c(1, 1, 1, 0)
durs <- c(3, 100, 9, 10)
inf.probs <- rep(probs, durs)
act.rates <- rep(acts, durs)

# A plot of the results of that interaction
par(mfrow = c(1,1))
plot(inf.probs, type = "S", ylim = c(0, 1), lwd = 2)
lines(act.rates, type = "S", col = 2, lwd = 2)
legend(1, 0.8, legend = c("inf.probs", "act.rates"),
       lwd = 3, col = 1:2, bty = "n")

