
##
## Tutorial 2: Dynamic Network Modeling with STERGMs in EpiModel
## Day 1 | Network Models and HIV/STI with EpiModel | Harvard 2017"
##


# Set the seed to replicate the results
set.seed(0)

# Load the R package
library(EpiModel)

# Set up the network
net1 <- network.initialize(100, directed = FALSE)

# Input the dissolution model coefficients
coef.diss.1 <- dissolution_coefs(~offset(edges), 90)
coef.diss.1

# Fit the model with netest
fit1 <- netest(net1,
               formation = ~edges,
               target.stats = 20,
               coef.diss = coef.diss.1)

# Summary of the model fit
summary(fit1)

# netdx simulates dynamic networks given the model fit
sim1 <- netdx(fit1, nsteps = 1000, nsims = 10,
              keep.tedgelist = TRUE)

# Numerical comparison of diagnostics
sim1

# Plot the diagnostics for the formation model
plot(sim1, type = "formation")

# Plot the durational diagnostics for the dissolution model
plot(sim1, type = "duration")

# Plot the (unbiased) dissolutional diagnostics
plot(sim1, type = "dissolution")

# Examine the timed edgelist
sim1$tedgelist[[1]][1:5,]

# Analysis on the timed edgelist
tel <- sim1$tedgelist[[1]]
hist(tel$duration)
mean(tel$duration[tel$onset < 100])
sum(tel$terminus.censored == TRUE)
plot(tel$onset, tel$terminus)
table(c(tel$head,tel$tail))
hist(table(c(tel$head,tel$tail)))

# Increase the population size
net2 <- network.initialize(1000, directed = FALSE)
fit2 <- netest(net2,
   formation = ~edges,
   target.stats = 200,
   coef.diss = dissolution_coefs(~offset(edges), 90))
sim2 <- netdx(fit2, nsteps = 1000, nsims = 10, keep.tedgelist = TRUE)
plot(sim2, type = "formation")
plot(sim2, type = "duration")
plot(sim2, type = "dissolution")

# Two-race MSM model
n <- 500
net3 <- network.initialize(n, directed = FALSE)
net3 %v% "race" <- c(rep("B", n/2), rep("W", n/2))
net3

# Formation formula
form.formula.3 <- ~edges + nodematch("race") + degree(0) + concurrent
target.stats.3 <- c(0.9*n/2, (0.9*n/2)*(5/6), 0.36*n, 0.18*n)

# Dissolution formula
diss.formula.3 <- ~offset(edges) + offset(nodematch("race"))

# Look at the supported dissolution models
?dissolution_coefs

# Fit the model
fit3 <- netest(net3,
             formation = form.formula.3,
             target.stats = target.stats.3,
             coef.diss = dissolution_coefs(~offset(edges) + offset(nodematch("race")),
                                           c(200, 100)))

# netdx to simulate the dynamic network
sim3 <- netdx(fit3, nsteps = 1000, nsims = 10, keep.tedgelist = TRUE)

# Print the results
sim3

# Plot the formation model stats
plot(sim3, type = "formation")

# Plot the dissolution model stats (sorry!)
plot(sim3, type = "duration")

# An alternative approach
race <- net3 %v% "race"
tel3 <- sim3$tedgelist[[1]]
mean(tel3$duration[(race[tel3$tail] != race[tel3$head]) & tel3$onset < 100])
mean(tel3$duration[(race[tel3$tail] == race[tel3$head]) & tel3$onset < 100])

# The edges dissolution approximation
fit4 <- netest(net3,
       formation = form.formula.3,
       target.stats = target.stats.3,
       coef.diss = dissolution_coefs(~offset(edges)+offset(nodematch("race")),
                                     c(20, 10)))

# Run the diagnostics and examine the output
sim4 <- netdx(fit4, nsteps = 1000, nsims = 10, keep.tedgelist = TRUE)

sim4
plot(sim4, type = "formation")

# Use the full STERGM estimation method
fit5 <- netest(net3,
         formation = form.formula.3,
         target.stats = target.stats.3,
         coef.diss = dissolution_coefs(~offset(edges) + offset(nodematch("race")),
                                       c(20, 10)),
         edapprox = FALSE)

# Compare the performance here
sim5 <- netdx(fit5, nsteps = 1000, nsims = 10, keep.tedgelist = TRUE)

plot(sim5, type = "formation")

race <- net3 %v% "race"
tel5 <- sim5$tedgelist[[1]]
mean(tel5$duration[(race[tel5$tail] != race[tel5$head]) & tel5$onset < 100])
mean(tel5$duration[(race[tel5$tail] == race[tel5$head]) & tel5$onset < 100])
