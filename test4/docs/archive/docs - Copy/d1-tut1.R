
##
## Tutorial 1: Exponential Random Graph Models (ERGMs) using statnet
## Day 1 | Network Models and HIV/STI with EpiModel | Harvard 2017"
##


# Load necessary R packages
library(statnet)
library(ergm)
library(sna)
library(coda)

# latest versions: ergm 3.7.1 and network 1.13.0 (as of 6/1/2017)
sessionInfo()

# To reproduce exact results
set.seed(0)

# tells us the datasets in our packages
data(package = "ergm")


# flomarriage model -------------------------------------------------------

# loads flomarriage and flobusiness data
data(florentine)

# Let's look at the flomarriage network properties
flomarriage

# Setup a 2 panel plot (for later)
par(mfrow = c(1, 2))

# Plot the flomarriage network
plot(flomarriage, main = "Florentine Marriage", cex.main = 0.8)

# Look at the $g(y)$ statistic for this model
summary(flomarriage~edges)

# Estimate the model
flomodel.01 <- ergm(flomarriage~edges)

# The fitted model object
summary(flomodel.01)

# Look at the g(y) stats for this model
summary(flomarriage ~ edges + triangle)
flomodel.02 <- ergm(flomarriage ~ edges + triangle)
summary(flomodel.02)

# this has the class ergm
class(flomodel.02)

# the ERGM object contains lots of components.
names(flomodel.02)

# you can extract/inspect individual components
flomodel.02$coef

# %v% references vertex attributes
wealth <- flomarriage %v% "wealth"
wealth

# summarize the distribution of wealth
summary(wealth)

# network plot with vertex size proportional to wealth
plot(flomarriage, vertex.cex = wealth/25,
     main = "Florentine marriage by wealth", cex.main = 0.8)

# observed statistics for the model
summary(flomarriage ~ edges + nodecov("wealth"))

flomodel.03 <- ergm(flomarriage ~ edges + nodecov("wealth"))
summary(flomodel.03)


# Add Health example ------------------------------------------------------

data(faux.mesa.high)
mesa <- faux.mesa.high

mesa

# Back to 1-panel plots
par(mfrow = c(1, 1))
plot(mesa, vertex.col = "Grade")
legend("bottomleft", fill = 7:12,
       legend = paste("Grade", 7:12), cex = 0.75)

fauxmodel.01 <- ergm(mesa ~ edges + nodematch("Grade", diff = TRUE) +
                       nodematch("Race", diff = TRUE))
summary(fauxmodel.01)

# Frequencies of race
table(mesa %v% "Race")
mixingmatrix(mesa, "Race")

summary(mesa ~edges + nodematch("Grade", diff = TRUE) +
          nodematch("Race", diff = TRUE))

help("ergm-terms")

summary(flobusiness ~ edges + degree(1))
fit <- ergm(flobusiness ~ edges + degree(1))
mcmc.diagnostics(fit)

fit <- ergm(flobusiness ~ edges + degree(1),
            control = control.ergm(MCMC.interval = 1))

data("faux.magnolia.high")
magnolia <- faux.magnolia.high
plot(magnolia, vertex.cex = 0.5)
summary(magnolia~edges+triangle)

fit <- ergm(magnolia ~ edges + triangle)

fit <- ergm(magnolia ~ edges + triangle,
            control = control.ergm(MCMLE.maxit = 2))
mcmc.diagnostics(fit)

fit <- ergm(magnolia ~ edges + gwesp(0.25, fixed = TRUE) +
              nodematch("Grade") + nodematch("Race") + nodematch("Sex"),
            control = control.ergm(MCMC.samplesize = 50000, MCMC.interval = 1000),
            verbose = TRUE)

mcmc.diagnostics(fit)

flomodel.03.sim <- simulate(flomodel.03, nsim = 10)
class(flomodel.03.sim)
summary(flomodel.03.sim)
length(flomodel.03.sim)
flomodel.03.sim[[1]]
plot(flomodel.03.sim[[1]], label = flomodel.03.sim[[1]] %v% "vertex.names")

flomodel.03.gof <- gof(flomodel.03 ~ degree + esp + distance)
flomodel.03.gof
plot(flomodel.03.gof)

mesamodel.02 <- ergm(mesa ~ edges)
mesamodel.02.gof <- gof(mesamodel.02 ~ degree + esp + distance)
plot(mesamodel.02.gof)


# egocentric network data -------------------------------------------------

ego.net <- network.initialize(500, directed = FALSE)
ego.net %v% "sex" <- c(rep(0, 250), rep(1, 250))

# node distn
ego.deg <- c(180, 245, 60, 15)

# adjusted tie distn
ego.mixmat <- matrix(c(164, 44, 26, 176)/2,
                     nrow = 2, byrow = TRUE)

ego.edges <- sum(ego.mixmat)
ego.sexmatch <- ego.mixmat[1, 1] + ego.mixmat[2, 2]

ego.target.stats <- c(ego.edges, ego.sexmatch)
ego.target.stats

ego.fit <- ergm(ego.net ~ edges + nodematch("sex"),
                target.stats = ego.target.stats)

summary(ego.fit)

ego.sim1 <- simulate(ego.fit)
plot(ego.sim1, vertex.cex = 0.65, vertex.col = "sex")

rbind(sim = summary(ego.sim1 ~ degree(c(0:3))), obs = ego.deg)
mixingmatrix(ego.sim1, "sex")
ego.mixmat

ego.sim100 <- simulate(ego.fit, nsim = 100)
ego.sim100

sim.stats <- attr(ego.sim100, "stats")
rbind(sim = colMeans(sim.stats), obs = ego.target.stats)

matplot(1:nrow(sim.stats), sim.stats,
        pch = c("e", "m", "0", "+"), cex = 0.65,
        main = "100 simulations from ego.fit model",
        sub = "(default settings)",
        xlab = "Replicate", ylab = "frequency")
abline(h = ego.target.stats, col = c(1:4))

ego.sim100 <- simulate(ego.fit, nsim = 100,
                       control = control.simulate.ergm(MCMC.interval = 10000))
sim.stats <- attr(ego.sim100, "stats")
matplot(1:nrow(sim.stats), sim.stats,
        pch = c("e", "m"), cex = 0.65,
        main = "100 simulations from ego.fit model",
        sub = "(MCMC.interval=10000)",
        xlab = "Replicate", ylab = "frequency")
abline(h = ego.target.stats, col = 1:2)

sim.fulldeg <- summary(ego.sim100 ~ degree(0:10))
colnames(sim.fulldeg) <- paste("deg", 0:10, sep = "")
sim.fulldeg[1:5, ]

sim.deg <- cbind(sim.fulldeg[, 1:3], apply(sim.fulldeg[, 4:11], 1, sum))
colnames(sim.deg) <- c(colnames(sim.fulldeg)[1:3], "degree3+")
rbind(sim = colMeans(sim.deg), obs = ego.deg)

matplot(1:nrow(sim.deg), sim.deg, pch = as.character(0:3), cex = 0.5,
        main = "Comparing ego.sims to non-targeted degree frequencies",
        sub = "(only total edges targeted)",
        xlab = "Replicate", ylab = "Frequencies")
abline(h = c(180, 245, 60, 15), col = c(1:4))

ego.isolates <- ego.deg[1]
ego.target.stats <- c(ego.edges, ego.sexmatch, ego.isolates)
ego.fit <- ergm(ego.net ~ edges + nodematch("sex") + degree(0),
                target.stats = ego.target.stats)
summary(ego.fit)

ego.sim100 <- simulate(ego.fit, nsim = 100,
                       control = control.simulate.ergm(MCMC.interval = 10000))
sim.stats <- attr(ego.sim100, "stats")
rbind(sim = colMeans(sim.stats), obs = ego.target.stats)

sim.fulldeg <- summary(ego.sim100 ~ degree(0:10))
sim.deg <- cbind(sim.fulldeg[, 1:3], apply(sim.fulldeg[ ,4:11], 1, sum))
colnames(sim.deg) <- c(colnames(sim.fulldeg)[1:3], "degree3+")
rbind(sim = colMeans(sim.deg), obs = ego.deg)

matplot(1:nrow(sim.deg), sim.deg, pch = as.character(0:3), cex = 0.5,
        main = "Comparing ego.sims to non-targeted degree frequencies",
        sub = "(only 0, 2+ and total edges targeted)",
        xlab = "Replicate", ylab = "Frequencies")
abline(h = c(180, 245, 60, 15), col = c(1:4))

