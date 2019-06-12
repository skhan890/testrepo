
##
## Tutorial 3: Extending EpiModel, A Different Demographic Approach
## Day 3 | Network Models and HIV/STI with EpiModel | Harvard 2017"
##

library("EpiModel")

aging <- function(dat, at) {

  ## Attributes
  if (at == 2) {
    n <- sum(dat$attr$active == 1)
    dat$attr$age <- sample(18:49, n, replace = TRUE)
  } else {
    dat$attr$age <- dat$attr$age + 1/12
  }

  ## Summary statistics
  if (at == 2) {
    dat$epi$meanAge <- rep(mean(dat$attr$age, na.rm = TRUE), 2)
  } else {
    dat$epi$meanAge[at] <- mean(dat$attr$age, na.rm = TRUE)
  }

  return(dat)
}

ages <- 18:49
death.rates <- 1/(70*12 - ages*12)
par(mar = c(3.2, 3.2, 1, 1), mgp = c(2, 1, 0))
plot(ages, death.rates, pch = 20, xlab = "age", ylab = "Death Risk")

dfunc <- function(dat, at) {

  # Parameters
  idsElig <- which(dat$attr$active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  # Processes
  if (nElig > 0) {
    ages <- dat$attr$age[idsElig]
    life.expt <- dat$param$life.expt
    death.rates <- pmin(1, 1/(life.expt*12 - ages*12))
    vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    # Update nodal attributes on attr and networkDynamic object
    if (nDeaths > 0) {
      dat$attr$active[idsDeaths] <- 0
      dat$attr$exitTime[idsDeaths] <- at
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDeaths, deactivate.edges = TRUE)
    }
  }

  # Summary statistics
  if (at == 2) {
    dat$epi$d.flow <- c(0, nDeaths)
  } else {
    dat$epi$d.flow[at] <- nDeaths
  }

  return(dat)
}

bfunc <- function(dat, at) {

  # Variables
  growth.rate <- dat$param$growth.rate
  exptPopSize <- dat$epi$num[1] * (1 + growth.rate * at)
  n <- network.size(dat$nw)
  tea.status <- dat$control$tea.status

  numNeeded <- exptPopSize - sum(dat$attr$active == 1)
  if (numNeeded > 0) {
    nBirths <- rpois(1, numNeeded)
  } else {
    nBirths <- 0
  }
  if (nBirths > 0) {
    dat$nw <- add.vertices(dat$nw, nv = nBirths)
    newNodes <- (n + 1):(n + nBirths)
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
  }

  # Update attributes
  if (nBirths > 0) {
    dat$attr$active <- c(dat$attr$active, rep(1, nBirths))
    dat$attr$status <- c(dat$attr$status, rep("s", nBirths))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nBirths))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nBirths))
    dat$attr$age <- c(dat$attr$age, rep(18, nBirths))
    if (tea.status == TRUE) {
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                          value = 0, onset = at,
                                          terminus = Inf, v = newNodes)
    }
  }

  # Summary statistics
  if (at == 2) {
    dat$epi$b.flow <- c(0, nBirths)
  } else {
    dat$epi$b.flow[at] <- nBirths
  }

  return(dat)
}

nw <- network.initialize(500, directed = FALSE)
est <- netest(nw, formation = ~edges, target.stats = 150,
              coef.diss = dissolution_coefs(~offset(edges), 60, mean(death.rates)))

param <- param.net(inf.prob = 0.15, growth.rate = 0.00083, life.expt = 70)
init <- init.net(i.num = 50)

control <- control.net(type = "SI", nsims = 5, nsteps = 250,
                       deaths.FUN = dfunc, births.FUN = bfunc, aging.FUN = aging,
                       depend = TRUE, save.network = FALSE)

mod <- netsim(est, param, init, control)

mod

par(mfrow = c(1,2))
plot(mod, main = "State Prevalences")
plot(mod, main = "State Sizes", sim.lines = TRUE,
     qnts = FALSE, mean.smooth = FALSE)

par(mfrow = c(1, 2))
plot(mod, y = "num", main = "Population Size", ylim = c(0, 1000))
plot(mod, y = "meanAge", main = "Mean Age", ylim = c(18, 70))

par(mfrow = c(1, 2))
plot(mod, y = "d.flow", mean.smooth = TRUE, qnts = 1, main = "Deaths")
plot(mod, y = "b.flow", mean.smooth = TRUE, qnts = 1, main = "Births")

