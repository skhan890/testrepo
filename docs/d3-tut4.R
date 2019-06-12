
##
## Tutorial 4: Extending EpiModel, Migration
## Day 3 | Network Models and HIV/STI with EpiModel | Harvard 2017"
##

## ----setup, message = FALSE----------------------------------------------
library("EpiModel")

## ----Ex2omigFunc---------------------------------------------------------
outmigrate <- function(dat, at) {

  # Variables
  active <- dat$attr$active
  status <- dat$attr$status
  rates <- dat$param$omig.rates

  # Process
  idsElig <- which(active == 1)
  nElig <- length(idsElig)

  nMig <- 0
  if (nElig > 0) {
    ratesElig <- rates[as.numeric(status == "i") + 1]
    vecMig <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecMig) > 0) {
      idsMig <- idsElig[vecMig]
      nMig <- length(idsMig)

      dat$attr$active[idsMig] <- 0
      dat$attr$exitTime[idsMig] <- at
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsMig, deactivate.edges = TRUE)
    }
  }

  # Summary statistics
  if (at == 2) {
    dat$epi$omig.flow <- c(0, nMig)
  } else {
    dat$epi$omig.flow[at] <- nMig
  }

  return(dat)
}

## ----Ex2imigFunc---------------------------------------------------------
inmigrate <- function(dat, at) {

  # Variables
  nw <- dat$nw
  n <- network.size(nw)
  active <- dat$attr$active

  exp.inmig <- dat$param$exp.inmig
  risk <- dat$param$imig.risk
  tea.status <- dat$control$tea.status

  # Add Nodes
  nMig <- 0
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  if (nElig > 0) {
    nMig <- rpois(1, exp.inmig)
    if (nMig > 0) {
      dat$nw <- add.vertices(dat$nw, nv = nMig)
      newNodes <- (n + 1):(n + nMig)
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf,
                                  v = newNodes)
    }
  }

  # Update attributes
  if (nMig > 0) {
    dat$attr$active <- c(dat$attr$active, rep(1, nMig))
    newStatus <- rbinom(nMig, 1, risk)
    newStatus <- ifelse(newStatus == 1, "i", "s")
    if (tea.status == TRUE) {
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus", value = newStatus,
                                          onset = at, terminus = Inf, v = newNodes)
    }
    dat$attr$status <- c(dat$attr$status, newStatus)
    infTime <- ifelse(newStatus == "i", at, NA)
    dat$attr$infTime <- c(dat$attr$infTime, infTime)
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nMig))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nMig))
  }

  # Summary statistics
  if (at == 2) {
    dat$epi$imig.flow <- c(0, nMig)
  } else {
    dat$epi$imig.flow[at] <- nMig
  }

  return(dat)
}

## ----Ex2nwParam----------------------------------------------------------
nw <- network.initialize(500, directed = FALSE)
est <- netest(nw, formation = ~edges, target.stats = 150,
              coef.diss = dissolution_coefs(~offset(edges), 60, 0.0046))

## ----Ex2epiParam---------------------------------------------------------
param <- param.net(inf.prob = 0.2, rec.rate = 0.02,
                   omig.rates = c(0.005, 0.001), exp.inmig = 3, imig.risk = 0.1)
init <- init.net(i.num = 50)

## ----Ex2epiControl-------------------------------------------------------
control <- control.net(type = "SIS", nsims = 5, nsteps = 500,
                       omig.FUN = outmigrate, imig.FUN = inmigrate,
                       depend = TRUE, save.network = FALSE)

## ----Ex2runNetsim, cache = TRUE, results = "hide"------------------------
mod <- netsim(est, param, init, control)

## ----Ex2print------------------------------------------------------------
mod
plot(mod)

## ----Ex2plot-------------------------------------------------------------
par(mfrow = c(1, 2))
plot(mod, y = "omig.flow", main = "Out Migration")
plot(mod, y = "imig.flow", main = "In Migration")

