
##
## Lab: HIV/STI Transmission Models with EpiModelHIV
## Day 5 | Network Models and HIV/STI with EpiModel | Harvard 2017"
##

# install.packages("EpiModel")
# install devtools if necessary, install.packages("devtools")
devtools::install_github("statnet/tergmLite")
devtools::install_github("statnet/EpiModelHIV", ref = "prep-sti")

library("EpiModelHIV")

data(st)
st

data(est)
est

param <- param_msm(st,
                  race.method = 1,
                  last.neg.test.B.int = 301,
                  last.neg.test.W.int = 315,
                  mean.test.B.int = 301,
                  mean.test.W.int = 315,
                  testing.pattern = "memoryless",
                  test.window.int = 21,

                  tt.traj.B.prob = c(0.077, 0.000, 0.356, 0.567),
                  tt.traj.W.prob = c(0.052, 0.000, 0.331, 0.617),

                  tx.init.B.prob = 0.092,
                  tx.init.W.prob = 0.127,
                  tx.halt.B.prob = 0.0102,
                  tx.halt.W.prob = 0.0071,
                  tx.reinit.B.prob = 0.00066,
                  tx.reinit.W.prob = 0.00291,

                  max.time.off.tx.full.int = 520 * 7,
                  max.time.on.tx.part.int = 52 * 15 * 7,
                  max.time.off.tx.part.int = 520 * 7,
                  vl.acute.rise.int = 45,
                  vl.acute.peak = 6.886,
                  vl.acute.fall.int = 45,
                  vl.set.point = 4.5,
                  vl.aids.onset.int = 520 * 7,
                  vl.aids.int = 52 * 2 * 7,
                  vl.fatal = 7,
                  vl.full.supp = 1.5,
                  vl.part.supp = 3.5,
                  full.supp.down.slope = 0.25,
                  full.supp.up.slope = 0.25,
                  part.supp.down.slope = 0.25,
                  part.supp.up.slope = 0.25,

                  b.B.rate = 1e-3 / 7,
                  b.W.rate = 1e-3 / 7,
                  birth.age = 18,
                  b.method = "fixed",

                  URAI.prob = 0.0082 * 1.09,
                  UIAI.prob = 0.0031 * 1.09,
                  acute.rr = 6,
                  circ.rr = 0.4,
                  condom.rr = 0.295,

                  disc.outset.main.B.prob = 0.685,
                  disc.outset.main.W.prob = 0.889,
                  disc.at.diag.main.B.prob = 1,
                  disc.at.diag.main.W.prob = 1,
                  disc.post.diag.main.B.prob = 0,
                  disc.post.diag.main.W.prob = 0,
                  disc.outset.pers.B.prob = 0.527,
                  disc.outset.pers.W.prob = 0.828,
                  disc.at.diag.pers.B.prob = 1,
                  disc.at.diag.pers.W.prob = 1,
                  disc.post.diag.pers.B.prob = 0,
                  disc.post.diag.pers.W.prob = 0,
                  disc.inst.B.prob = 0.445,
                  disc.inst.W.prob = 0.691,

                  circ.B.prob = 0.874,
                  circ.W.prob = 0.918,

                  ccr5.B.prob = c(0, 0.034),
                  ccr5.W.prob = c(0.021, 0.176),
                  ccr5.heteroz.rr = 0.3,

                  num.inst.ai.classes = 1,
                  base.ai.main.BB.rate = 0.17,
                  base.ai.main.BW.rate = 0.26,
                  base.ai.main.WW.rate = 0.23,
                  base.ai.pers.BB.rate = 0.11,
                  base.ai.pers.BW.rate = 0.16,
                  base.ai.pers.WW.rate = 0.14,
                  ai.scale = 1.15,

                  cond.main.BB.prob = 0.38,
                  cond.main.BW.prob = 0.10,
                  cond.main.WW.prob = 0.15,
                  cond.pers.always.prob = 0.216,
                  cond.pers.BB.prob = 0.26,
                  cond.pers.BW.prob = 0.26,
                  cond.pers.WW.prob = 0.26,
                  cond.inst.always.prob = 0.326,
                  cond.inst.BB.prob = 0.27,
                  cond.inst.BW.prob = 0.27,
                  cond.inst.WW.prob = 0.27,
                  cond.always.prob.corr = 0.5,
                  cond.rr.BB = 1,
                  cond.rr.BW = 1,
                  cond.rr.WW = 1,
                  cond.diag.main.beta = -0.67,
                  cond.discl.main.beta = -0.85,
                  cond.diag.pers.beta = -0.67,
                  cond.discl.pers.beta = -0.85,
                  cond.diag.inst.beta = -0.67,
                  cond.discl.inst.beta = -0.85,

                  vv.iev.BB.prob = 0.42,
                  vv.iev.BW.prob = 0.56,
                  vv.iev.WW.prob = 0.49,

                  prep.start = Inf,
                  prep.elig.model = "base",
                  prep.class.prob = c(0.211, 0.07, 0.1, 0.619),
                  prep.class.hr = c(1, 0.69, 0.19, 0.05),
                  prep.coverage = 0,
                  prep.cov.method = "curr",
                  prep.cov.rate = 1,
                  prep.tst.int = 90,
                  prep.risk.int = 182,
                  prep.risk.reassess = TRUE,

                  rcomp.prob = 0,
                  rcomp.adh.groups = 0:3,
                  rcomp.main.only = FALSE,
                  rcomp.discl.only = FALSE,

                  rgc.tprob = 0.357698,
                  ugc.tprob = 0.248095,
                  rct.tprob = 0.321597,
                  uct.tprob = 0.212965,

                  rgc.sympt.prob = 0.076975,
                  ugc.sympt.prob = 0.824368,
                  rct.sympt.prob = 0.103517,
                  uct.sympt.prob = 0.885045,

                  rgc.asympt.int = 35.11851 * 7,
                  ugc.asympt.int = 35.11851 * 7,
                  gc.tx.int = 2 * 7,
                  gc.ntx.int = NA,

                  rct.asympt.int = 44.24538 * 7,
                  uct.asympt.int = 44.24538 * 7,
                  ct.tx.int = 2 * 7,
                  ct.ntx.int = NA,

                  gc.prob.cease = 0,
                  ct.prob.cease = 0,

                  gc.sympt.prob.tx = 0.90,
                  ct.sympt.prob.tx = 0.85,
                  gc.asympt.prob.tx = 0,
                  ct.asympt.prob.tx = 0,

                  prep.sti.screen.int = 182,
                  prep.sti.prob.tx = 1,
                  prep.continue.stand.tx = TRUE,

                  sti.cond.rr = 0.3,

                  hiv.rgc.rr = 2.780673,
                  hiv.ugc.rr = 1.732363,
                  hiv.rct.rr = 2.780673,
                  hiv.uct.rr = 1.732363,
                  hiv.dual.rr = 0.2)

init <- init_msm(nwstats = st,
                 prev.B = 0.253,
                 prev.W = 0.253,
                 prev.ugc = 0.005,
                 prev.rgc = 0.005,
                 prev.uct = 0.013,
                 prev.rct = 0.013)

control <- control_msm(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        initialize.FUN = initialize_msm,
                        aging.FUN = aging_msm,
                        deaths.FUN = deaths_msm,
                        births.FUN = births_msm,
                        test.FUN = test_msm,
                        tx.FUN = tx_msm,
                        prep.FUN = prep_msm,
                        progress.FUN = progress_msm,
                        vl.FUN = vl_msm,
                        aiclass.FUN = NULL,
                        roleclass.FUN = NULL,
                        resim_nets.FUN = simnet_msm,
                        disclose.FUN = disclose_msm,
                        acts.FUN = acts_msm,
                        condoms.FUN = condoms_msm,
                        riskhist.FUN = riskhist_msm,
                        position.FUN = position_msm,
                        trans.FUN = trans_msm,
                        stitrans.FUN = sti_trans,
                        stirecov.FUN = sti_recov,
                        stitx.FUN = sti_tx,
                        prev.FUN = prevalence_msm,
                        verbose.FUN = verbose_msm,
                        save.nwstats = FALSE,
                        verbose = TRUE,
                        verbose.int = 1)

sim <- netsim(est, param, init, control)

dat <- initialize_msm(est, param, init, control, s = 1)

# debugonce(progress_msm)
for (at in 2:100) {
  dat <- aging_msm(dat, at)       ## <1 ms
  dat <- deaths_msm(dat, at)      ## 4 ms
  dat <- births_msm(dat, at)      ## 6 ms
  dat <- test_msm(dat, at)        ## 2 ms
  dat <- tx_msm(dat, at)          ## 3 ms
  dat <- prep_msm(dat, at)        ## 2 ms
  dat <- progress_msm(dat, at)    ## 2 ms
  dat <- vl_msm(dat, at)          ## 3 ms
  dat <- simnet_msm(dat, at)      ## 53 ms
  dat <- disclose_msm(dat, at)    ## 1 ms
  dat <- acts_msm(dat, at)        ## 1 ms
  dat <- condoms_msm(dat, at)     ## 2 ms
  dat <- riskhist_msm(dat, at)    ## 4 ms
  dat <- position_msm(dat, at)    ## 1 ms
  dat <- trans_msm(dat, at)       ## 1 ms
  dat <- sti_trans(dat, at)       ## 4 ms
  dat <- sti_recov(dat, at)       ## 3 ms
  dat <- sti_tx(dat, at)          ## 2 ms
  dat <- prevalence_msm(dat, at)  ## 1 ms
  cat(at, ".", sep = "")
}

data(st)
param <- param_msm(st,
                  prep.start = 100, # change this > 26 but < nsteps
                  prep.elig.model = "base",
                  prep.class.prob = c(0.211, 0.07, 0.1, 0.619),
                  prep.class.hr = c(1, 0.69, 0.19, 0.05),
                  prep.coverage = 0, # change this between 0 and 1
                  prep.tst.int = 90,
                  prep.risk.int = 182,
                  prep.risk.reassess = TRUE,
                  rcomp.prob = 0,
                  rcomp.adh.groups = 0:3,
                  rcomp.main.only = FALSE,
                  rcomp.discl.only = FALSE,
                  prep.sti.screen.int = 182,
                  prep.sti.prob.tx = 1,
                  prep.continue.stand.tx = TRUE)
init <- init_msm(nwstats = st)
control <- control_msm(simno = 1,
                       nsteps = 200, # change this to be longer if necessary
                       nsims = 1, # bump this up
                       ncores = 1,
                       verbose = TRUE)
data(est)
sim <- netsim(est, param, init, control)
sim

# Here are some plots you might be interested in
plot(sim, y = "ir100")
plot(sim, y = "i.prev")
plot(sim, y = "ir100.gc")
plot(sim, y = "ir100.ct")

# Extract the data and conduct some data analysis
df <- as.data.frame(sim)
names(df)

# Cumulative incidence of NG per 100 person-years at risk
sum(df$ir100.gc, na.rm = TRUE)

# Average incidence rate
mean(df$ir100.gc, na.rm = TRUE)

