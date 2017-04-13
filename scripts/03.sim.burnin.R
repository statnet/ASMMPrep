
library("EpiModelHPC")
library("EpiModelAdol")

rm(list = ls())

load("scenarios/adol/est/nwstats.rda")
load("scenarios/adol/est/fit.rda")



param <- param.adol(nwstats = st, prep.elig.model = "none")


init <- init.adol(nwstats = st)
control <- control.adol(nsteps = 1040, nsims = 1, prevfull = TRUE, save.transmat = TRUE,
                        save.other = "attr", save.network = TRUE,
                        verbose = TRUE, verbose.int = 1)

sim.counts <- netsim(est, param, init, control)


save(sim.counts, file = "scenarios/adol/est/sim.counts.rda")




