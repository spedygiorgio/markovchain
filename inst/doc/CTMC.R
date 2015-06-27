## ----ctmcInit, echo = TRUE, message=FALSE, warning=FALSE-----------------
library(markovchain)
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3,
                       1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, 
                 byrow = byRow, generator = gen, 
                 name = "Molecular Transition Model")      

## ----ctmcRandom0, echo = TRUE, message=FALSE, warning=FALSE--------------
statesDist <- c(0.8, 0.2)
rctmc(n = 3, ctmc = molecularCTMC, initDist = statesDist)

## ----ctmcRandom1, echo = TRUE, message=FALSE, warning=FALSE--------------
statesDist <- c(0.8, 0.2)
rctmc(n = Inf, ctmc = molecularCTMC, initDist = statesDist, T = 1)

