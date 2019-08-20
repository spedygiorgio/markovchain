#FITTING MARKOV CHAINS

#esempi di check fitting markov chains x tests
ciao <- c("a", "a", "b", "b", "a", NA, "b", "a", "b", "a", "a")
simpleMcCiaoFit <- markovchainFit(ciao)
data(rain)
checksAlofiRawTransitions <- createSequenceMatrix(rain$rain)

#check by matrix

data(holson)
myHolson <- as.matrix(holson[,-1])
rownames(myHolson) <- holson$id
checkmarkovchainFitList <- markovchainListFit(data=myHolson)
#devtools::use_data(simpleMcCiaoFit,checksAlofiRawTransitions,checkmarkovchainFitList,internal = TRUE,overwrite = TRUE)

#STEADY STATE ANALYSIS

mc1 <- as(matrix(c(0.5, 0.5, 0, 0,
                   0.3, 0.3, 0.4, 0,
                   0, 0.5, 0.5, 0,
                   0, 0, 0.5, 0.5),
                 byrow = TRUE, nrow = 4), 
          "markovchain")


mc2 <- as(matrix(c(0, 0.5, 0.50, 0, 0, 0,
                   0, 0, 0.99, 0, 0, 0.01,
                   0, 0.8, 0, 0.2, 0, 0,
                   0, 0, 0, 0, 1, 0,
                   0, 0, 0, 1, 0, 0,
                   0, 0, 0, 0, 0, 1), nrow = 6, byrow = TRUE),
          "markovchain")


mc3 <- as(matrix(c(0.4, 0.6, 0, 0, 0,
                   0.5, 0.5, 0, 0, 0,
                   0, 0, 0.3, 0.7, 0,
                   0, 0, 0.5, 0.4, 0.1,
                   0, 0, 0, 0.8, 0.2), nrow = 5, ncol = 5, byrow = TRUE),
          "markovchain")


statesNames <- letters[1:5]
mc4 <- new("markovchain", 
           transitionMatrix = matrix(c(0, 1/3, 0, 2/3, 0,
                                       1/2, 0, 0, 0, 1/2,
                                       0, 0, 1/2, 1/2, 0,
                                       0, 0, 1/2, 1/2, 0,
                                       0, 0, 0, 0, 1), 
                                     nrow = 5, ncol = 5, byrow = T),
           name = "Mathematica MC", 
           states = statesNames)


mc5 <- matrix(c(0.5, 0.5, 0, 0,
                0.5, 0.5, 0, 0,
                1/3, 1/6, 1/6, 1/3,
                0, 0, 0, 1), 
              nrow = 4,ncol = 4, byrow = T)
mc5 <- as(mc5, "markovchain")

mcRain <- markovchainFit(data = rain$rain)$estimate

generateSteadyStates <- function(mc) {
  list(object = mc,
       steadyStates = steadyStates(mc)
  )
}

MCs <- list(mc1, mc2, mc3, mc4, mc5, mcRain)
tMCs <- lapply(MCs, t)
knownSteadyStatesMCs <- append(MCs, tMCs)

M <- matlab::zeros(5, 5)
M[1,1] <- M[5,5] <- 1
M[2,1] <- M[2,3] <- 1/2
M[3,2] <- M[3,4] <- 1/2
M[4,2] <- M[4,5] <- 1/2

mcHitting <- new("markovchain", transitionMatrix = M)
#SAVING

usethis::use_data(simpleMcCiaoFit, checksAlofiRawTransitions, checkmarkovchainFitList, knownSteadyStatesMCs,
                  mcHitting, internal = TRUE, overwrite = TRUE)