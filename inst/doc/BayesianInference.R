## ----loadAndDoExample, echo = TRUE, message=FALSE, warning=FALSE---------
library(markovchain)
weatherStates <- c("sunny", "cloudy", "rain")
byRow <- TRUE
weatherMatrix <- matrix(data = c(0.7, 0.2, 0.1, 
                                 0.3, 0.4, 0.3, 
                                 0.2, 0.4, 0.4), 
                        byrow = byRow, nrow = 3, 
                        dimnames = list(weatherStates, weatherStates))
mcWeather <- new("markovchain", states = weatherStates, 
                 byrow = byRow, transitionMatrix = weatherMatrix, 
                 name = "Weather")      
weathersOfDays <- rmarkovchain(n = 365, object = mcWeather, t0 = "sunny")

## ----MAPFit--------------------------------------------------------------
hyperMatrix<-matrix(c(1, 1, 2, 
                      3, 2, 1,
                      2, 2, 3), 
                    nrow = 3, byrow = TRUE,
                    dimnames = list(weatherStates,weatherStates))
markovchainFit(weathersOfDays[1:200], method = "map", 
               confidencelevel = 0.92, hyperparam = hyperMatrix)
predictiveDistribution(weathersOfDays[1:200], 
                       weathersOfDays[201:365],hyperparam = hyperMatrix)

## ----MAPFit2-------------------------------------------------------------
hyperMatrix2<- hyperMatrix[c(2,3,1), c(2,3,1)]
markovchainFit(weathersOfDays[1:200], method = "map", 
               confidencelevel = 0.92, hyperparam = hyperMatrix2)
predictiveDistribution(weathersOfDays[1:200], 
                       weathersOfDays[201:365],hyperparam = hyperMatrix2)

## ----inferHyperparam-----------------------------------------------------
inferHyperparam(transMatr = weatherMatrix, scale = c(10, 10, 10))

## ----inferHyperparam2----------------------------------------------------
inferHyperparam(data = weathersOfDays[1:15])

## ----inferHyperparam3----------------------------------------------------
hyperMatrix3 <- inferHyperparam(transMatr = weatherMatrix, scale = c(10, 10, 10))
hyperMatrix3 <- hyperMatrix3$scaledInference

hyperMatrix4 <- inferHyperparam(data = weathersOfDays[1:15])
hyperMatrix4 <- hyperMatrix4$dataInference

## ----MAPandMLE-----------------------------------------------------------
data(preproglucacon, package = "markovchain")
preproglucacon <- preproglucacon[[2]]
MLEest <- markovchainFit(preproglucacon, method = "mle")
MAPest <- markovchainFit(preproglucacon, method = "map")
MLEest$estimate
MAPest$estimate

