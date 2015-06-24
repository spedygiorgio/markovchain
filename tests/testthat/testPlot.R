library(markovchain)

P <- matrix(c(0, 0.5, 0.5,
              0.5, 0, 0.5, 
              0.5, 0.5, 0), byrow=T, ncol=3)
P
k <- new("markovchain", transitionMatrix=P)
plot(k) 

P <- matrix(c(0, 0.5, 0.5,
              0.5, 0, 0.5, 
              0.4, 0.6, 0), byrow=T, ncol=3)
P
k <- new("markovchain", transitionMatrix=P)
plot(k) 

P <- matrix(c(0, 0.5, 0.5,
              0.5, 0, 0.5, 
              0.5, 0.4, 0.1), byrow=T, ncol=3)
P
k <- new("markovchain", transitionMatrix=P)
# plot(k) 

P <- matrix(c(0.1, 0.4, 0.5,
              0.4, 0.1, 0.5, 
              0.5, 0.3, 0.2), byrow=T, ncol=3)
P
k <- new("markovchain", states = c("1", "2", "3"), transitionMatrix=P, name = "test")
# plot(k) 

mcWeather <- new("markovchain", states = c("sunny", "cloudy", "rain"),
                  transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
                                                       0.3, 0.4, 0.3,
                                                       0.2, 0.45, 0.35), byrow = T, nrow = 3),
                  name = "Weather")
# mcWeather
# plot(mcWeather)
