#library(markovchain)
require(diagram)
require(DiagrammeR)
# P <- matrix(c(0, 0.5, 0.5,
#               0.5, 0, 0.5, 
#               0.5, 0.5, 0), byrow=T, ncol=3)
# P
# k <- new("markovchain", transitionMatrix=P)
# plot(k) 

# P <- matrix(c(0, 0.5, 0.5,
#               0.5, 0, 0.5, 
#               0.4, 0.6, 0), byrow=T, ncol=3)
# P
# k <- new("markovchain", transitionMatrix=P)
# plot(k) 

# P <- matrix(c(0, 0.5, 0.5,
#               0.5, 0, 0.5, 
#               0.5, 0.4, 0.1), byrow=T, ncol=3)
# P
# k <- new("markovchain", transitionMatrix=P)
# plot(k) 

# P <- matrix(c(0.1, 0.4, 0.5,
#               0.4, 0.1, 0.5, 
#               0.5, 0.3, 0.2), byrow=T, ncol=3)
# P
# k <- new("markovchain", states = c("1", "2", "3"), transitionMatrix=P, name = "test")
# plot(k) 

mcWeather <- new("markovchain", states = c("sunny", "cloudy", "rain"),
                  transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
                                                       0.3, 0.4, 0.3,
                                                       0.2, 0.45, 0.35), byrow = T, nrow = 3),
                  name = "Weather")
mcWeather
plot(mcWeather)
plot(mcWeather, package = "diagram", box.size = 0.06, main = "Weather transition matrix")
plot(mcWeather, package = "DiagrammeR", label ="Weather transition matrix", labelloc="t")

#   curves <- matrix(nrow = ncol(mat), ncol = ncol(mat), 0)
#   curves[3, 1] <- curves[1, 6] <- -0.35
#   curves[4, 6] <- curves[6, 4] <- curves[5, 6] <- curves[6, 5] <- 0.08
#   curves[3, 6] <- 0.35
#   plotmat(mat, pos = c(3, 2, 1), curve = curves,
#             name = colnames(mat), lwd = 1, box.lwd = 2,
#             cex.txt = 0.8, box.cex = 0.8, box.size = 0.08,
#             arr.length = 0.5, box.type = "circle", box.prop = 1,
#             shadow.size = 0.01, self.cex = 0.6, my = -0.075, mx = -0.01,
#             relsize = 0.9, self.shiftx = c(0, 0, 0.125, -0.12, 0.125, 0),
#             self.shifty = 0, main = "Diagram")
