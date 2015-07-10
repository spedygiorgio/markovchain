library(diagram)

.plotdiagram <- function(object) {
  mat = object@transitionMatrix
  # print(mat)
  # plotmat(mat)
  res <- plotmat(mat, box.size = 0.06)
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
  return (res)
}

# mcWeather <- new("markovchain", states = c("sunny", "cloudy", "rain"),
#                  transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
#                                                     0.3, 0.4, 0.3,
#                                                     0.2, 0.45, 0.35), byrow = T, nrow = 3),
#                  name = "Weather")
# mcWeather
# .plotdiagram(mcWeather)
