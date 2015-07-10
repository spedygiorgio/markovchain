library(diagram)
library(DiagrammeR)

.plotdiagram <- function(object) {
  mat <- object@transitionMatrix
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

.plotDiagrammeR <- function(object) {
  mat <- object@transitionMatrix
  
  res <- DiagrammeR("
  graph LR;
             A(1)-->|0.5|B(2)
              A-->|0.5|A
              B-->|0.1|B
             B-->|0.3|C(3)
              C-->|0.9|B;
             
             style A fill:#A2EB86, stroke:#04C4AB, stroke-width:2px;
             style B fill:#FFF289, stroke:#FCFCFF, stroke-width:2px, stroke-dasharray: 4, 4;
             style C fill:#FFA070, stroke:#FF5E5E, stroke-width:2px;
             ")
  
  res <- grViz("
  digraph circles {
        
        # a 'graph' statement
        graph [overlap = true, fontsize = 10]
        
        # several 'node' statements
#         node [shape = box,
#         fontname = Helvetica]
#         A; B; C; D; E; F
        
        node [shape = circle,
        fixedsize = true,
        width = 0.9] // sets as circles
        1; 2; 3; 4; #5; 6; 7; 8
        
        # several 'edge' statements
        1->1 [label = 0.80] 1->2 [label = 0.1] 1->3 [label = 0.1] 
        2->1 [label = 0.75] 2->2 [label = 0.25]
        3->2 [label = 0.5] 3->3 [label='0.4'] 3->4 [label = 0.1]
        4->4 [label = 1.0] 
        
        # A->1 B->2 B->3 B->4 C->A
#         1->D E->A 2->4 1->5 1->F
#         E->6 4->6 5->7 6->7 3->8
  }
  ")
  return (res)
}

# mcWeather <- new("markovchain", states = c("sunny", "cloudy", "rain"),
#                  transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
#                                                     0.3, 0.4, 0.3,
#                                                     0.2, 0.45, 0.35), byrow = T, nrow = 3),
#                  name = "Weather")
# mcWeather
# # .plotdiagram(mcWeather)
# .plotDiagrammeR(mcWeather)
