library(diagram)
library(DiagrammeR)

.plotdiagram <- function(object, options = NULL) {
  mat <- object@transitionMatrix
  res <- plotmat(mat, options)
  return (res)
}

.plotDiagrammeR <- function(object, options = NULL) {
  mat <- object@transitionMatrix
  nodes <- ''
  for(i in 1:nrow(mat)) {
    nodes<- paste0(nodes, i, "; ")
  }
  # print(nodes)
  edges <- ''
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      # 1->1 [label = 0.80]
      edges <- paste0(edges, i, "->", j, " [label = ", mat[i,j], "] ")
    }
  }
  # print(edges)
#   res <- DiagrammeR("
#   graph LR;
#              A(1)-->|0.5|B(2)
#               A-->|0.5|A
#               B-->|0.1|B
#              B-->|0.3|C(3)
#               C-->|0.9|B;
#              
#              style A fill:#A2EB86, stroke:#04C4AB, stroke-width:2px;
#              style B fill:#FFF289, stroke:#FCFCFF, stroke-width:2px, stroke-dasharray: 4, 4;
#              style C fill:#FFA070, stroke:#FF5E5E, stroke-width:2px;
#              ")
  
  res <- grViz(paste0("
  digraph circles {
        # a 'graph' statement
        graph [overlap = true, fontsize = 10]

        node [shape = circle,
        fixedsize = true,
        width = 0.9] // sets as circles
        ", nodes, "
        
        # several 'edge' statements
        ", edges, " 
  }
  "))
  
#   res <- grViz("
#   digraph circles {
#         # a 'graph' statement
#         graph [overlap = true, fontsize = 10]
#         
#         # several 'node' statements
# #         node [shape = box,
# #         fontname = Helvetica]
# #         A; B; C; D; E; F
#         
#         node [shape = circle,
#         fixedsize = true,
#         width = 0.9] // sets as circles
#         1; 2; 3; 4; 
#         
#         # several 'edge' statements
#         1->1 [label = 0.80] 1->2 [label = 0.1] 1->3 [label = 0.1] 
#         2->1 [label = 0.75] 2->2 [label = 0.25]
#         3->2 [label = 0.5] 3->3 [label='0.4'] 3->4 [label = 0.1]
#         4->4 [label = 1.0] 
#   }
#   ")
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
