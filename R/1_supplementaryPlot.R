# library(diagram)
# library(DiagrammeR)

.plotdiagram <- function(object, ...) {
  mat <- object@transitionMatrix
  #res <- plotmat(mat, ...)
  diagram::plotmat(mat, ...)
  #return (res)
}

.plotDiagrammeR <- function(object, ...) {
  mat <- object@transitionMatrix
  names <- rownames(mat)
  nodes <- ''
  for(i in 1:nrow(mat)) {
    nodes<- paste0(nodes, names[i], "; ")
  }
  edges <- ''
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      edges <- paste0(edges, names[i], "->", names[j], " [label = ", mat[i,j], "] ")
    }
  }
  
  dots <- list(...)
  args <- ""
  for(name in names(dots)) {
    args <- paste0(args, name, "=\"", dots[[name]], "\" ")
  }
  # print(args)
  res <- DiagrammeR::grViz(paste0("
  digraph circles {
        graph [overlap = true, fontsize = 10]

        node [shape = circle,
        fixedsize = true,
        width = 0.9] // sets as circles
        ", nodes, "
        
        ", edges, args,"
// labelfontsize = 20 labelloc='t' label ='Weather transition matrix'
  }
  "))
  
#   res <- grViz("
#   digraph circles {
#         # a 'graph' statement
#         graph [overlap = true, fontsize = 10]
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
# .plotdiagram(mcWeather, box.size = 0.06)
# .plotDiagrammeR(mcWeather, label ="Weather transition matrix", labelloc="t")
# plot(mcWeather, package = "DiagrammeR", label ="Weather transition matrix")
