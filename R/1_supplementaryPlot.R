# plot a diagram using diagram for a markovchain object
.plotdiagram <- function(object, ...) {
  if(class(object) == "markovchain"){
    mat <- object@transitionMatrix
  } else if(class(object) == "ctmc"){
    mat <- object@generator
  }
  if(object@byrow == FALSE) {
    mat <- t(mat)
  }
  
  # pass the matrix as columnwise fashion
  diagram::plotmat(t(mat),relsize = 0.75, ...)
}

# plot a diagram using DiagrammeR for a markovchain object
.plotDiagrammeR <- function(object, ...) {
  if(class(object) == "markovchain"){
  mat <- object@transitionMatrix
  } else if(class(object) == "ctmc"){
    mat <- object@generator
  }
  names <- rownames(mat)
  
  # names of nodes
  nodes <- ''
  for(i in 1:nrow(mat)) {
    nodes <- paste0(nodes, names[i], "; ")
  }
  
  # store edges
  edges <- ''
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      edges <- paste0(edges, names[i], "->", names[j], " [label = ", mat[i,j], "] ")
    }
  }
  
  # extract extra parameter
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
  
  return (res)
}

#  How to do plotting?
#  mcWeather <- new("markovchain", states = c("sunny", "cloudy", "rain"),
#                   transitionMatrix = matrix(data = c(0.70,  0.2,  0.1,
#                                                       0.3,  0.4,  0.3,
#                                                       0.2, 0.45, 0.35), byrow = T, nrow = 3),
#                   name = "Weather")
#  mcWeather
#  .plotdiagram(mcWeather, box.size = 0.06)
#  .plotDiagrammeR(mcWeather, label ="Weather transition matrix", labelloc = "t")
#   plot(mcWeather, package = "DiagrammeR", label = "Weather transition matrix")