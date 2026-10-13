# plot a diagram using diagram for a markovchain object
.plotdiagram <- function(object, ...) {
  # Works for discrete (markovchain) and continuous (ctmc) objects; the two classes store the
  # rates/probabilities in different slots (transitionMatrix vs generator).
  if(is(object,"markovchain")){
    mat <- object@transitionMatrix
    # Communicating classes = groups of states that can reach each other. Each class gets a
    # random colour so that the picture shows the structure of the chain.
    list <- .communicatingClassesRcpp(object)
    sections <- length(list)
    colorList <- grDevices::colors()
    colorList <- sample(colorList,sections)
    # Default colour is white; overwritten below for every state of every class.
    colorvector <- rep("white",length(object@states))
    for(i in 1:length(list)){
      part <- list[[i]]
      for(j in 1:length(part)){
        colorvector[match(part[j],object@states)] <- colorList[i]
      }
    }
  } else if(is(object,"ctmc")){
    mat <- object@generator
    colorvector <- rep("white",length(object@states))
  }
  # diagram::plotmat() expects the matrix the other way round, so make it row-oriented first.
  if(object@byrow == FALSE) {
    mat <- t(mat)
  }
  
  
  # diagram is only 'suggested': if it is missing we warn and return NULL instead of failing.
  if (!requireNamespace("diagram", quietly = TRUE)) {
    warning("package diagram is required to draw this plot and is not available: NULL is returned",
            call. = FALSE)
    return(invisible(NULL))
  } else{
    diagram::plotmat(t(mat),relsize = 0.75,box.col = colorvector, ...)  
  }
  
}

# plot a diagram using DiagrammeR for a markovchain object
.plotDiagrammeR <- function(object, ...) {
  # DiagrammeR draws the graph from a text description in the DOT language (see graphExport.R).
  if(is(object,"markovchain")){
  mat <- object@transitionMatrix
  } else if(is(object,"ctmc")){
    mat <- object@generator
  }
  # State names become node labels
  names <- rownames(mat)
  
  # Build the DOT text: first one node per state ...
  # names of nodes
  nodes <- ''
  for(i in 1:nrow(mat)) {
    nodes <- paste0(nodes, names[i], "; ")
  }
  
  # ... then one edge (arrow) per pair (i, j), labelled with its probability/rate.
  # store edges
  edges <- ''
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      edges <- paste0(edges, names[i], "->", names[j], " [label = ", mat[i,j], "] ")
    }
  }
  
  # Extra arguments given by the user (label=, labelloc=, ...) are appended as DOT attributes.
  # extract extra parameter
  dots <- list(...)
  args <- ""
  for(name in names(dots)) {
    args <- paste0(args, name, "=\"", dots[[name]], "\" ")
  }
  
  # print(args)
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
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
  } else {
    warning("package DiagrammeR is required to draw this plot and is not available: NULL is returned",
            call. = FALSE)
    return(invisible(NULL))
  }
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