#' @title Markov Chain class
#' @name markovchain-class
#' @aliases markovchain-class *,markovchain,markovchain-method
#'   *,markovchain,matrix-method *,markovchain,numeric-method
#'   *,matrix,markovchain-method *,numeric,markovchain-method
#'   ==,markovchain,markovchain-method !=,markovchain,markovchain-method
#'   absorbingStates,markovchain-method transientStates,markovchain-method
#'   recurrentStates,markovchain-method transientClasses,markovchain-method
#'   recurrentClasses,markovchain-method communicatingClasses,markovchain-method
#'   steadyStates,markovchain-method meanNumVisits,markovchain-method
#'   is.regular,markovchain-method is.irreducible,markovchain-method
#'   is.accessible,markovchain,character,character-method
#'   is.accessible,markovchain,missing,missing-method
#'   absorptionProbabilities,markovchain-method
#'   meanFirstPassageTime,markovchain,character-method 
#'   meanFirstPassageTime,markovchain,missing-method
#'   meanAbsorptionTime,markovchain-method
#'   meanRecurrenceTime,markovchain-method
#'   conditionalDistribution,markovchain-method hittingProbabilities,markovchain-method
#'   canonicForm,markovchain-method coerce,data.frame,markovchain-method
#'   coerce,markovchain,data.frame-method coerce,table,markovchain-method
#'   coerce,markovchain,igraph-method coerce,markovchain,matrix-method
#'   coerce,markovchain,sparseMatrix-method coerce,sparseMatrix,markovchain-method
#'   coerce,matrix,markovchain-method coerce,msm,markovchain-method
#'   coerce,msm.est,markovchain-method coerce,etm,markovchain-method
#'   dim,markovchain-method initialize,markovchain-method
#'   names<-,markovchain-method plot,markovchain,missing-method
#'   predict,markovchain-method print,markovchain-method
#'   show,markovchain-method summary,markovchain-method
#'   sort,markovchain-method t,markovchain-method
#'   [,markovchain,ANY,ANY,ANY-method ^,markovchain,numeric-method
#' @description The S4 class that describes \code{markovchain} objects.
#' 
#' @param states Name of the states. Must be the same of \code{colnames} and \code{rownames} of the transition matrix
#' @param byrow TRUE or FALSE indicating whether the supplied matrix 
#'   is either stochastic by rows or by columns
#' @param transitionMatrix Square transition matrix
#' @param name Optional character name of the Markov chain
#' 
#' @section Creation of objects:
#' 
#' Objects can be created by calls of the form \code{new("markovchain", states, byrow, transitionMatrix, ...)}.
#' 
#' @section Methods:
#' 
#' \describe{
#'    \item{*}{\code{signature(e1 = "markovchain", e2 = "markovchain")}: multiply two \code{markovchain} objects}
#'    \item{*}{\code{signature(e1 = "markovchain", e2 = "matrix")}: markovchain by matrix multiplication}
#'    \item{*}{\code{signature(e1 = "markovchain", e2 = "numeric")}: markovchain by numeric vector multiplication }
#'    \item{*}{\code{signature(e1 = "matrix", e2 = "markovchain")}: matrix by markov chain}
#'    \item{*}{\code{signature(e1 = "numeric", e2 = "markovchain")}: numeric vector by \code{markovchain} multiplication   }
#'    \item{[}{\code{signature(x = "markovchain", i = "ANY", j = "ANY", drop = "ANY")}: ... }
#'    \item{^}{\code{signature(e1 = "markovchain", e2 = "numeric")}: power of a \code{markovchain} object}
#'    \item{==}{\code{signature(e1 = "markovchain", e2 = "markovchain")}: equality of two \code{markovchain} object}
#'    \item{!=}{\code{signature(e1 = "markovchain", e2 = "markovchain")}: non-equality of two \code{markovchain} object}
#'    \item{absorbingStates}{\code{signature(object = "markovchain")}: method to get absorbing states }
#'    \item{canonicForm}{\code{signature(object = "markovchain")}: return a \code{markovchain} object into canonic form }
#'    \item{coerce}{\code{signature(from = "markovchain", to = "data.frame")}: coerce method from markovchain to \code{data.frame}}
#'    \item{conditionalDistribution}{\code{signature(object = "markovchain")}: returns the conditional probability of subsequent states given a state}
#'    \item{coerce}{\code{signature(from = "data.frame", to = "markovchain")}: coerce method from \code{data.frame} to \code{markovchain}}
#'    \item{coerce}{\code{signature(from = "table", to = "markovchain")}: coerce method from \code{table} to \code{markovchain} }
#'    \item{coerce}{\code{signature(from = "msm", to = "markovchain")}: coerce method from \code{msm} to \code{markovchain} }
#'    \item{coerce}{\code{signature(from = "msm.est", to = "markovchain")}: coerce method from \code{msm.est} (but only from a Probability Matrix) to \code{markovchain} }
#'    \item{coerce}{\code{signature(from = "etm", to = "markovchain")}: coerce method from \code{etm} to \code{markovchain} }
#'    \item{coerce}{\code{signature(from = "sparseMatrix", to = "markovchain")}: coerce method from \code{sparseMatrix} to \code{markovchain} }
#'    \item{coerce}{\code{signature(from = "markovchain", to = "igraph")}: coercing to \code{igraph} objects }
#'    \item{coerce}{\code{signature(from = "markovchain", to = "matrix")}: coercing to \code{matrix} objects }
#'    \item{coerce}{\code{signature(from = "markovchain", to = "sparseMatrix")}: coercing to \code{sparseMatrix} objects }
#'    \item{coerce}{\code{signature(from = "matrix", to = "markovchain")}: coercing to \code{markovchain} objects from \code{matrix} one }
#'    \item{dim}{\code{signature(x = "markovchain")}: method to get the size}
#'    \item{names}{\code{signature(x = "markovchain")}: method to get the names of states}
#'    \item{names<-}{\code{signature(x = "markovchain", value = "character")}: method to set the names of states}
#'    \item{initialize}{\code{signature(.Object = "markovchain")}: initialize method }
#'    \item{plot}{\code{signature(x = "markovchain", y = "missing")}: plot method for \code{markovchain} objects }
#'    \item{predict}{\code{signature(object = "markovchain")}: predict method }
#'    \item{print}{\code{signature(x = "markovchain")}: print method. }
#'    \item{show}{\code{signature(object = "markovchain")}: show method. }
#'    \item{sort}{\code{signature(x = "markovchain", decreasing=FALSE)}: sorting the transition matrix. }
#'    \item{states}{\code{signature(object = "markovchain")}: returns the names of states (as \code{names}. }
#'    \item{steadyStates}{\code{signature(object = "markovchain")}: method to get the steady vector. }
#'    \item{summary}{\code{signature(object = "markovchain")}: method to summarize structure of the markov chain }
#'    \item{transientStates}{\code{signature(object = "markovchain")}: method to get the transient states. }
#'    \item{t}{\code{signature(x = "markovchain")}: transpose matrix }
#'    \item{transitionProbability}{\code{signature(object = "markovchain")}: transition probability }
#' }
#' 
#' @references 
#' A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
#' 
#' @author Giorgio Spedicato
#' @note  
#' \enumerate{
#' \item \code{markovchain} object are backed by S4 Classes.
#' \item Validation method is used to assess whether either columns or rows totals to one. 
#' Rounding is used up to \code{.Machine$double.eps * 100}. If state names are not properly
#' defined for a probability  \code{matrix}, coercing to \code{markovhcain} object leads 
#' to overriding states name with artificial "s1", "s2", ... sequence. In addition, operator
#' overloading has been applied for \eqn{+,*,^,==,!=} operators.
#' }
#' 
#' @seealso \code{\link{markovchainSequence}},\code{\link{markovchainFit}}
#' 
#' @examples
#' #show markovchain definition
#' showClass("markovchain")
#' #create a simple Markov chain
#' transMatr<-matrix(c(0.4,0.6,.3,.7),nrow=2,byrow=TRUE)
#' simpleMc<-new("markovchain", states=c("a","b"),
#'               transitionMatrix=transMatr, 
#'               name="simpleMc")
#' #power
#' simpleMc^4
#' #some methods
#' steadyStates(simpleMc)
#' absorbingStates(simpleMc)
#' simpleMc[2,1]
#' t(simpleMc)
#' is.irreducible(simpleMc)
#' #conditional distributions
#' conditionalDistribution(simpleMc, "b")
#' #example for predict method
#' sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
#' mcFit<-markovchainFit(data=sequence)
#' predict(mcFit$estimate, newdata="b",n.ahead=3)
#' #direct conversion
#' myMc<-as(transMatr, "markovchain")
#' 
#' #example of summary
#' summary(simpleMc)
#' \dontrun{plot(simpleMc)}
#' 
#' @keywords classes
#' 
#' @export
setClass(
  # Class name
  "markovchain",
  # Define the slots
  slots = list(states = "character", byrow = "logical",
  transitionMatrix = "matrix", name = "character"),
  # Set the default values for the slots
  prototype = list(
    states = c("a", "b"), 
    byrow = TRUE,
    transitionMatrix = matrix(
      data = c(0, 1, 1, 0), 
      nrow = 2, 
      byrow = TRUE, 
      dimnames = list(c("a", "b"), c("a", "b"))), 
    name = "Unnamed Markov chain")
)

# Initializing method for markovchain objects
setMethod(
  "initialize",
  signature(.Object = "markovchain"),
  function (.Object, states, byrow, transitionMatrix, name, ...) {
    # Put the standard markovchain
    if (missing(transitionMatrix)) {
      transitionMatrix <- matrix(
        data = c(0, 1, 1, 0),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("a", "b")))
    }
    
    rowNames <- rownames(transitionMatrix)
    colNames <- colnames(transitionMatrix)
    
    # Check names of transition matrix
    # if all names are missing it initializes them to "1", "2", ....
    if (all(is.null(rowNames), is.null(colNames)) == TRUE) {
      if (missing(states)) {
        numRows <- nrow(transitionMatrix)
        stateNames <- as.character(seq(1:numRows))
      } else {
        stateNames <- states
      }
      
      rownames(transitionMatrix) <- stateNames
      colnames(transitionMatrix) <- stateNames
      
    # Fix when rownames null
    } else if (is.null(rowNames)) {
      rownames(transitionMatrix) <- colNames
    # Fix when colnames null
    } else if (is.null(colNames)) {
      colnames(transitionMatrix) <- rowNames
    # Fix when different
    } else if (! setequal(rowNames, colNames)) {
      colnames(transitionMatrix) <- rowNames
    }
    
    if (missing(states))
      states <- rownames(transitionMatrix)
    
    if (missing(byrow))
      byrow <- TRUE
    
    if (missing(name))
      name <- "Unnamed Markov chain"
    
    callNextMethod(
      .Object,
      states = states,
      byrow = byrow,
      transitionMatrix = transitionMatrix,
      name = name,
      ...
    )
  }
)

#' @title Non homogeneus discrete time Markov Chains class
#' @name markovchainList-class
#' @aliases [[,markovchainList-method dim,markovchainList-method
#'   predict,markovchainList-method print,markovchainList-method
#'   show,markovchainList-method
#' @description A class to handle non homogeneous discrete Markov chains
#' 
#' @param markovchains Object of class \code{"list"}: a list of markovchains
#' @param name Object of class \code{"character"}: optional name of the class
#' 
#' @section Objects from the Class:
#'
#'   A \code{markovchainlist} is a list of \code{markovchain} objects. They can
#'   be used to model non homogeneous discrete time Markov Chains, when
#'   transition probabilities (and possible states) change by time.
#' @section Methods:
#' \describe{
#' \item{[[}{\code{signature(x = "markovchainList")}: extract the
#' i-th \code{markovchain} }
#' \item{dim}{\code{signature(x = "markovchainList")}: number 
#' of \code{markovchain} underlying the matrix }
#' \item{predict}{\code{signature(object = "markovchainList")}: predict 
#' from a \code{markovchainList} }
#' \item{print}{\code{signature(x = "markovchainList")}: prints the list 
#'   of markovchains }
#' \item{show}{\code{signature(object = "markovchainList")}: same as \code{print} }
#' }
#' 
#' @references 
#' A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
#' 
#' @author  Giorgio Spedicato
#' 
#' @note 
#' The class consists in a list of \code{markovchain} objects. 
#' It is aimed at working with non homogeneous Markov chains.
#' 
#' @seealso \code{\linkS4class{markovchain}}
#' @examples
#' showClass("markovchainList")
#' #define a markovchainList
#' statesNames=c("a","b")
#' 
#' mcA<-new("markovchain",name="MCA", 
#'          transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),
#'                           byrow=TRUE, nrow=2, 
#'                           dimnames=list(statesNames,statesNames))
#'         )
#'                                                            
#' mcB<-new("markovchain", states=c("a","b","c"), name="MCB",
#'          transitionMatrix=matrix(c(0.2,0.5,0.3,0,1,0,0.1,0.8,0.1),
#'          nrow=3, byrow=TRUE))
#'  
#' mcC<-new("markovchain", states=c("a","b","c","d"), name="MCC",
#'          transitionMatrix=matrix(c(0.25,0.75,0,0,0.4,0.6,
#'                                    0,0,0,0,0.1,0.9,0,0,0.7,0.3), 
#'                                  nrow=4, byrow=TRUE)
#' )
#' mcList<-new("markovchainList",markovchains=list(mcA, mcB, mcC), 
#'            name="Non - homogeneous Markov Chain")
#' 
#' @keywords classes
#' 
#' @export
setClass(
  "markovchainList",
  slots = list(
    markovchains = "list",
    name = "character")
)

# Verifies whether a markovchainList object is valid or not
# A markovchainList is valid iff all the slots are markovchain objects
# Returns true if the markovchainList is valid, the indexes of the
# wrong slots otherwise
setValidity(
  "markovchainList",
  function(object) {
    check <- FALSE
    markovchains <- object@markovchains
    
    classes <- sapply(markovchains, class)
    nonMarkovchain <- which(classes != "markovchain")
    errors <- sapply(nonMarkovchain, function(i) {
      paste(i, "-th element class is not 'markovchain'")
    })
    
    if (length(errors) == 0) TRUE else errors
  }
)

# generic method to print out states

#' @name states
#' 
#' @title Defined states of a transition matrix
#' 
#' @description This method returns the states of a transition matrix.
#' 
#' @param object A discrete \code{markovchain} object
#' @return The character vector corresponding to states slot.
#' 
#' @references A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
#' 
#' @author Giorgio Spedicato
#' 
#' @seealso \code{\linkS4class{markovchain}}
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#'                 matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
#'                 byrow = TRUE, dimnames=list(statesNames,statesNames)),
#'                 name = "A markovchain Object" 
#' )
#' states(markovB)
#' names(markovB)
#' 
#' @rdname states
#' 
#' @export
setGeneric("states", function(object) standardGeneric("states"))

#' @rdname states
#' @title states
setMethod(
  "states",
  "markovchain", 
  function(object) {
    object@states
  }
)

#' @title Returns the states for a Markov chain object
#'
#' @param x object we want to return states for
#' 
#' @rdname names
setMethod(
  "names",
  "markovchain", 
  function(x) {
    x@states
  }
)

#' @title Method to retrieve name of markovchain object  
#' 
#' @name name
#' 
#' @description This method returns the name of a markovchain object
#' 
#' @param object A markovchain object
#' @rdname getName
#' @author Giorgio Spedicato, Deepak Yadav
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#'                 matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
#'                 byrow = TRUE, dimnames=list(statesNames,statesNames)),
#'                 name = "A markovchain Object" 
#' )
#' name(markovB)
#' 
#' @export
setGeneric("name", function(object) standardGeneric("name"))


#' @rdname getName
setMethod(
  "name", 
  "markovchain", 
  function(object) {
    object@name
})

#' @title Method to set name of markovchain object
#' 
#' @name name<-
#' 
#' @description This method modifies the existing name of markovchain object
#' 
#' @param object A markovchain object
#' @param value New name of markovchain object
#' @rdname setName
#' @author Giorgio Spedicato, Deepak Yadav
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#'                 matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
#'                 byrow = TRUE, dimnames=list(statesNames,statesNames)),
#'                 name = "A markovchain Object" 
#' )
#' name(markovB) <- "dangerous mc"
#' 
#' @export
setGeneric("name<-", function(object, value) standardGeneric("name<-"))

#' @rdname setName
setMethod(
  "name<-", 
  "markovchain", 
  function(object, value) {
    object@name <- value
    object
  }
)

setMethod(
  "names<-", 
  "markovchain", 
  function(x, value) {
    rownames(x@transitionMatrix) <- value
    colnames(x@transitionMatrix) <- value
    x@states <- value
    x
  }
)


#' @exportMethod dim
setGeneric("dim")

# Generic methods to get the dim of a markovchain and markovchainList
setMethod(
  "dim",
  "markovchain", 
  function(x) {
    nrow(x@transitionMatrix)
  }
)

setMethod(
  "dim",
  "markovchainList", 
  function(x) {
    length(x@markovchains)
  }
)


# method  to set the validity of a markovchain object
setValidity(
  "markovchain",
  function(object) {
    errors <- character()
    transitionMatrix <- object@transitionMatrix
    states           <- object@states
    
    if (length(setdiff(states, unique(states))) > 0) {
      msg    <- "Error! States must be unique!"
      errors <- c(errors, msg)
    }
    # Performs a set of checks. If any error arises, it ends up concatenated to errors
    
    # Check all values of transition matrix belongs to [0, 1]
    maybeProbabilities <- sapply(as.numeric(transitionMatrix), .isProbability)
    
    if (any(maybeProbabilities) == FALSE) {
      msg    <- "Error! Some elements of transitionMatrix are not probabilities"
      errors <- c(errors, msg)
    }
    
    # Check whether matrix is square matrix or not
    if (nrow(transitionMatrix) != ncol(transitionMatrix)) {
      msg    <- "Error! transitionMatrix is not a square matrix"
      errors <- c(errors, msg)
    }
    
    if (!.checkMatrix(transitionMatrix, object@byrow)) {
      msg <- paste(
               paste("Error!", 
               ifelse(object@byrow, "Rows", "Cols")),
               "of transition matrix do not sum to one"
             )
      errors <- c(errors, msg)
    }
    
    # Check whether column names or rows names equal to state names or not
    if (! setequal(colnames(transitionMatrix), states)) {
      msg    <- "Error! Colnames of transitionMatrix do not match states"
      errors <- c(errors, msg)
    }
    if (! setequal(rownames(transitionMatrix), states)) {
      msg    <- "Error! Rownames of transitionMatrix do not match states"
      errors <- c(errors, msg)
    }
    
    if (length(errors) > 0) errors else TRUE
  }
)


# generic method to extract transition probability
# from state t0 to state t1

#' @name transitionProbability
#' @title Function to get the transition probabilities from initial 
#'        to subsequent states.
#' @description This is a convenience function to get transition probabilities.
#' 
#' @param object A \code{markovchain} object.
#' @param t0 Initial state.
#' @param t1 Subsequent state.
#' 
#' @references A First Course in Probability (8th Edition), 
#'             Sheldon Ross, Prentice Hall 2010
#' 
#' @return Numeric Vector  
#' 
#' @author Giorgio Spedicato
#' @seealso \code{\linkS4class{markovchain}}
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#'                 matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
#'                 byrow = TRUE, dimnames=list(statesNames,statesNames)),
#'                name = "A markovchain Object" 
#' )    
#' transitionProbability(markovB,"b", "c")
#' @rdname transitionProbability
#'      
#' @exportMethod transitionProbability
setGeneric("transitionProbability", function(object, t0, t1) standardGeneric("transitionProbability"))

#' @rdname transitionProbability
setMethod("transitionProbability", "markovchain", 
  function(object, t0, t1) {
    fromState <- which(object@states == t0)
    toState <- which(object@states == t1)
    out <- ifelse(object@byrow == TRUE, object@transitionMatrix[fromState, toState] , 
                  object@transitionMatrix[toState, fromState])
    return(out)
  }
)

#  print, plot and show methods

.showInt <- function(object, verbose = TRUE) {
	
  # find the direction
  if (object@byrow == TRUE) {
	  direction <- "(by rows)" 
	} else {
	  direction <- "(by cols)" 
	}
  
	if (verbose == TRUE) {
	  cat(object@name, "\n A ", dim(object), "- dimensional discrete Markov Chain defined by the following states: \n",
	      paste(states(object), collapse=", "), "\n The transition matrix ", 
	      direction, " is defined as follows: \n")
	}
  
	print(object@transitionMatrix)
	cat("\n")
}

#' @exportMethod show
setGeneric("show")

# show methods for markovchain and markovchain list objects 
setMethod("show", "markovchain",
  function(object){
    .showInt(object)
  }
)

setMethod("show", "markovchainList",
  function(object) {
    cat(object@name, " list of Markov chain(s)", "\n")
    for(i in 1:length(object@markovchains)) {
      cat("Markovchain ",i,"\n")
      show(object@markovchains[[i]])
    }
  }
)

#' @exportMethod print
setGeneric("print")

# print methods
setMethod("print", "markovchainList", function(x) show(x))
setMethod("print", "markovchain",
  function(x){
   object <- x
   .showInt(object, verbose = FALSE)
  }
)

.getNet <- function(object, round = FALSE) {
 	
  # function to get the absorbency matrix to plot and export to igraph
	#
 	# Args: 
	# object: a markovchain object
	# round: boolean to round
 	#
	# Returns:
	#
	# a graph adjacency
  
	if (object@byrow == FALSE) {
	  object <- t(object)
	}
  
	matr <- object@transitionMatrix*100
	if(round == TRUE) {
	  matr <- round(matr, 2)
	}
	
	net <- graph.adjacency(adjmatrix = matr, weighted = TRUE, mode = "directed")
	return(net)
}

getColorVector <- function(object){
  list <- .communicatingClassesRcpp(object)
  sections <- length(list)
  colorList <- grDevices::colors()
  colorList <- sample(colorList,sections)
  colorvector <- rep("white",length(object@states))
  for(i in 1:length(list)){
    part <- list[[i]]
    for(j in 1:length(part)){
      colorvector[match(part[j],object@states)] <- colorList[i]
    }
  }
  return(colorvector)
}


#' @exportMethod plot
setGeneric("plot")

# Plot methods for markovchain objects

# plot method from stat5
setMethod("plot", signature(x = "markovchain", y = "missing"),
  function(x, y, package = "igraph", ...) {
    switch(package,
     diagram = {
       if (requireNamespace("diagram", quietly = TRUE)) {
         .plotdiagram(object = x, ...)
       } else {
         netMc <- .getNet(object = x, round = TRUE)
         edgeLabel <- round(E(netMc)$weight / 100, 2)
         plot.igraph(x = netMc, edge.label = edgeLabel, ...)
       }
     },
     
     DiagrammeR = {
       if (requireNamespace("DiagrammeR", quietly = TRUE)) {
         .plotDiagrammeR(object = x, ...)
       } else {
         netMc <- .getNet(object = x, round = TRUE)
         edgeLabel <- round(E(netMc)$weight / 100, 2)
         plot.igraph(x = netMc, edge.label = edgeLabel, ...)
       }
     },
     {
       netMc <- .getNet(object = x,round = TRUE)
       edgeLabel <- round(E(netMc)$weight / 100, 2)
       plot.igraph(x = netMc, edge.label = edgeLabel, ...)
    })
  }
)


##################################################AS METHODS#########################

.checkMatrix <- function(matr, byrow = TRUE, verbose = FALSE) {
	
  # firstly, check size
	if (ncol(matr) != nrow(matr)) {
		if(verbose) stop("Error! Not a rectangular matrix")
		return(FALSE)
	}
  
  # secondly, check is stochastic
  isStochastic <- .isStochasticMatrix(matr, byrow)
  
  if (!isStochastic) {
    if (verbose)
	    stop("Error! Either rows or cols should sum to 1")
  
	  return(FALSE)
  }
  
	# if all test are passed
	return(TRUE)
}

# Internal function to return a markovchain object given a matrix
.matrix2Mc <- function(from) {
	
  # whether given matrix is a transition matrix or not
  # if it is then how probabilities are stored
  # row-wise or columnwise
  
	byrow <- TRUE
	checkByRows <- .checkMatrix(from, byrow = byrow)
	
	if(!checkByRows) {
	  byrow <- FALSE
		checkByCols <- .checkMatrix(from, byrow = byrow)
		
		if(!checkByCols) {
		  #error could be either in rows or in cols
		  if (any(colSums(from) != 1)) cat("columns sums not equal to one are:", which(colSums(from) != 1),"\n")
		  if (any(rowSums(from) != 1)) cat("row sums not equal to one are:", which(rowSums(from) != 1),"\n")
		  stop("Error! Not a transition matrix")	
		}
	}
	
	# extract states names
	if(byrow) {
	  namesCandidate <- rownames(from) 
	} else {
	  namesCandidate <- colnames(from)
	}
	
	# if states names is not there create it s1, s2, s3, ....
	if(is.null(namesCandidate)) {
		namesCandidate <- paste("s", 1:nrow(from), sep = "")
	}
	
	# create markovchain object
	out <- new("markovchain", transitionMatrix = from, states = namesCandidate, byrow = byrow)
	
	invisible(out)
}

#' @exportMethod coerce
NULL

# coerce matrix to markovchain object using internal method
# example: as("some matrix", "markovchain")
setAs(from = "matrix", to = "markovchain", def = .matrix2Mc)

# Function to transform a markovchain into a data.frame
# Args:
# from: a markovchain object
#
# returns:
# a data.frame
.mc2Df <- function(from) {
  
  # number of rows or columns
	nr <- nrow(from@transitionMatrix)
	for(i in 1:nr){
		for(j in 1:nr){
			t0 <- from@states[i]
			t1 <- from@states[j]
			prob <- transitionProbability(object = from, t0 = t0, t1 = t1)
			#cope with the new default of R 4.0 (5-3-2020)
			rowDf <- data.frame(t0 = t0, t1 = t1, prob = prob,stringsAsFactors = TRUE )
			
			# go to else part if first row of data frame is generated
			if(exists("outDf")) {
			  outDf <- rbind(outDf, rowDf)
			} else {
			  outDf <- rowDf
			}
		}
	}
	
	return(outDf)
}

# method to convert(coerce) from markovchain to data.frame
setAs(from = "markovchain", to = "data.frame", def = .mc2Df)

# method to find the column which stores transition probability
.whichColProb <- function(df) {
	
  # column number which stores transition probability
  out <- 0
  
  # check for validity of data frame
	if(ncol(df) > 3) {
	  warning("Warning! More than three columns. Only the first three will be used")
	}
  
	if(ncol(df) < 3) {
	  stop("Error! Three columns needed")
	}
	
	for(i in 1:ncol(df)) {
	    
	  # when found the first numeric and probability col
			if((is(df[, i], "numeric")) & (all(sapply(df[, i], .isProbability) == TRUE))) {
					out <- i
					break
			}
	}
  
	return(out)
}

# Function to convert from a data.frame containing initial, ending 
#    and probability columns to a proper markovchain object
#
# Args:
# from: a data.frame
#
# Returns:
# A markovchain object 

.df2Mc <- function(from) {
	
	statesNames <- unique(from[, 1])
	colProb <- .whichColProb(from) # what is the use
	
	# transition matrix
	prMatr <- zeros(length(statesNames))
	rownames(prMatr) <- statesNames
	colnames(prMatr) <- statesNames
	
	
	for(i in 1:nrow(from)) {
		idRow <- which(statesNames == from[i, 1]) # assume first col from
		idCol <- which(statesNames == from[i, 2]) # assume second col to
		prMatr[idRow, idCol] <- from[i, 3]        # assume third col t-probability
	}
	
 	out <- new("markovchain", transitionMatrix = prMatr)
	return(out)
}

# method to convert(coerce) data frame to markovchain object 
setAs(from = "data.frame", to = "markovchain", def = .df2Mc)


# example
# data <- data.frame(from = c("a", "a", "b", "b", "b", "b"), 
#                      to = c("a", "b", "b", "b", "b", "a"))
# 
# from <- table(data)
# .table2Mc(from)

.table2Mc <- function(from) {
	
  # check whether table has square dimension or not
	if(dim(from)[1] != dim(from)[2]) {
	  stop("Error! Table is not squared")
	}
  
  # rows ond columns name should be same
	if(!setequal(rownames(from),colnames(from))) {
	  stop("Error! Rows not equal to coulumns")
	}
  
  temp <- unclass(as.matrix(from))
  
  # make same sequence of col / row
	fromMatr <- temp[, order(rownames(temp))]
	
	# obtain transition matrix
	outMatr <- fromMatr / rowSums(fromMatr)
	
	out <- new("markovchain", states = rownames(temp), 
	           transitionMatrix = outMatr, byrow=TRUE)
	
	return(out)
}

# coerce table to markovchain object
setAs(from = "table", to = "markovchain", def = .table2Mc)


# function from msm to markovchain
# msm is a package. Use this package to create msm object.
# see how to create msm object using ?msm

.msm2Mc <- function(from) {
  if(requireNamespace(package='msm', quietly = TRUE)) {
  temp <- msm::pmatrix.msm(from)
  prMatr <- unclass(as.matrix(temp))
  out <- new("markovchain", transitionMatrix = prMatr)
  } else {
    out <- NULL
    print("msm unavailable")
  }
  return(out)
}


# coerce msm object to markovchain object
setClass("msm")
setAs(from = "msm", to = "markovchain", def = .msm2Mc)


# function for msm.est to mc. Assume a probability matrix given
.msmest2Mc <- function(from) {
  
  if (is.matrix(from)) {
    # central estimate
    pMatr <- from 
  }
    
  if (is.list(from)) {
    # central estimate
    pMatr <- from[[1]] 
  }
    
  out <- new("markovchain", transitionMatrix = as(pMatr, "matrix"))
  
  return(out)
}

# coerce ms.est to markovchain object
setClass("msm.est")
setAs(from = "msm.est", to = "markovchain", def = .msmest2Mc)


# function from etm to markovchain
.etm2Mc<-function(from) {
  
  # data frame consists of  'from' and 'to' column
  df <- from$trans
  
  # name of states
  elements <- from$state.names
  # number of unique states
  nelements <- length(elements)
  
  # temporary t-matrix
  prMatr <- zeros(nelements)
  dimnames(prMatr) <- list(elements, elements)
  
  # populate t-matrix
  for(i in 1:dim(df)[1]) {
    r <- df[i, ] # each row one by one
    stateFrom <- r$from
    stateTo <- r$to
    prMatr[stateFrom, stateTo] <- prMatr[stateFrom, stateTo] + 1
  }
  
  # convert freq-matrix to trans-matrix
  rsums <- rowSums(prMatr)
  prMatr <- prMatr / rsums
  
  # take care of rows with all entries 0
  if(any(rsums == 0)) {
    indicesToBeSanitized <- which(rsums == 0)
    
    for(i in indicesToBeSanitized) {
      for(j in 1:nelements) {
        prMatr[i, j] <- 1 / nelements
      }
    }
  }
  
  # create markovchain object
  out <- new("markovchain", transitionMatrix = prMatr)
  return(out)
}

# coerce etm object to markovchain object
setClass("etm")
setAs(from = "etm", to = "markovchain", def = .etm2Mc)


#sparse matrix from Matrix package
.sparseMatrix2markovchain<-function(from){
  temp<-as(from,"matrix")
  out <- as(temp, "markovchain")
  return(out)
}

.markovchain2sparseMatrix<-function(from){
  temp<-as(from,"matrix")
  out <- as(temp, "sparseMatrix")
  return(out)
}


setAs(from = "sparseMatrix", to = "markovchain", def = .sparseMatrix2markovchain)
setAs(from = "markovchain", to = "sparseMatrix", def = .markovchain2sparseMatrix)



# functions and methods to return a matrix
.mc2matrix <- function(from) {
	out <- from@transitionMatrix
	return(out)
}

# coerce markovchain object to matrix(transition)
setAs(from = "markovchain", to = "matrix", def = .mc2matrix)


# functions and methods to return a matrix
.mc2igraph <- function(from) {
  
  # convert the markovchain to data.frame
	temp <- .mc2Df(from=from) 
	
	# convert the data frame to igraph graph
	# need to set only non zero weights
	out <- graph.data.frame(d=temp[temp$prob>0,]) 
	return(out)
}

# coerce markovchain object to igraph
setClass("igraph")
setAs(from = "markovchain", to = "igraph", def = .mc2igraph)


#' @exportMethod t
setGeneric("t")


# transposing method for markovchain objects
setMethod("t", "markovchain", 
  function(x) { 
    out <- new("markovchain", byrow = !x@byrow, 
               transitionMatrix = t(x@transitionMatrix))
    
    return(out)
  } 
)

#' @exportMethod *
setGeneric("*")

# function to multiplicate two markov chains
#
# Args:
# e1: first markovchain
# e2: second markov chain
#
# Returns:
# if feasible, a markovchain where the transition matrix is e1*e2

setMethod("*", c("markovchain", "markovchain"),
  function(e1, e2) {
    
    # compare states of markovchains
    if(!setequal(e1@states, e2@states)) {
      warning("Warning! Different states")
    }
    
    # dimension must be equal
    if(!setequal(dim(e1@transitionMatrix), dim(e2@transitionMatrix))) {
      stop("Error! Different size")
    }

    # both must be either row wise or col wise
    if(!(e1@byrow == e2@byrow)) {
      stop("Error! Both transition matrix should be defined either by row or by column")
    }

    newStates <- e1@states
    newTransMatr <- e1@transitionMatrix %*% e2@transitionMatrix
    byRow <- e1@byrow
    # multiplicated matrix takes the first matrix's name
    mcName <- e1@name 
    
    out<-new("markovchain", states = newStates, transitionMatrix = newTransMatr, 
             byrow = byRow, name = mcName)
    
    return(out)
  }
)

# methods implemented for multiplication of markovchain object with 
# matrix, 1-D vector, and vice-versa

setMethod("*", c("matrix", "markovchain"),
  function(e1, e2) {
    out <- e1 %*% e2@transitionMatrix
    return(out)
  }
)

setMethod("*", c("markovchain", "matrix"),
  function(e1, e2) {
    out <- e1@transitionMatrix %*% e2
    return(out)
  }
)

setMethod("*", c("numeric", "markovchain"), 
  function(e1, e2) {
    if(length(e1) != dim(e2)) {
      stop("Error! Uncompatible dimensions")
    } else {
      out <- e1 %*% e2@transitionMatrix
    }
    
    return(out)
  }
)

setMethod("*", c("markovchain", "numeric"), 
  function(e1, e2) {
    if(length(e2) != dim(e1)) {
      stop("Error! Uncompatible dimensions")
    } else {
      out <- e1@transitionMatrix %*% e2
    }
    
     return(out)
  }
)

#' @exportMethod ==
setGeneric("==")

# compare two markovchain object
setMethod("==", c("markovchain", "markovchain"), 
  function(e1, e2) {
    out <- .approxEqualMatricesRcpp(e1@transitionMatrix, e2@transitionMatrix)
    return(out)
  }
)

#' @exportMethod !=
setGeneric("!=")

setMethod("!=", c("markovchain", "markovchain"), 
  function(e1, e2) {
    out <- FALSE
    out <- !(e1 == e2)
    return(out)
  }
)

#'@exportMethod ^
setGeneric("^")

# markovchain raise to some power
# this method is O(n³ log(m)) where n = {num cols (= rows) of e1} and m = e2
setMethod("^", c("markovchain", "numeric"), 
  function(e1, e2) {
    out <- new("markovchain", states = e1@states, byrow = e1@byrow,
               transitionMatrix = e1@transitionMatrix %^% e2,
               name = paste(e1@name, "^", e2, sep = "")
              )
    
    return(out)
  }
)

#' @exportMethod [
setGeneric("[")

# methods to directly access transition matrix elements
setMethod("[", signature(x = "markovchain", i = "ANY", j = "ANY"), 
  function(x, i, j) {
    out <- x@transitionMatrix[i, j]
    return(out)
  }
)

#' @exportMethod [[
setGeneric("[[")

# methods to directly access markovchain objects composing a markovchainList object
setMethod("[[", signature(x = "markovchainList", i = "ANY"),
  function(x, i) {
    out <- x@markovchains[[i]]
    return(out)
  }
)

# transition probabilty vector from a given state

#' @title \code{conditionalDistribution} of a Markov Chain
#' 
#' @name conditionalDistribution
#' 
#' @description It extracts the conditional distribution of the subsequent state, 
#'              given current state.
#' 
#' @param object A \code{markovchain} object.
#' @param state Subsequent state.
#' 
#' @author Giorgio Spedicato, Deepak Yadav
#' 
#' @return A named probability vector
#' @references A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
#' 
#' @seealso \code{\linkS4class{markovchain}}
#' 
#' @examples 
#' # define a markov chain
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix = 
#'                matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1),nrow = 3, 
#'                       byrow = TRUE, dimnames = list(statesNames, statesNames)))
#'                       
#' conditionalDistribution(markovB, "b")                       
#' 
#' @exportMethod conditionalDistribution
setGeneric("conditionalDistribution", function(object, state) standardGeneric("conditionalDistribution"))
setMethod("conditionalDistribution", "markovchain",
  function(object, state) {
    # get the states names
    stateNames <- states(object) 
    
    # number of unique states
    out <- numeric(length(stateNames))
    
    # states are assumed to be sorted
    index2Take <- which(stateNames == state) 
    
    if(object@byrow == TRUE) {
      out <- object@transitionMatrix[index2Take, ]
    } else {
      out <- object@transitionMatrix[, index2Take]
    }
    
    # names the output and returs it
    names(out) <- stateNames
    
    return(out) 
  }
)
		  
# Function to get the mode of a probability vector
# 
# Args:
# probVector: the probability vector
# ties: specifies if ties are to be sampled, otherwise more than one element is returned
#
# Returns:
# the name of the model element

.getMode <- function(probVector, ties = "random") {
	maxIndex <- which(probVector == max(probVector))
	temp <- probVector[maxIndex] # index of maximum probabilty
	
	if((ties == "random") & (length(temp) > 1)) {
	  out <- sample(temp, 1) 
	} else {
	  out <- temp
	}
	
	return(names(out))
}


#' @exportMethod predict
setGeneric("predict")

# predict method for markovchain objects
# given initial state return a vector of next n.ahead states

setMethod("predict", "markovchain", 
  function(object, newdata, n.ahead = 1) {
    # identify the last state
    lastState <- newdata[length(newdata)]
    out <- character()
    
    for(i in 1:n.ahead) {
      # cyclically determine the most probable subsequent state from the conditional distribution
      newState <- .getMode(probVector = conditionalDistribution(object, lastState), ties = "random") 
      out <- c(out, newState)
      lastState <- newState
    }
    
    return(out)
  }
)

# predict method for markovchainList objects
setMethod("predict", "markovchainList",
  function(object, newdata, n.ahead = 1, continue = FALSE) {
    # object a markovchainList
    # newdata = the actual data 
    # n.ahead = how much ahead 
    # continue = veryfy if that lasts
      
    # allocate output
    out <- character() 
    actualPos <- length(newdata) 
    lastState <- newdata[actualPos] # take last position
    for(i in 1:n.ahead) {
      newPos <- actualPos + i - 1
      if(newPos <= dim(object)) {
        newState <- predict(object = object[[newPos]], newdata = lastState, n.ahead = 1)
        out <- c(out, newState)
        lastState <- newState
      } else {
          if(continue == TRUE) {
            newState <- predict(object = object[[dim(object)]], newdata = lastState, n.ahead = 1)
            out <- c(out, newState)
            lastState <- newState
          } else break;
        }
    }

    return(out)
  }
)

#sort method for markovchain objects

setGeneric("sort", function(x, decreasing=FALSE, ...) standardGeneric("sort"))

setMethod("sort", signature(x="markovchain"), function(x, decreasing = FALSE) {
    #get matrix and state names 2 be sorted
   
    matr2besorted<-x@transitionMatrix 
    if (x@byrow) 
      states2besorted <- rownames(matr2besorted) 
    else 
      states2besorted <- colnames(matr2besorted)
    
    #sorting
    sort_index<-order(states2besorted,decreasing = decreasing)
    
    #reallocating
    matr_sorted<-matr2besorted[sort_index,sort_index]
    states_sorted<-states2besorted[sort_index]
    
    out<-x
    
    out@transitionMatrix<-matr_sorted
    out@states<-states_sorted
    
    return(out)
  }
)


# method to get stationary states

#' @name steadyStates
#' @title Stationary states of a \code{markovchain} object
#' 
#' @description This method returns the stationary vector in matricial form of a markovchain object.
#' @param object A discrete \code{markovchain} object
#' 
#' @return A matrix corresponding to the stationary states
#' 
#' @references A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
#' @author Giorgio Spedicato
#' @seealso \code{\linkS4class{markovchain}}
#' 
#' @note The steady states are identified starting from which eigenvectors correspond 
#'       to identity eigenvalues and then normalizing them to sum up to unity. When negative values are found 
#'       in the matrix, the eigenvalues extraction is performed on the recurrent classes submatrix.
#'       
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#'                 matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
#'                 byrow = TRUE, dimnames=list(statesNames,statesNames)),
#'                name = "A markovchain Object" 
#' )       
#' steadyStates(markovB)
#' 
#' @rdname steadyStates
#' @exportMethod steadyStates
setGeneric("steadyStates", function(object) standardGeneric("steadyStates"))
