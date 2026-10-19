#' @title Continuous time Markov chain class
#' @name ctmc-class
#' @aliases dim,ctmc-method
#' @aliases initialize,ctmc_method
#' @aliases states,ctmc-method
#' @aliases steadyStates,ctmc-method
#' @aliases plot,ctmc,missing-method
#' @description The S4 class of continuous time Markov chains, described by
#'   their generator (rate) matrix.
#' 
#' @slot states Names of the states, equal to the row and column names of the
#'   generator.
#' @slot byrow Logical: are the rates out of each state in its row (rows sum
#'   to zero) or in its column?
#' @slot generator The square generator matrix: non-negative off-diagonal
#'   rates, with rows (or columns) summing to zero.
#' @slot name Optional name of the chain.
#' 
#' @section Methods:
#' 
#' \describe{
#' \item{dim}{\code{signature(x = "ctmc")}: method to get the size} 
#' \item{initialize}{\code{signature(.Object = "ctmc")}: initialize
#'   method }
#' \item{states}{\code{signature(object = "ctmc")}: states method. }
#' \item{steadyStates}{\code{signature(object = "ctmc")}: method to get the
#'   steady state vector. } 
#' \item{plot}{\code{signature(x = "ctmc", y = "missing")}: plot method 
#'   for \code{ctmc} objects }
#' }
#' 
#' @references
#' Introduction to Stochastic Processes with Applications in the Biosciences
#' (2013), David F. Anderson, University of Wisconsin at Madison. Sai Bhargav
#' Yalamanchi, Giorgio Spedicato
#' 
#' @note The validity method checks that the rows (or columns) of the
#'   generator sum to zero, rounded to the fifth decimal. If the generator has
#'   no state names, the states are named "s1", "s2", ...
#' @seealso \code{\link{generatorToTransitionMatrix}},\code{\link{rctmc}}
#'
#' @examples
#' energyStates <- c("sigma", "sigma_star")
#' byRow <- TRUE
#' gen <- matrix(data = c(-3, 3,
#'                        1, -1), nrow = 2,
#'               byrow = byRow, dimnames = list(energyStates, energyStates))
#' molecularCTMC <- new("ctmc", states = energyStates, 
#'                      byrow = byRow, generator = gen, 
#'                      name = "Molecular Transition Model")
#' steadyStates(molecularCTMC)
#' \donttest{plot(molecularCTMC)}
#' 
#' @keywords classes
#' 
#' @export
# S4 class 'ctmc' = continuous time Markov chain. It is described by a GENERATOR matrix Q
# (off-diagonal entries = jump rates, diagonal = minus the total rate of leaving the state,
# so every row sums to zero), instead of a matrix of probabilities.
setClass("ctmc",
         slots = list(states = "character",
                      byrow = "logical",
                      generator = "matrix",
                      name = "character"),
         prototype = list(
           states = c("a", "b"),
           byrow = TRUE,
           generator = matrix(data = c(-1, 1, 1, -1), byrow = TRUE, nrow = 2,
                              dimnames = list(c("a", "b"), c("a", "b"))),
           name = "Unnamed CTMC")
)

# initialize() is called by new(): it fills in default values for the missing arguments and
# repairs the dimension names before the object is built.
setMethod("initialize",
          signature(.Object = "ctmc"),
          function (.Object, states, byrow, generator,name,...) {
            # put the standard markovchain 
            if(missing(generator)) generator=matrix(data=c(-1, 1, 1, -1), #create a naive matrix
                                                    nrow=2,
                                                    byrow=TRUE, 
                                                    dimnames=list(c("a", "b"), c("a", "b"))
            )
            
            # check names of transition matrix, filling in missing/inconsistent
            # row or column names and resolving the states vector accordingly
            # Shared helper (utils.R): makes row names, column names and states consistent.
            filled <- .fillDimNames(generator, if (missing(states)) NULL else states)
            generator <- filled$matrix
            states <- filled$states
            if(missing(byrow)) byrow=TRUE #set byrow as true by default
            if(missing(name)) name="Unnamed Markov chain"  #generic name to the object
            callNextMethod(.Object, states = states, byrow = byrow, generator=generator,name=name,...)
          }
)

# Accessors: states(x) returns the state names, dim(x) the number of states.
#returns states of the ctmc
setMethod("states", "ctmc", 
          function(object) {
            out <- object@states
            return(out)
          }
)

#returns states of the ctmc
setMethod("dim", "ctmc", 
          function(x) {
            out <- nrow(x@generator)
            return(out)
          }
)

# Validity: runs every time an object is created and returns TRUE or an error message.
# Here it checks that the matrix is a generator (rows, or columns, sum to 0), is square and
# that its names match the states.
setValidity("ctmc",
            function(object) {
              check<-NULL
              # performs a set of check whose results are saved in check
              if (!.isGenRcpp(object@generator)) check <- "Error! Not a generator matrix" 
              if (object@byrow) {
                if(any(round(rowSums(object@generator),5)!=0)) check <- "Error! Row sums not equal to zero"
              } else {
                if(any(round(colSums(object@generator),5)!=0)) check <- "Error! Col sums not equal to zero"
              } 
              if (nrow(object@generator)!=ncol(object@generator)) check <- "Error! Not squared matrix" #check if squalre matrix
              if (!setequal(colnames(object@generator),object@states)) check <- "Error! Colnames <> states" #checks if 
              if (!setequal(rownames(object@generator),object@states)) check <- "Error! Rownames <> states"
              if ( is.null(check) ) return(TRUE) else return(check)
            }
)

# Finds the eigenvector of the eigenvalue 1 of a matrix. For a stochastic matrix this eigenvector is
# proportional to the stationary distribution. Returns NULL (with a warning) if there is not exactly
# one such eigenvalue, or if the vector is complex.
.ctmcEigen<-function(matr, transpose=TRUE)
{
  # Function to extract eigenvalues, core of get steady states 
  #
  # Args:
  # matr: the matrix to extract
  # transpose:  boolean indicating whether the matrx shall be transpose
  #
  # Results:
  # a matrix / vector
  # Row-stochastic matrices need the transpose: the stationary distribution is a LEFT eigenvector.
  if (transpose) tMatr <- t(matr) else tMatr <- matr #trasposing
  eigenResults <- eigen(x=tMatr,symmetric=FALSE) #perform the eigenvalue extraction
  # Pick the eigenvalues equal to 1 (rounded to 3 digits to absorb numerical noise).
  onesIndex <- which(round(eigenResults$values,3)==1) #takes the one eigenvalue
  #do the following: 1:get eigenvectors whose eigenvalues==1
  #2: normalize
  if (length(onesIndex)==0) {
    warning("No eigenvalue = 1 found - the embedded Markov Chain must be irreducible, recurrent")
    return(NULL)
  }
  if(length(onesIndex) > 1){
    warning("Eigenvalue = 1 multiplicity > 1! - the embedded Markov Chain must be irreducible, recurrent")
    return(NULL)
  }
  if (transpose)
  {
    eigenTake <- as.matrix(t(eigenResults$vectors[,onesIndex])) 
    if(rowSums(Im(eigenTake)) != 0){
      warning("Eigenvector corresponding to largest eigenvalue has a non-zero imaginary part - the embedded Markov Chain must be irreducible, recurrent")
      return(NULL)
    }
    out <- eigenTake
  } else {
    eigenTake <- as.matrix(eigenResults$vectors[,onesIndex]) 
    if(colSums(Im(eigenTake)) != 0){
      warning("Eigenvector corresponding to largest eigenvalue has a non-zero imaginary part - the embedded Markov Chain must be irreducible, recurrent")
      return(NULL)
    }
    out <- eigenTake
  }
  return(out)
}

# Steady state of a ctmc: (1) take the embedded jump chain (the discrete chain that records only the
# visited states), (2) find its stationary vector, (3) divide each entry by the rate of leaving the
# state (the longer a chain stays in a state, the more weight that state has) and renormalise.
setMethod("steadyStates", "ctmc", 
          function(object) {
            transposeYN <- object@byrow
            transMatr <- generatorToTransitionMatrix(object@generator, byrow = object@byrow)
            out<-.ctmcEigen(matr=transMatr, transpose=transposeYN) 
            if(is.null(out)) {
              warning("Warning! No steady state")
              return(NULL)
            }
            if (transposeYN) { 
              colnames(out) <- object@states
            } else {
              rownames(out) <- object@states
            }
            out <- - out / diag(object@generator)
            if (transposeYN) {
              out <- out / rowSums(out)
            }
            else{
              out <- out / colSums(out)
            }
            return(out)
          }
)



# Builds an igraph graph from the generator, for plotting (weights are multiplied by 100 to keep
# the digits when igraph rounds them).
# internal function for plotting ctmc object using igraph
.getNetctmc <- function(object, round = FALSE) {
  # function to get the graph adjacency object to plot and export to igraph
  #
  # Args: 
  # object: a ctmc object
  # round: boolean to round
  #
  # Returns:
  #
  # a graph adjacency
  
  # The graph has an edge i -> j for every positive rate. The generator is
  # taken by rows (transposed if stored by columns), and its diagonal, the
  # negative total exit rates, is left out: as self-loops with negative
  # weights it only made igraph warn and drop all the weights in the layout.
  generator <- if (object@byrow) object@generator else t(object@generator)
  diag(generator) <- 0
  matr <- generator * 100
  if (round) {
    matr <- round(matr, 2)
  }
  
  net <- graph.adjacency(adjmatrix = matr, weighted = TRUE, mode = "directed")
  return(net)
}


# plot(): the package argument chooses the drawing library; if the chosen one is not installed
# we fall back on igraph, which is always available.
setMethod("plot",signature(x="ctmc",y="missing"),
          function(x,y,package = "igraph",...){
            switch(package,
                   diagram = {
                     if (requireNamespace("diagram", quietly = TRUE)) {
                       .plotdiagram(object = x, ...)
                     } else {
                       netMc <- .getNetctmc(object = x, round = TRUE)
                       edgeLabel <- round(E(netMc)$weight / 100, 2)
                       plot.igraph(x = netMc, edge.label = edgeLabel, ...)
                     }
                   },
                   
                   DiagrammeR = {
                     if (requireNamespace("DiagrammeR", quietly = TRUE)) {
                       .plotDiagrammeR(object = x, ...)
                     } else {
                       netMc <- .getNetctmc(object = x, round = TRUE)
                       edgeLabel <- round(E(netMc)$weight / 100, 2)
                       plot.igraph(x = netMc, edge.label = edgeLabel, ...)
                     }
                   },
                   {
                     netMc <- .getNetctmc(object = x,round = TRUE)
                     edgeLabel <- round(E(netMc)$weight / 100, 2)
                     plot.igraph(x = netMc, edge.label = edgeLabel, ...)
                   })
          }
          )


# 'ictmc' = imprecise CTMC: a generator Q whose rates may vary within the interval given in `range`.

#' An S4 class for representing Imprecise Continuous Time Markovchains
#' 
#' @slot states a vector of states present in the ICTMC model
#' @slot Q matrix representing the generator demonstrated in the form of variables
#' @slot range a matrix that stores values of range of variables
#' @slot name name given to ICTMC
#' 
ictmc <- setClass("ictmc",
                  slots = list(states = "character", Q = "matrix",
                               range = "matrix", name = "character")
                  )




# Default values of the arguments, as for the ctmc class.
setMethod("initialize",
          signature(.Object = "ictmc"),
          function (.Object, states, Q, range, name, ...) {
            
            
            if(missing(Q)) Q=matrix(data=c(-1, 1, 1, -1), #create a naive matrix
                                                    nrow=2,
                                                    byrow=TRUE, 
                                                    dimnames=list(c("n", "y"), c("n", "y"))
            )
            
            
            if(missing(range)) range = matrix(c(1/52, 3/52, 1/2, 2),
                                              nrow = 2,
                                              byrow = 2)
            #if all names are missing it initializes them to "1", "2",...
            filled <- .fillDimNames(Q, if (missing(states)) NULL else states)
            Q <- filled$matrix
            states <- filled$states
            
            if(missing(name)) name="Unnamed imprecise CTMC"  #generic name to the object
            callNextMethod(.Object, states = states, Q = Q, range=range,name=name,...)
          }
          
)




# Validity of an ictmc: Q must be a generator and each interval [low, high] must be well ordered
# and not negative.
setValidity("ictmc",
            function(object) {
              check<-NULL
              # performs a set of check whose results are saved in check
              if (!.isGenRcpp(object@Q)) check <- "Error! Not a generator matrix" 
              if(any(round(rowSums(object@Q),5)!=0)) check <- "Error! Row sums not equal to zero"
              if ( nrow(object@Q) != ncol(object@Q )) check <- "Error! Not squared matrix" #check if square matrix
              if ( !setequal(colnames(object@Q),object@states )) check <- "Error! Colnames <> states" 
              if ( !setequal(rownames(object@Q),object@states )) check <- "Error! Rownames <> states"
              
              
              if(nrow(object@range) != nrow(object@Q) && ncol(object@range) != 2) check <- "Error! dimension of range matrix not correct."
              for(i in 1:nrow(object@Q)){
                if( object@range[i,1] > object@range[i,2] ){
                  check <- "Error, improper values set in range matrix."
                }
                  
              }
              if(min(object@range) < 0){
                check <- "Error, values in the range matrix should be greater than zero."
              }
              
              
              if ( is.null(check) ) return(TRUE) else return(check)
            }
)






