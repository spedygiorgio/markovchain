# define Markov Chain class
setClass("markovchain", # class name
         
         # Define the slots
         slots = list(states = "character", byrow = "logical",
                      transitionMatrix = "matrix", name = "character"),
         
         # Set the default values for the slots
         prototype = list(states = c("a","b"), byrow = TRUE, 
                          transitionMatrix = matrix(data = c(0,1,1,0),
                                                    nrow = 2, byrow = TRUE, 
                                                    dimnames = list(c("a","b"), c("a","b"))),  
                          name = "Unnamed Markov chain")
)

# initializing method for markovchain objects
setMethod("initialize",
          signature(.Object = "markovchain"),
          function (.Object, states, byrow, transitionMatrix, name, ...) {
            
            # put the standard markovchain 
            if(missing(transitionMatrix)) {
              transitionMatrix <- matrix(data = c(0, 1, 1, 0),
                                         nrow = 2,
                                         byrow = TRUE, 
                                         dimnames = list(c("a","b"), c("a","b"))
              ) 
            }
            
            # check names of transition matrix
            # if all names are missing it initializes them to "1", "2", ....
            
            if(all(is.null(rownames(transitionMatrix)), is.null(colnames(transitionMatrix))) == TRUE) { 
              
              if(missing(states)) {
                nr <- nrow(transitionMatrix)
                stateNames <- as.character(seq(1:nr))
              } else {stateNames <- states}
              
              rownames(transitionMatrix) <- stateNames
              colnames(transitionMatrix) <- stateNames
              
            } else if(is.null(rownames(transitionMatrix))) { # fix when rownames null
              
              rownames(transitionMatrix) <- colnames(transitionMatrix)
              
            } else if(is.null(colnames(transitionMatrix))) { # fix when colnames null
              
              colnames(transitionMatrix) <- rownames(transitionMatrix)
              
            } else if(!setequal(rownames(transitionMatrix), colnames(transitionMatrix)))  {
              
              colnames(transitionMatrix) <- rownames(transitionMatrix) # fix when different
              
            }
            
            if(missing(states)) {
              states <- rownames(transitionMatrix)
            }
            
            if(missing(byrow)) {
              byrow <- TRUE
            }
            
            if(missing(name)) {
              name <- "Unnamed Markov chain" 
            }
            
            callNextMethod(.Object, states = states, byrow = byrow, 
                           transitionMatrix = transitionMatrix, name = name, ...)
          }
)

# define Markov Chain List class
setClass("markovchainList", 
         slots = list(markovchains = "list", 
		                  name = "character")
)

# verifies if a markovchainList object is valid or not
setValidity("markovchainList",
		         function(object) {
		           check <- FALSE 
		           for(i in 1:length(object@markovchains)) {
			            if(class(object@markovchains[[i]]) != "markovchain") {
			              # All elements in the list should be a markovchain object 
			              check <- "Error! All elements should be of class 'markovchain'" 
			            }
		           }
		           
		           if(check == FALSE) check <- TRUE
		           return(check)
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
#'                name = "A markovchain Object" 
#' )
#' states(markovB)
#' 
#' @rdname states
#' 
#' @export
setGeneric("states", function(object) standardGeneric("states"))

#' @rdname states
setMethod("states","markovchain", 
          function(object) {
            out <- object@states
            return(out)
          }
)

#' @title Method to retrieve name of markovchain object  
#' 
#' @name name
#' 
#' @description This method returns the name of markovchain object
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
#'                name = "A markovchain Object" 
#' )
#' name(markovB)
#' 
#' @export
setGeneric("name", function(object) standardGeneric("name"))


#' @rdname getName
setMethod("name", "markovchain", function(object) {
  out <- object@name
  return(out)
})

#' @title Method to set name of markovchain object
#' 
#' @name name<-
#' 
#' @description This method modify the existing name of markovchain object
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
#'                name = "A markovchain Object" 
#' )
#' name(markovB) <- "dangerous mc"
#' 
#' @export
setGeneric("name<-", function(object, value) standardGeneric("name<-"))

#' @rdname setName
setMethod("name<-", "markovchain", 
          function(object, value) {
            object@name <- value
            object
          }
)

# adding a method names: to get names
setMethod("names","markovchain", 
          function(x) {
            out <- x@states
            return(out)
          }
)

# adding a method names: to set names
setMethod("names<-", "markovchain", 
          function(x, value) {
            rownames(x@transitionMatrix) <- value
            colnames(x@transitionMatrix) <- value
            x@states <- value
            return(x)
          }
)


# generic methods to get the dim of a markovchain and markovchainList

setMethod("dim", "markovchain", 
		function(x) {
			out <- nrow(x@transitionMatrix)
			return(out)
		}
)

setMethod("dim", "markovchainList", 
		function(x) {
			  out <- length(x@markovchains)
			  return(out)
		  }
)


# method  to set the validity of a markovchain object
setValidity("markovchain",
		function(object) {
			check <- FALSE
			
			# performs a set of check whose results are saved in check
			
			# check all values of transition matrix belongs to [0, 1]
			if (any(sapply(as.numeric(object@transitionMatrix),.isProbRcpp)) == FALSE) {
			  check <- "Error! Some elements are not probabilities"
			}
			
			# rows sum or columns sum = 1
			if (object@byrow == TRUE) {
			  
			  # absolute difference
			  absdiff <- abs(1-zapsmall(rowSums(object@transitionMatrix)))
				
			  if(any(absdiff > .Machine$double.eps*100)) {
			    pos2check<-which(absdiff > .Machine$double.eps*100)
				  check <- paste("Error! Row sums not equal to one","check positions:",pos2check)
				}
			} else {
			  
			  # absolute difference
			  absdiff <- abs(1-zapsmall(colSums(object@transitionMatrix)))
			  
			  if(any(absdiff > .Machine$double.eps*100)) {
			    pos2check<-which(absdiff > .Machine$double.eps*100)
			    check <- paste("Error! Row sums not equal to one","check positions:",pos2check)
				}
			}
			
			
			# check whether matrix is square amtrix or not
			if (nrow(object@transitionMatrix) != ncol(object@transitionMatrix)) {
			  check <- "Error! Not squared matrix" #check if squalre matrix
			}
      
			# check whether column names or rows names equal to state names or not
			if (!setequal(colnames(object@transitionMatrix), object@states)) {
			  check <- "Error! Colnames <> states" 
			}
      if (!setequal(rownames(object@transitionMatrix), object@states)) {
        check <- "Error! Rownames <> states"
      }
			
      if (check == FALSE) {
        return(TRUE)
      }  else {
        return(check)
      }
		}
)

# matr : matrix
# transpose : boolean indicating whether the matrix shall be transposed or not
# output : a matrix / vector

.mcEigen <- function(matr, transpose = TRUE) {
  
  if (transpose) {
    tMatr <- t(matr) 
  } else {
    tMatr <- matr
  }
  
  # perform the eigenvalue extraction
  eigenResults <- eigen(x = tMatr, symmetric = FALSE) 
  
  # takes the one eigenvalues
  onesIndex <- which( sapply (eigenResults$values, function(e){ isTRUE(all.equal( as.complex(e),1+0i)) } ))
  
  # do the following: 
  # 1 : get eigenvectors whose eigenvalues == 1
  # 2 : normalize
  
  if (length(onesIndex) == 0) {
    warning("No eigenvalue = 1 found")
    return(NULL)
  }
  
  # Gives always a norm-based order to eigenvectors
  eigenvectors <- as.matrix( eigenResults$vectors[, onesIndex] )

  if (transpose == TRUE) {
    eigenTake <- as.matrix(t(eigenvectors)) 
    out <- eigenTake / rowSums(eigenTake) # normalize
  } else {
    eigenTake <- as.matrix(eigenvectors) 
    out <- eigenTake / colSums(eigenTake) # normalize
  }
  
  # subset the eigenvectors
  # normalize
  # take the real part: need to be sanitized
	# @DEEPAK: later we have to see and optimize this part. I am not sure taking
	#       the real part is most appropriate.
  
  out <- Re(out)
  return(out)
}

# method to get stationary states

#' @name steadyStates
#' @title Stationary states of a \code{markovchain} objeect
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
#' @export
setGeneric("steadyStates", function(object) standardGeneric("steadyStates"))

#' @rdname steadyStates
setMethod("steadyStates","markovchain", 
		function(object) {
			transposeYN <- FALSE
			if(object@byrow == TRUE) {
			  transposeYN <- TRUE		
			}
      
			out <- .mcEigen(matr = object@transitionMatrix, transpose = transposeYN)
			
			# if any element negative
			if (min(out)<0) {
			  warning("Negative elements in steady states, working on closed classes submatrix")
			  if(object@byrow==TRUE) myObject=object else myObject=t(object)
			  out <- .steadyStatesByRecurrentClasses(object=myObject)
			}
      
			if(is.null(out)) {
				warning("Warning! No steady state")
				return(NULL)
			} else{
			  # order vectors lexicographically
			  out <- .mcLexSort(out)
		    if (object@byrow==FALSE) out<-t(out)
		  }
			
			if(transposeYN == TRUE) { 
				colnames(out) <- object@states
			} else {
				rownames(out) <- object@states
			}

      return(out)
    }
)




#' @title steadyStatesByRecurrent classes
#' 
#' @description Function to extract steady states when needed
#' Only recurrent closed classes are considered 
#' @author Christope Dutang and Giorgio Spedicato
#' @return A matrix
.steadyStatesByRecurrentClasses<-function(object) {
  #inizialization
  M<-object@transitionMatrix
  #transpose bycol matrices
  if (object@byrow==FALSE) M <- t(M)
  namesSequence<-names(object)
  #characterizin recurrent classes
  recClasses<-recurrentClasses(object)
  numRecClasses<-length(recClasses)
  recurrentClassesNames<-unlist(recClasses)
  #extracting recurrent classes
  Msub <- M[rownames(M) %in% recurrentClassesNames, colnames(M) %in% recurrentClassesNames]
  
  out<-matrix(0, nrow=numRecClasses, ncol = dim(object))
  colnames(out)<-names(object)
  #getting their steady states, calculating first indexes of eigenvalues equal to 1
  onesIndex <- which( sapply (eigen(Msub)$values, function(e){ isTRUE(all.equal( as.complex(e), 1+0i)) }) )
  partialOutput<-t(eigen(Msub)$vectors[, onesIndex]) / colSums(eigen(Msub)$vectors[, onesIndex])
  colnames(partialOutput)<-recurrentClassesNames
  #allocating to their columns
  out[,colnames(out) %in% recurrentClassesNames]<-partialOutput
  return(out)
}


# generic function to extract absorbing states

#' @rdname absorbingStates
#' 
#' @export
setGeneric("absorbingStates", function(object) standardGeneric("absorbingStates"))
setMethod("absorbingStates", "markovchain", 
          function(object) {
			      out <- character()
			      matr <- object@transitionMatrix
			      transposeYN <- FALSE
			      if(object@byrow == TRUE) {
			        transposeYN <- TRUE
			      }
			      
			      steadyStates <- .mcEigen(matr = matr, transpose = transposeYN)
			      if(is.null(steadyStates)) {
			        return(character(0))
			      }
			      
			      # identify which states are absorbing if they are diagonally one
			      if(transposeYN == TRUE) {
			        maxCols <- apply(steadyStates, 2, "max")
			      } else {
			        maxCols <- apply(steadyStates, 1, "max")  
			      }
			      
			      index <- which(maxCols == 1)
			      if(length(index) > 0) {
			        out <- object@states[index] 
			      }
			      
			      return(out)
          } 
)

# generic method to extract transient states

#' @rdname absorbingStates
#' 
#' @export
setGeneric("transientStates", function(object) standardGeneric("transientStates"))

#' @rdname absorbingStates
setMethod("transientStates", "markovchain", 
	       	function(object) {
			      out <- character()
			      
			      # make byRow = true for the matrix
			      if(object@byrow == TRUE) {
			        matr <- object@transitionMatrix
			      } else {
			        matr <- t(object@transitionMatrix)
			      }
			      
			      temp <- .commclassesKernelRcpp(matr)
			      index <- which(temp$v == FALSE)
			      if(length(index) > 0) {
			        out <- names(temp$v[index])
			      }
			      
			      return(out)
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
#' @export
setGeneric("transitionProbability", function(object, t0, t1) standardGeneric("transitionProbability"))

#' @rdname transitionProbability
setMethod("transitionProbability", "markovchain", 
	        function(object, t0, t1) {
		        out <- numeric(1)
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



# method to convert into canonic form : a markovchain object
# TODO: check meaning of this function

#' @rdname absorbingStates
#' 
#' @export
setGeneric("canonicForm", function(object) standardGeneric("canonicForm"))
setMethod("canonicForm", "markovchain",
          function(object) {
            # Obtain the canonical form Q of a stochastic matrix P
            P <- object@transitionMatrix
            
            # Uses the internal function commclassesKernelRcpp
            comclasList <- .commclassesKernelRcpp(P)
            
            # vu is a row vector of 0s and 1s. vu(i) = 1 if
            # the class C(i) is closed, and 0 otherwise
            vu <- comclasList$v
			
            # find index of closed communicating classes
            u <- matlab::find(vu == TRUE)
            
            # find index of open communicating classes
            w <- matlab::find(vu == FALSE)
            
            # Cmatr(i,j) is 1 if and only if j is in the
            # communicating class of i.
            Cmatr <- comclasList$C
            
            # R is now the set of representatives of closed classes
            # Each closed class has a unique representative in R.
            R <- numeric()
            while(length(u) > 0) {
              # everytime add a unique closed communicating classes index
              R <- c(R, u[1])
              
              # remove the duplicate communication classes as u[1]
              vu <- as.logical(vu * (Cmatr[u[1], ] == FALSE))
              
              # rest communicating classes index are hidden inside u
              u <- find(vu == TRUE);
            }
            
            # we have now a permutation p of indices, p, that
            # gives the new stochastic matrix Q.
            p <- numeric()
            for (i in 1:length(R))
            {
              a <- find(Cmatr[R[i], ])
              p <- c(p,a)
            }
            
            # append open communicating classes index
            p <- c(p, w)
            
            # extract canonical form out of given matrix using 
            # permutation of indexes calculated above
            Q <- P[p, p]
            
            out <- new("markovchain", transitionMatrix = Q, name = object@name)
            return(out)
          }
)

.canonicForm <- function(object) {
  # Obtain the canonical form Q of a stochastic matrix P
  P <- object@transitionMatrix
  
  # Uses the internal function commclassesKernelRcpp
  comclasList <- .commclassesKernelRcpp(P)
  
  # vu is a row vector of 0s and 1s. vu(i) = 1 if
  # the class C(i) is closed, and 0 otherwise
  vu <- comclasList$v
  
  # find index of closed communicating classes
  u <- matlab::find(vu == TRUE)
  
  # find index of open communicating classes
  w <- matlab::find(vu == FALSE)
  
  # Cmatr(i,j) is 1 if and only if j is in the
  # communicating class of i.
  Cmatr <- comclasList$C
  
  # R is now the set of representatives of closed classes
  # Each closed class has a unique representative in R.
  R <- numeric()
  while(length(u) > 0) {
    # everytime add a unique closed communicating classes index
    R <- c(R, u[1])
    
    # remove the duplicate communication classes as u[1]
    vu <- as.logical(vu * (Cmatr[u[1], ] == FALSE))
    
    # rest communicating classes index are hidden inside u
    u <- find(vu == TRUE);
  }
  
  # we have now a permutation p of indices, p, that
  # gives the new stochastic matrix Q.
  p <- numeric()
  for (i in 1:length(R))
  {
    a <- find(Cmatr[R[i], ])
    p <- c(p,a)
  }
  
  # append open communicating classes index
  p <- c(p, w)
  
  # extract canonical form out of given matrix using 
  # permutation of indexes calculated above
  Q <- P[p, p]
  
  out <- new("markovchain", transitionMatrix = Q, name = object@name)
  return(out)
}

# summary method for markovchain class
# lists: closed, transient classes, irreducibility, absorbint, transient states
setMethod("summary", signature(object = "markovchain"),
		      function(object){
			      
		        # list of closed, recurrent and transient classes
		        outs <- .summaryKernelRcpp(object)
			      
		        # display name of the markovchain object
			      cat(object@name," Markov chain that is composed by:", "\n")
			      
			      # number of closed classes
			      check <- length(outs$closedClasses)
			      
			      cat("Closed classes:","\n")
			      
			      # display closed classes
			      if(check == 0) cat("NONE", "\n") else {
				      for(i in 1:check) cat(outs$closedClasses[[i]], "\n")
			      }
			      
			      # number of recurrent classes
			      check <- length(outs$recurrentClasses)
			
			      cat("Recurrent classes:", "\n")
			      
			      # display recurrent classes
			      if(check == 0) cat("NONE", "\n") else {
			          cat("{")
			          cat(outs$recurrentClasses[[1]], sep = ",")
			          cat("}")
			          if(check > 1) {
			            for(i in 2:check) {
			              cat(",{")
			              cat(outs$recurrentClasses[[i]], sep = ",")
			              cat("}")
			            }
			          }
			          cat("\n")
			      }
			      
			      # number of transient classes
			      check <- length(outs$transientClasses)
			      
			      cat("Transient classes:","\n")
			
			      # display transient classes
			      if(check == 0) cat("NONE", "\n") else {
			          cat("{")
			          cat(outs$transientClasses[[1]], sep = ",")
			          cat("}")
			          if(check > 1) { 
			            for(i in 2:check) {
			              cat(",{")
			              cat(outs$transientClasses[[i]], sep = ",")
			              cat("}")
			            }
			          }
			          cat("\n")
			      }
			
			      # bool to say about irreducibility of markovchain
			      irreducibility <- is.irreducible(object)
			      
			      if(irreducibility) 
			        cat("The Markov chain is irreducible", "\n") 
			      else cat("The Markov chain is not irreducible", "\n")
			      
			      # display absorbing states
			      check <- absorbingStates(object)
			      if(length(check) == 0) check <- "NONE"
			      cat("The absorbing states are:", check )
			      cat("\n")
			      
			      # return outs
			      # useful when user will assign the value returned
			      invisible(outs) 
          }
)

##################################################AS METHODS#########################

.checkMatrix <- function(matr, byrow = TRUE, verbose = FALSE) {
	
  # first check: size
	if (dim(matr)[1] != dim(matr)[2]) {
		if(verbose) stop("Error! Rectangular matrix")
		return(FALSE)
	}
	
	# second check: all elements are probs
	for(i in 1:nrow(matr)) {
		for(j in 1:ncol(matr)){
			if(!(.isProbRcpp(matr[i, j]))) {
				if(verbose) stop("Error! Non probabilities")
				  return(FALSE)
			}
		}
	}
	
	# third check: either columns or rows sum to one
  # to perform only one check 	
	if(byrow == FALSE) {
	  matr <- t(matr) 
	}
	
  # calculate row's sum
	check <- rowSums(matr)
	
	for( i in 1:length(check)) {
	  if (abs(1-check[i]) > .Machine$double.eps) {
	    if(verbose) {
	      myMessage<-paste("Error! Either rows or cols should sum to 1","check state",i)
	      stop(myMessage) 
	    }
	    return(FALSE)
	  }
	}
	
	# if all test are passed
	return(TRUE)
}

# Internal function to return a markovchain object given a matrix
.matrix2Mc <- function(from) {
	
  # whether given matrix is a transition matrix or not
  # if it is then how probabilities are stored
  # row-wise or columnwise
  
	byrow <- FALSE
	checkByRow <- .checkMatrix(from, byrow = TRUE)
	
	if(checkByRow) {
	  byrow <- TRUE
	} else  {
		checkByCols <- .checkMatrix(from, byrow = FALSE)
		if(!checkByCols) {
		  #error could be either in rows or in cols
		  if (any(colSums(from)!=1)) cat("columns sums not equal to one are:",which(colSums(from)!=1),"\n")
		  if (any(rowSums(from)!=1)) cat("columns sums not equal to one are:",which(rowSums(from)!=1),"\n")
		  stop("Error! Not a probability matrix")	
		}
	}
	
	# extract states names
	if(byrow == TRUE) {
	  namesCandidate <- rownames(from) 
	} else {
	  namesCandidate <- colnames(from)
	}
	
	# if states names is not there create it s1, s2, s3, ....
	if(is.null(namesCandidate)) {
		namesCandidate <- character()
		for(i in 1:nrow(from)) {
		  namesCandidate <- c(namesCandidate, paste("s", i, sep = "")) 
		}
	}
	
	# create markovchain object
	out <- new("markovchain", transitionMatrix = from, states = namesCandidate, byrow = byrow)
	
	invisible(out)
}


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
			rowDf <- data.frame(t0 = t0, t1 = t1, prob = prob)
			
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
			if((class(df[, i]) == "numeric") & (all(sapply(df[, i], .isProbRcpp) == TRUE))) {
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
  temp <- msm::pmatrix.msm(from)
  prMatr <- unclass(as.matrix(temp))
  out <- new("markovchain", transitionMatrix = prMatr)
  return(out)
}

# coerce msm object to markovchain object
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
  prMatr <- matlab::zeros(nelements)
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
setAs(from = "markovchain", to = "igraph", def = .mc2igraph)


# transposing method for markovchain objects
setMethod("t", "markovchain", 
		      function(x) { 
			      out <- new("markovchain", byrow = !x@byrow, 
			                 transitionMatrix = t(x@transitionMatrix))
			      
			      return(out)
		      } 
)


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

# compare two markovchain object
setMethod("==", c("markovchain", "markovchain"),
          function(e1, e2) {
            out <- FALSE
            out <- identical(e1@transitionMatrix, e2@transitionMatrix)
            return(out)
          }
)

setMethod("!=", c("markovchain", "markovchain"),
          function(e1, e2) {
            out <- FALSE
            out <- !(e1 == e2)
            return(out)
          }
)

# markovchain raise to some power
setMethod("^", c("markovchain", "numeric"),
          function(e1, e2) {
            out <- new("markovchain", states = e1@states, byrow = e1@byrow,
                       transitionMatrix = e1@transitionMatrix %^% e2,
                       name = paste(e1@name, "^", e2, sep = "")
                      )
            
            return(out)
          }
)


# methods to directly access transition matrix elements
setMethod("[", signature(x = "markovchain", i = "ANY", j = "ANY"),
          function(x, i, j) {
            out <- x@transitionMatrix[i, j]
            return(out)
          }
)

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
#' @export
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
		      definition = function(object, newdata, n.ahead = 1, continue = FALSE) {
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

# Wrapper for a function to lexicographically sort the rows of a matrixx
# m : matrix
.mcLexSort <- function(m) {
  matrix(unlist(.lexicographical_sort(m)), nrow=nrow(m), byrow = T)
}


#sort method for markovchain objects

setGeneric("sort", function(x, decreasing=FALSE) standardGeneric("sort"))

setMethod("sort", signature(x="markovchain"), function(x, decreasing=FALSE){
  
  #get matrix and state names 2 be sorted
 
  matr2besorted<-x@transitionMatrix 
  if (x@byrow==TRUE) states2besorted<-rownames(matr2besorted) else states2besorted<-colnames(matr2besorted)
  
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
