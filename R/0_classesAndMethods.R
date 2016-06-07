
# define Markov Chain class

setClass("markovchain", #class name
  representation(states = "character", byrow = "logical",
  transitionMatrix = "matrix", name = "character"),
  prototype(states = c("a","b"), byrow = TRUE, # prototypizing
  transitionMatrix=matrix(data = c(0,1,1,0),
  nrow=2, byrow=TRUE, dimnames=list(c("a","b"), c("a","b"))),
  name="Unnamed Markov chain")
)

# define Markov Chain class

setClass("markovchainList", 
		representation(markovchains = "list", 
		name = "character")
)

# verifies if a markovchainList object is valid function is valid

setValidity("markovchainList",
		function(object){
		check <- NULL
		for(i in length(object@markovchains))
		{
			if(class(object@markovchains[[i]]) != "markovchain") check <- "Error! All elements should be of class 'markovchain'" # All elemeents in the list should be a markovchain element
		}
		if(is.null(check)) check <- TRUE
		return(check)
	}
)

#initializing method for markovchain objects

setMethod("initialize",
		signature(.Object = "markovchain"),
		function (.Object, states, byrow, transitionMatrix,name,...) 
		{
			# put the standard markovchain 
			if(missing(transitionMatrix)) transitionMatrix<-matrix(data=c(0,1,1,0), #create a naive matrix
						nrow=2,
						byrow=TRUE, 
						dimnames=list(c("a","b"), c("a","b"))
			)
			
		# check names of transition matrix
			if(all(is.null(rownames(transitionMatrix)), is.null(colnames(transitionMatrix)))==TRUE) { #if all names are missing it initializes them to "1", "2",...
				if(missing(states)) {
					nr<-nrow(transitionMatrix)
					stateNames<-as.character(seq(1:nr))
				} else {stateNames<-states}
				
			rownames(transitionMatrix)<-stateNames
			colnames(transitionMatrix)<-stateNames
	 	} else if(is.null(rownames(transitionMatrix))) { #fix when rownames null
		  rownames(transitionMatrix)<-colnames(transitionMatrix)
		} else if(is.null(colnames(transitionMatrix))) { #fix when colnames null
			colnames(transitionMatrix)<-rownames(transitionMatrix)
		} else if(!setequal(rownames(transitionMatrix),colnames(transitionMatrix)))  colnames(transitionMatrix)=rownames(transitionMatrix) #fix when different
		if(missing(states)) states=rownames(transitionMatrix) #assign
		if(missing(byrow)) byrow=TRUE #set byrow as true by default
    if(missing(name)) name="Unnamed Markov chain"
		callNextMethod(.Object, states = states, byrow = byrow, transitionMatrix=transitionMatrix,name=name,...)
    }
)


# .isProb<-function(prob)
# {
# 	if (class(prob)!="numeric") return(FALSE)
# 	if (prob<0 | prob >1) return(FALSE)
# 	return(TRUE)
# }


# generic method to print out states

setGeneric("states", function(object) standardGeneric("states"))
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
#' @export
setGeneric("name<-", function(object, value) standardGeneric("name<-"))

#' @rdname setName
setMethod("name<-", "markovchain", 
          function(object, value) {
            object@name <- value
            object
          }
)

#adding a method names: to get names
setMethod("names","markovchain", 
          function(x) {
            out <- x@states
            return(out)
          }
)

#adding a method names: to set names
setMethod("names<-","markovchain", 
          function(x,value) {
            x@states<-value
            rownames(x@transitionMatrix)<-value
            colnames(x@transitionMatrix)<-value
            return(x)
          }
)


 # generic methods to get the dim of a markovchain and markovchainList

setMethod("dim","markovchain", 
		function(x) {
			out <- nrow(x@transitionMatrix)
			return(out)
		}
)

setMethod("dim","markovchainList", 
		function(x) {
			  out <- length(x@markovchains)
			  return(out)
		  }
)


 # method  to set the validity of a markovchain object

setValidity("markovchain",
		function(object) {
			check<-NULL
			# performs a set of check whose results are saved in check
			if (any(sapply(as.numeric(object@transitionMatrix),.isProbRcpp))==FALSE) check <- "Error! Some elements are not probabilities" #checks if probability
			if (object@byrow==TRUE) {
				if(any(round(rowSums(object@transitionMatrix),5)!=1)) check <- "Error! Row sums not equal to one"
			} else {
				if(any(round(colSums(object@transitionMatrix),5)!=1)) check <- "Error! Col sums not equal to one"
			} #checks if col sums not equal to one
			if (nrow(object@transitionMatrix)!=ncol(object@transitionMatrix)) check <- "Error! Not squared matrix" #check if squalre matrix
            if (!setequal(colnames(object@transitionMatrix),object@states)) check <- "Error! Colnames <> states" #checks if 
            if (!setequal(rownames(object@transitionMatrix),object@states)) check <- "Error! Rownames <> states"
			if ( is.null(check) ) return(TRUE) else return(check)
		}
)

.mcEigen<-function(matr, transpose=TRUE)
{
  # Function to extract eigenvalues, core of get steady states 
  #
  # Args:
  # matr: the matrix to extract
  # transpose:  boolean indicating whether the matrx shall be transpose
  #
  # Results:
  # a matrix / vector
  if (transpose) tMatr <- t(matr) else tMatr <- matr #trasposing
  eigenResults <- eigen(x=tMatr,symmetric=FALSE) #perform the eigenvalue extraction
  onesIndex <- which(round(eigenResults$values,3)==1) #takes the one eigenvalue
  #do the following: 1:get eigenvectors whose eigenvalues==1
  #2: normalize
  if (length(onesIndex)==0) {
    warning("No eigenvalue = 1 found")
    return(NULL)
  }
  if (transpose==TRUE)
  {
    eigenTake <- as.matrix(t(eigenResults$vectors[,onesIndex])) 
    out <- eigenTake/rowSums(eigenTake) 
  } else {
    eigenTake <- as.matrix(eigenResults$vectors[,onesIndex]) 
    out <- eigenTake/colSums(eigenTake)
  }
  # subset the eigenvectors
  # normalize
  # take the real part: need to be sanitized
	#@TAE: later we have to see and optimize this part. I am not sure taking
	#the real part is the most appropriate.
  out <- Re(out)
  return(out)
}

#method to get stationary states
setGeneric("steadyStates", function(object) standardGeneric("steadyStates"))
setMethod("steadyStates","markovchain", 
		function(object) {
			transposeYN <- FALSE
			if(object@byrow==TRUE) transposeYN <- TRUE		
            out<-.mcEigen(matr=object@transitionMatrix, transpose=transposeYN) #wrapper for .mcEigen
            if(is.null(out)) {
				warning("Warning! No steady state")
				return(NULL)
            }
			if(transposeYN==TRUE) { 
				colnames(out) <- object@states
			} else {
				rownames(out) <- object@states
			}
            #if(nrow(out)==1) out<-as.numeric(out)
            return(out)
          }
)


# generic function to extract absorbing states

#' @rdname absorbingStates
#' 
#' @export
setGeneric("absorbingStates", function(object) standardGeneric("absorbingStates"))
setMethod("absorbingStates","markovchain", 
          function(object) {
			  out <- character()
			  matr <- object@transitionMatrix #extract the byrow transition matrix
			  transposeYN <- FALSE
			  if(object@byrow==TRUE) transposeYN <- TRUE
			  steadyStates <- .mcEigen(matr=matr, transpose=transposeYN) #checkk
			  if(is.null(steadyStates)) return(character(0))
			  #identify which states are absorbing if they are diagonally one
			  if(transposeYN==TRUE) maxCols <- apply(steadyStates, 2, "max") else maxCols <- apply(steadyStates, 1, "max")  
			  index <- which(maxCols==1)
			  if(length(index)>0) out <- object@states[index] 
			  return(out)
          }
)

# generic method to extract transient states

#' @rdname absorbingStates
#' 
#' @export
setGeneric("transientStates", function(object) standardGeneric("transientStates"))
setMethod("transientStates","markovchain", 
		function(object) {
			out <- character()
			matr <- object@transitionMatrix #extract the byrow transition matrix
			temp <- .commclassesKernelRcpp(matr)
			index <- which(temp$v==FALSE)
			if(length(index)>0) out <- names(temp$v[index])
			return(out)
		}
)

#generic method to extract transition probability
setGeneric("transitionProbability", function(object, t0, t1) standardGeneric("transitionProbability"))
setMethod("transitionProbability","markovchain", 
	function(object, t0, t1) {
		out <- numeric(1)
		fromState <- which(object@states==t0)
		toState <- which(object@states==t1)
		out<-ifelse(object@byrow==TRUE, object@transitionMatrix[fromState, toState] , object@transitionMatrix[toState, fromState])
		return(out)
	}
)

#print plot show methods

.showInt <- function(object, verbose=TRUE)
{
	if (object@byrow==TRUE) direction="(by rows)" else direction="(by cols)"
	if (verbose==TRUE) cat(object@name,"\n A ",dim(object),"- dimensional discrete Markov Chain characterized by following states: \n",
	                       paste(states(object),collapse=", "), "\n The transition matrix  ", 
	                       direction," is defined as follows: \n")
	print(object@transitionMatrix)
	cat("\n")
}


#show methods for markovchain and markovchain list objects 
setMethod("show","markovchain", #markovchain
          function(object){
          .showInt(object)
          }
)

setMethod("show", "markovchainList",
          function(object){
		  cat(object@name, " list of Markov chain(s)","\n")
          for(i in 1:length(object@markovchains)) 
          {
            cat("Markovchain ",i,"\n")
            show(object@markovchains[[i]])
          }
          }
)

#print methods

setMethod("print","markovchainList",function(x) show(x))
setMethod("print","markovchain", #metodo print
          function(x){
           object <- x
		   .showInt(object,verbose=FALSE)
          }
)

.getNet<-function(object,round=FALSE)
{
 	# function to get the absorbency matrix to plot and export to igraph
	#
 	# Args: 
	# object: a markovchain object
	# round: boolean to round
 	#
	# Returns:
	#
	# a graph adjacenty
	if (object@byrow==FALSE) object <- t(object)
	#with this fails 
	#matr <- Matrix(data=object@transitionMatrix, sparse = TRUE)*100 #need to multiply by 100
	#with this works. 
	matr<-object@transitionMatrix*100
	if(round==TRUE) matr <- round(matr,2)
	net <- graph.adjacency(adjmatrix=matr, weighted=TRUE, mode="directed")
	#net<-graph_from_adjacency_matrix(adjmatrix = matr,weighted = TRUE,mode="directed")
	return(net)
}


 # Plot methods for markovchain objects



#plot method from stat5
setMethod("plot", signature(x="markovchain", y="missing"),
		function(x, y, package="igraph",...){
		  switch(package,
		         diagram = {
		           if (requireNamespace("diagram", quietly = TRUE)) {
		             .plotdiagram(object=x,...)
		           } else {
		            # cat("diagram unavailable, using standard method","\n")
		             netMc <- .getNet(object=x,round=TRUE)
		             edgeLabel <- round(E(netMc)$weight/100,2)
		             plot.igraph(x=netMc,edge.label=edgeLabel, ...)
		           }
		           
		         },
		         DiagrammeR= {
		           if (requireNamespace("DiagrammeR", quietly = TRUE)) {
		             .plotDiagrammeR(object=x,...)
		           } else {
		         #    cat("DiagrammeR unavailable, using standard method","\n")
		             netMc <- .getNet(object=x,round=TRUE)
		             edgeLabel <- round(E(netMc)$weight/100,2)
		             plot.igraph(x=netMc,edge.label=edgeLabel, ...)
		           }
		          
		         },
		         {
		           netMc <- .getNet(object=x,round=TRUE)
		           edgeLabel <- round(E(netMc)$weight/100,2)
		           plot.igraph(x=netMc,edge.label=edgeLabel, ...)
		         })
		}
)


#@TAE: create an internal function that does this. Check also if the canonic form function 
#is appropriate

 # method to convert into canonic form a markovchain object
 # TODO: check meaninsg of this function

#' @rdname absorbingStates
#' 
#' @export
setGeneric("canonicForm",function(object) standardGeneric("canonicForm"))
setMethod("canonicForm","markovchain",
          function(object)
          {
            P <- object@transitionMatrix
            comclasList <- .commclassesKernelRcpp(P)
            vu <- comclasList$v
			
            u <- find(vu==TRUE)
            w <- find(vu==FALSE)
            
            Cmatr <- comclasList$C
            R <- numeric()
            while(length(u)>0)
            {
              R <- c(R,u[1])
              vu <- as.logical(vu*(Cmatr[u[1],]==FALSE));
              u <- find(vu==TRUE);
            }
            p <- numeric()
            for (i in 1:length(R))
            {
              a <- find(Cmatr[R[i],])
              p <- c(p,a)
            }
            p <- c(p, w)
            Q <- P[p,p]
            out<-new("markovchain",transitionMatrix=Q,name=object@name)
            return(out)
          }
)

.canonicForm<-function(object)
{
  P <- object@transitionMatrix
  comclasList <- .commclassesKernelRcpp(P)
  vu <- comclasList$v
  
  u <- find(vu==TRUE)
  w <- find(vu==FALSE)
  
  Cmatr <- comclasList$C
  R <- numeric()
  while(length(u)>0)
  {
    R <- c(R,u[1])
    vu <- as.logical(vu*(Cmatr[u[1],]==FALSE));
    u <- find(vu==TRUE);
  }
  p <- numeric()
  for (i in 1:length(R))
  {
    a <- find(Cmatr[R[i],])
    p <- c(p,a)
  }
  p <- c(p, w)
  Q <- P[p,p]
  out<-new("markovchain",transitionMatrix=Q,name=object@name)
  return(out)
}

# summary method for markovchain class
# lists: closed, transient classes, irreducibility, absorbint, transient states
setMethod("summary", signature(object="markovchain"),
		function(object){
			outs <- .summaryKernelRcpp(object)
			cat(object@name," Markov chain that is composed by:","\n")
			check <- length(outs$closedClasses)
			cat("Closed classes:","\n")
			if(check==0) cat("NONE","\n") else {
				for(i in 1:check) cat(outs$closedClasses[[i]],"\n")
			}
			check <- length(outs$recurrentClasses)
			cat("Recurrent classes:","\n")
			if(check==0) cat("NONE","\n") else {
			  cat("{")
			  cat(outs$recurrentClasses[[1]], sep=",")
			  cat("}")
			  if(check > 1) {
			    for(i in 2:check) {
			      cat(",{")
			      cat(outs$recurrentClasses[[i]],sep=",")
			      cat("}")
			    }
			  }
			  cat("\n")
			}
			check <- length(outs$transientClasses)
			cat("Transient classes:","\n")
			if(check==0) cat("NONE","\n") else {
				# for(i in 1:check) cat(outs$transientClasses[[i]],"\n")
			  cat("{")
			  cat(outs$transientClasses[[1]], sep=",")
			  cat("}")
			  if(check > 1) { 
			    for(i in 2:check) {
			      cat(",{")
			      cat(outs$transientClasses[[i]],sep=",")
			      cat("}")
			    }
			  }
			  cat("\n")
			}
			irreducibility <- is.irreducible(object)
			if(irreducibility) cat("The Markov chain is irreducible","\n") else cat("The Markov chain is not irreducible","\n")
			check <- absorbingStates(object)
			if(length(check)==0) check <- "NONE"
			cat("The absorbing states are:",check )
			cat("\n")
			invisible(outs) 
          }
)
##################################################AS METHODS#########################

.checkMatrix <- function(matr, byrow=TRUE, verbose=FALSE) {
	#first check: size
	if (dim(matr)[1]!=dim(matr)[2]) {
		if(verbose) stop("Error! Rectangular matrix")
		return(FALSE)
	}
	
	#second check: all elements are probs
	
	for(i in 1:nrow(matr)) {
		for(j in 1:ncol(matr)){
			if(!(.isProbRcpp(matr[i,j]))){
				if(verbose) stop("Error! Non probabilities")
				return(FALSE)}
		}
	}
	
	#third check: either columns or rows sum to one
		
	if(byrow==FALSE) matr <- t(matr) #to perform only one check
	
	check <- rowSums(matr)
	
	for( i in 1:length(check)) if (abs(1-check[i]) > .Machine$double.eps) {
			if(verbose) stop("Error! Either rows or cols should sum to 1")
			return(FALSE)
		}
	#if all test are passed
	return(TRUE)
}

.matrix2Mc<-function(from){
	
	#Internal function to return a markovchain object given a matrix
	
	byrow<-FALSE
	checkByRow <- .checkMatrix(from,byrow = TRUE)
	
	if(checkByRow) byrow=TRUE else  {
		checkByCols <- .checkMatrix(from,byrow = FALSE)
		if(!checkByCols) stop("Error! Not a probability matrix")	
	}
	
	if(byrow==TRUE) namesCandidate <- rownames(from) else namesCandidate <- colnames(from)
	if(is.null(namesCandidate)) {
		namesCandidate <- character()
		for(i in 1:nrow(from)) namesCandidate<-c(namesCandidate, paste("s",i,sep=""))
	}
	
	out<-new("markovchain",transitionMatrix=from, states=namesCandidate,byrow=byrow)
	
	invisible(out)

}

setAs(from="matrix", to="markovchain", def=.matrix2Mc)




.mc2Df<-function(from)
{
# Function to transform a markovchain into a data.frame
# Args:
#
# from: a markovchain object
#
# Returns:
# a data.frame
	nr<-nrow(from@transitionMatrix) #cycles for all rows and columns
	for(i in 1:nr){
		for(j in 1:nr){
			t0 <- from@states[i]
			t1 <- from@states[j]
			prob <- transitionProbability(object=from, t0=t0, t1=t1)
			rowDf <- data.frame(t0=t0, t1=t1, prob=prob)
			if(exists("outDf")) outDf <- rbind(outDf, rowDf) else outDf <- rowDf
		}
	}
	return(outDf)
}

#method to convert from markovchain to data.frame

setAs(from="markovchain", to="data.frame", def=.mc2Df)

.whichColProb<-function(df)
{
	out=0
	if(ncol(df)>3) warning("Warning! More than three column. Only the first three will be used")
	if(ncol(df)<3) stop("Error! Three column needed")
	
	for(i in 1:ncol(df))
		{
			if((class(df[,i])=="numeric")&(all(sapply(df[,i], .isProbRcpp)==TRUE))) #when found the first numeric and probability col
				{
					out=i
					break
				}
		}
	return(out)
}

.df2Mc<-function(from)
 {
	 #
	 # Function to convert from a data.frame containing initial, ending and probability columns to a proper markovchain object
	 #
	 # Args:
	 #
	 # from: a data.frame
	 #
	 # Returns:
	 #
	 # A markovchain object 
	statesNames <- unique(from[,1])
	colProb <- .whichColProb(from)
	prMatr <- zeros(length(statesNames))
	rownames(prMatr)<-statesNames
	colnames(prMatr)<-statesNames
	for(i in 1:nrow(from)) {
		idRow <- which(statesNames==from[i,1]) #assume first col from
		idCol <- which(statesNames==from[i,2]) #assume second row to
		prMatr[idRow,idCol]<-from[i,3]
	}
 	out<-new("markovchain", transitionMatrix=prMatr)
	return(out)
 }

#method to convert 
setAs(from="data.frame", to="markovchain", def=.df2Mc)


.table2Mc<-function(from)
{
	#checks
	if(dim(from)[1]!=dim(from)[2]) stop("Error! Table is not squared")
	if(!setequal(rownames(from),colnames(from))) stop("Error! Rows not equal to coulumns")
	#temp <- as.matrix(from) #compatibility issues raised by Kurt
  temp <- unclass(as.matrix(from)) #following Kurt's suggestion
	fromMatr <- temp[,order(rownames(temp))] #makes the same sequence of col / row
	outMatr <- fromMatr/rowSums(fromMatr)
	out <- new("markovchain",states=rownames(temp), transitionMatrix=outMatr, byrow=TRUE)
	return(out)
}

setAs(from="table", to="markovchain", def=.table2Mc)

#function from msm to markovchain

.msm2Mc<-function(from)
{
  temp <- msm::pmatrix.msm(from)
  prMatr <- unclass(as.matrix(temp))
  out<-new("markovchain", transitionMatrix=prMatr)
  return(out)
}

setAs(from="msm", to="markovchain", def=.msm2Mc)

#function for msm.est to mc. Assume a probability matrix given

.msmest2Mc<-function(from)
{
  
  if (is.matrix(from))
    pMatr <- from #central estimate
  if (is.list(from))
    pMatr <- from[[1]] #central estimate
  
  out<-new("markovchain", transitionMatrix=as(pMatr,"matrix")) #need force matrix
  return(out)
}

setAs(from="msm.est", to="markovchain", def=.msmest2Mc)


#function from etm to markovchain

.etm2Mc<-function(from)
{
  df<-from$trans
  elements<-from$state.names
  nelements<-length(elements)
  prMatr<-zeros(nelements)
  dimnames(prMatr) = list(elements, elements)
  for(i in 1:dim(df)[1]) {
    r<-df[i,]
    stateFrom<-r$from
    stateTo<-r$to
    prMatr[stateFrom, stateTo]<-prMatr[stateFrom, stateTo] + 1
  }
  rsums<-rowSums(prMatr)
  prMatr<-prMatr/rsums
  if(any(rsums == 0)) {
    indicesToBeSanitized<-which(rsums==0)
    for(i in indicesToBeSanitized) {
      for(j in 1:nelements) prMatr[i,j]<-1/nelements
    }
  }
  out<-new("markovchain", transitionMatrix=prMatr)
  return(out)
}

setAs(from="etm", to="markovchain", def=.etm2Mc)

#functions and methods to return a matrix

.mc2matrix<-function(from)
{
	out <- from@transitionMatrix
	return(out)
}

setAs(from="markovchain", to="matrix", def=.mc2matrix)

#functions and methods to return a matrix

.mc2igraph<-function(from)
{
	temp <- .mc2Df(from) #convert the markovchain to data.frame
	out <- graph.data.frame(temp) #convert the data frame to igraph graph
	return(out) #return
}

setAs(from="markovchain", to="igraph", def=.mc2igraph)


#transposing method for markovchain objects

setMethod("t", "markovchain", 
		function(x) {
			out <- new("markovchain", byrow=!x@byrow, transitionMatrix=t(x@transitionMatrix))
			return(out)
		} 
)

#multiplicationMethod
#by a markov chain (like a power by 2)
setMethod("*", c("markovchain", "markovchain"),
function(e1, e2) {
	
	#
	# function to multiplicate two markov chains
	#
	# Args:
	#
	# e1: first markovchain
	# e2: second markov chain
	#
	# Returns:
	# if feasible, a markovchain where the transition matrix i e1*e2
			
			if(!setequal(e1@states,e1@states)) warning("Warning! Different states")
			if(!setequal(dim(e1@transitionMatrix), dim(e2@transitionMatrix))) stop("Error! Different size")
			if(!(e1@byrow==e2@byrow)) stop("Error! Both transition matrix should be defined either by row or by column")
			newStates <- e1@states
			newTransMatr <- e1@transitionMatrix%*%e2@transitionMatrix
			byRow <- e1@byrow
			mcName <- e1@name #multiplicated matrix takes the first matrix's name
			out<-new("markovchain", states=newStates, transitionMatrix=newTransMatr, byrow=byRow, name=mcName)
			return(out)
}
)

#methods implemented for multplicating with matrices, vectors, etc...

setMethod("*", c("matrix", "markovchain"),
		function(e1, e2) {
			out <- e1%*%e2@transitionMatrix
			return(out)
		}
)

setMethod("*", c("markovchain","matrix"),
		function(e1, e2) {
			out <- e1@transitionMatrix%*%e2
			return(out)
}
)

setMethod("*", c("numeric","markovchain"),
		function(e1, e2) {
			if(length(e1)!=dim(e2)) stop("Error! Uncompatible dimensions") else out <- e1%*%e2@transitionMatrix
			return(out)
		}
)

setMethod("*", c("markovchain","numeric"),
		function(e1, e2) {
			if(length(e2)!=dim(e1)) stop("Error! Uncompatible dimensions") else out <- e1@transitionMatrix%*%e2
			return(out)
		}
)


setMethod("==", c("markovchain","markovchain"),
          function(e1, e2) {
            out <- FALSE
            out <-
              identical(e1@transitionMatrix, e2@transitionMatrix)
            return(out)
          })

setMethod("!=", c("markovchain","markovchain"),
          function(e1, e2) {
            out <- FALSE
            out <-
              (!(identical(
                e1@transitionMatrix, e2@transitionMatrix
              )))
            return(out)
          })

setMethod("^", c("markovchain", "numeric"),
          function(e1, e2) {
            out <-
              new(
                "markovchain", states = e1@states, byrow = e1@byrow,transitionMatrix = e1@transitionMatrix %^%
                  e2,
                name = paste(e1@name,"^",e2,sep = "")
              )
            return(out)
          })


#methods to directly access transition matrix elements

setMethod("[",
          signature(x = "markovchain", i = "ANY", j = "ANY"),
          function(x, i, j) {
            out <- x@transitionMatrix[i,j]
            return(out)
          })

#methods to directly access markovchain objects composing a markovchainList object


setMethod("[[",
		signature(x = "markovchainList", i = "ANY"),
		function(x, i) {
			out <- x@markovchains[[i]]
            return(out)
		}
)

setGeneric("conditionalDistribution", function(object,state) standardGeneric("conditionalDistribution"))
setMethod("conditionalDistribution","markovchain", #metodo plot
          function(object,state){
			stateNames <- states(object) #get the names
			out <- numeric(length(stateNames)) #allocater oiutvect
			index2Take <- which(stateNames==state) #states are assumed to be sorted
			if(object@byrow==TRUE) #returns the probability vector depending by sorting
			{
				out <- object@transitionMatrix[index2Take,]
			} else {
				out <- object@transitionMatrix[,index2Take]
			}	
			#names the output and returs it
			names(out)<-stateNames
			return(out) 
		}
)
		  

#geth the mode of a probability vector
.getMode<-function(probVector,ties="random")
{
	#
	# Function to get the mode of a probability vector
	# 
	# Args:
	#
	# probvector: the probability vector
	# ties: specifies if ties are to be sampled, otherwise more than one element is returned
	#
	# Returns:
	#
	# the name of the model element
	maxIndex <- which(probVector==max(probVector))
	temp <- probVector[maxIndex]
	if((ties=="random")&(length(temp)>1)) out <- sample(temp,1) else out <- temp
	return(names(out))
}

#predict method for markovchain objects

setMethod("predict","markovchain", 
		function(object,newdata,n.ahead=1) {
			lastState <- newdata[length(newdata)] #identify the last state
            out <- character()
            for(i in 1:n.ahead)
            {
              #cyclically determine the most probabile subsequent state from the conditional distribution
              newState <- .getMode(probVector=conditionalDistribution(object,lastState), ties="random") 
              out <- c(out,newState)
              lastState <- newState
            }
            return(out)
          }
)

setMethod("predict","markovchainList",
		definition = function(object,newdata,n.ahead=1,continue=FALSE) {
			#object a markovchainList, newdata=the actual data, n.ahead=how much ahead, continue=veryfi if thake last
			out<-character() #allocate output
			actualPos<-length(newdata) 
			lastState<-newdata[actualPos] #takes last position
			for(i in 1:n.ahead) #cycles
			{
				newPos <- actualPos + i - 1 #increments
				if(newPos<=dim(object)) #applies if within the length of markovchainList
				{
					newState <- predict(object=object[[newPos]], newdata = lastState, n.ahead=1)
					out <- c(out,newState)
					lastState <- newState
				} else {
					if(continue==TRUE) #applies if allowed to predict over last state
					{
						newState <- predict(object=object[[dim(object)]], newdata = lastState, n.ahead=1)
						out <- c(out,newState)
						lastState <- newState
					} else break;
				} #chiude else
			} #chiude for
			return(out)
		}

)