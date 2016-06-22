#' Function to generate a sequence of states from homogeneous Markov chains.
#' 
#' Provided any \code{markovchain} object, it returns a sequence of 
#' states coming from the underlying stationary distribution. 
#' 
#' @param n Sample size
#' @param markovchain \code{markovchain} object
#' @param t0 The initial state
#' @param include.t0 Specify if the initial state shall be used
#' @param useRCpp Boolean. Should RCpp fast implementation being used? Default is yes.
#' 
#' @details A sequence of size n is sampled.
#' 
#' @return A Character Vector
#' 
#' @references A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
#' 
#' @author Giorgio Spedicato
#' 
#' @seealso \code{\link{markovchainFit}}
#' 
#' @examples 
#' # define the markovchain object
#' statesNames <- c("a", "b", "c")
#' mcB <- new("markovchain", states = statesNames, 
#'    transitionMatrix = matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), 
#'    nrow = 3, byrow = TRUE, dimnames = list(statesNames, statesNames)))
#' 
#' # show the sequence
#' outs <- markovchainSequence(n = 100, markovchain = mcB, t0 = "a")
#' 
#' @export

markovchainSequence <-function (n, markovchain, t0 = sample(markovchain@states, 1),
                               include.t0 = FALSE, useRCpp = TRUE) {
  
  # check whether given initial state is possible state or not
  if (!(t0 %in% markovchain@states))
    stop("Error! Initial state not defined")
  
  # call to cpp implmentation of markovchainSequence
  if (useRCpp) {
    return(.markovchainSequenceRcpp(n, markovchain, t0, include.t0))
  }
  
  # R implementation of the function
  
  # create a sequence of size n initially not initialized
  chain <- rep(NA,n)
  
  # initial state
  state <- t0
  
  # populate the sequence
  for (i in 1:n) {
    # row probabilty corresponding to the current state
    rowProbs <- markovchain@transitionMatrix[which(markovchain@states == state), ]
    
    # select the next state
    outstate <- sample(size = 1, x = markovchain@states, prob = rowProbs)
    
    # store the new state
    chain[i] <- outstate
    
    # update the current state
    state <- outstate
  }
  
  # output
  out <- chain
  
  # whether to include initial state or not
  if (include.t0) {
    out <- c(t0, out)
  }
  
  return(out)
}


##################
# random sampler #
##################

# check if the subsequent states are included in the previous ones

# TODO: too strong contidion; should be changed by checking that
# all states that can be reached in one step at t-1 are named  
# in object[[t]]

# check the validity of non homogeneous markovchain list
# object is a list of markovchain object
.checkSequence <- function(object) {
  # assume non homogeneous markovchain list is valid
  out <- TRUE
  
  # list of one transition matrix implies valid
  if (length(object) == 1) {
    return(out) 
  }
  
  # if number of transition matrices are more than one  
  for (i in 2:length(object)) {
    
    # select the states which are reachable in one step
    if(object[[i - 1]]@byrow) {
      reachable <- (colSums(object[[i - 1]]@transitionMatrix) != 0)
    } else {
      reachable <- (rowSums(object[[i - 1]]@transitionMatrix) != 0)
    }
    
    # possible states in the previous markovchain object
    statesNm1 <- states(object[[i - 1]])[reachable]
    
    # possible states in the current markovchain object
    statesN <- states(object[[i]])
    
    # common states 
    intersection <- intersect(statesNm1, statesN)
    
    # condition to check whether statesNm1 is a subset of statesN or not
    if (setequal(intersection, statesNm1) == FALSE) {
      out <- FALSE
      break
    }
    
  }
  
  return(out)
}

#' Function to generate a sequence of states from homogeneous or non-homogeneous Markov chains.
#' 
#' Provided any \code{markovchain} or \code{markovchainList} objects, it returns a sequence of 
#' states coming from the underlying stationary distribution. 
#' 
#' @param n Sample size
#' @param object Either a \code{markovchain} or a \code{markovchainList} object
#' @param what It specifies whether either a \code{data.frame} or a \code{matrix} 
#'        (each rows represent a simulation) or a \code{list} is returned.
#' @param useRCpp Boolean. Should RCpp fast implementation being used? Default is yes.
#' @param parallel Boolean. Should parallel implementation being used? Default is yes.
#' @param num.cores Number of Cores to be used
#' @param ... additional parameters passed to the internal sampler
#' 
#' @details When a homogeneous process is assumed (\code{markovchain} object) a sequence is 
#' sampled of size n. When an non - homogeneous process is assumed,
#' n samples are taken but the process is assumed to last from the begin to the end of the 
#' non-homogeneous markov process.
#' 
#' @return Character Vector, data.frame, list or matrix
#' 
#' @references A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
#' 
#' @author Giorgio Spedicato
#' 
#' @note Check the type of input
#' 
#' @seealso \code{\link{markovchainFit}}
#' 
#' @examples 
#' # define the markovchain object
#' statesNames <- c("a", "b", "c")
#' mcB <- new("markovchain", states = statesNames, 
#'    transitionMatrix = matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), 
#'    nrow = 3, byrow = TRUE, dimnames = list(statesNames, statesNames)))
#' 
#' # show the sequence
#' outs <- rmarkovchain(n = 100, object = mcB, what = "list")
#' 
#' 
#' #define markovchainList object
#' statesName = c("a", "b", "c")
#' mcA <- new("markovchain", states = statesNames, transitionMatrix = 
#'    matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), nrow = 3, 
#'    byrow = TRUE, dimnames = list(statesNames, statesNames)))
#' mcB <- new("markovchain", states = statesNames, transitionMatrix = 
#'    matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), nrow = 3, 
#'    byrow = TRUE, dimnames = list(statesNames, statesNames)))
#' mcC <- new("markovchain", states = statesNames, transitionMatrix = 
#'    matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), nrow = 3, 
#'    byrow = TRUE, dimnames = list(statesNames, statesNames)))
#' mclist <- new("markovchainList", markovchains = list(mcA, mcB, mcC)) 
#' 
#' # show the list of sequence
#' rmarkovchain(100, mclist, "list")
#'      
#' @export

rmarkovchain <- function(n, object, what = "data.frame", useRCpp = TRUE, parallel = FALSE, num.cores = NULL, ...) {
  
  # check the class of the object
  if (class(object) == "markovchain") {
    out <- markovchainSequence(n = n, markovchain = object, useRCpp = useRCpp, ...)
    return(out)
  }
    
  if (class(object) == "markovchainList")
  {
    #######################################################
    if(useRCpp && !parallel) {
      
      # if include.t0 is not passed as extra argument then set include.t0 as false
      include.t0 <- list(...)$include.t0
      include.t0 <- ifelse(is.null(include.t0), FALSE, include.t0)
      
      # check whether initial state is passed or not
      t0 <- list(...)$t0
      if (is.null(t0)) t0 <- character()
      
      # call fast cpp function
      dataList <- .markovchainListRcpp(n, object@markovchains, include.t0, t0)
      
      # format in which results to be returned
      if (what == "data.frame") {
        out <- data.frame(iteration = dataList[[1]], values = dataList[[2]])
      }
        
      else {
        # output in matrix format
        # each row is an independent sequence
        out <- matrix(data = dataList[[2]], nrow = n, byrow = TRUE)
        
        # output in list format
        if (what == "list") {
          outlist <- list()
          for (i in 1:nrow(out))
            outlist[[i]] <- out[i, ]
          out <- outlist
        }
      } 
      return(out)
    }
    ##########################################################
    if(useRCpp && parallel) {
      
      # Calculate the number of cores
      # It's not good to use all cores
      no_cores <- max(1,parallel::detectCores() - 1)
      
      # number of cores specified should be less than or equal to maximum cores available
      if((! is.null(num.cores))  && num.cores <= no_cores + 1 && num.cores >= 1) {
        no_cores <- num.cores
      }
      
      RcppParallel::setThreadOptions(no_cores)
      
      # if include.t0 is not passed as extra argument then set include.t0 as false
      include.t0 <- list(...)$include.t0
      include.t0 <- ifelse(is.null(include.t0), FALSE, include.t0)
      
      # check whether initial state is passed or not
      t0 <- list(...)$t0
      if (is.null(t0)) t0 <- character()
      
      dataList <- .markovchainSequenceParallelRcpp(object, n, include.t0, t0)
      
      if(what == "list") return(dataList)
      
      # dimension of matrix to be returned
      nrow <- length(dataList)
      ncol <- length(dataList[[1]])   
      
      if(what == "matrix") {
        out <- matrix(nrow = nrow, ncol = ncol)
        for(i in 1:nrow) out[i, ] <- dataList[[i]]
        return(out)
      }
      
      iteration <- numeric()
      values <- character()
      
      # if what id data frame
      for(i in 1:nrow) {
        iteration <- append(iteration, rep(i, ncol))
        values <- append(values, dataList[[i]])
      }
      
      return(data.frame(iteration = iteration, values = values))
    }
    
    ##########################################################
    if(!useRCpp && parallel) {
      # if include.t0 is not passed as extra argument then set include.t0 as false
      include.t0 <- list(...)$include.t0
      include.t0 <- ifelse(is.null(include.t0), FALSE, include.t0)
      
      # check whether initial state is passed or not
      t0 <- list(...)$t0
      if (is.null(t0)) t0 <- character()
      
      dataList <- .markovchainSequenceParallel(n, object, t0, num.cores, include.t0)
      
      if(what == "list") return(dataList)
      
      # dimension of matrix to be returned
      nrow <- length(dataList)
      ncol <- length(dataList[[1]])   
      
      if(what == "matrix") {
        out <- matrix(nrow = nrow, ncol = ncol)
        for(i in 1:nrow) out[i, ] <- dataList[[i]]
        return(out)
      }
      
      iteration <- numeric()
      values <- character()
      
      # if what id data frame
      for(i in 1:nrow) {
        iteration <- append(iteration, rep(i, ncol))
        values <- append(values, dataList[[i]])
      }
      
      return(data.frame(iteration = iteration, values = values))
      
    }
    ##########################################################
    
    # store list of markovchain object in object
    object <- object@markovchains
    
    # check the validity of markovchainList object
    verify <- .checkSequence(object = object)
    
    # show warning if sequence is invalid
    if (!verify) {
      warning("Warning: some states in the markovchain sequences are not contained in the following states!")
    }
    
    # helper vector
    iteration <- numeric()
    values <- character()
    
    # create one sequence in each iteration
    for (i in 1:n) {
      
      # the first iteration may include initial state
      sampledValues <- markovchainSequence(n = 1, markovchain = object[[1]], ...)
      outIter <- rep(i, length(sampledValues))
      
      # number of markovchain objects are more than one
      if (length(object) > 1) {
        for (j in 2:length(object)) {
          pos2take <- length(sampledValues)
          # select new state of the sequence from the old state
          # t0 refers to the old state
          newVals <-markovchainSequence(n = 1, markovchain = object[[j]], t0 = sampledValues[pos2take]) 
          
          # update in every iteration
          outIter <- c(outIter, i)
          sampledValues <- c(sampledValues, newVals)
        }
      }
      
      # populate the helper vectors
      iteration <- c(iteration, outIter)
      values <- c(values, sampledValues)
    }
    
    # defining the output
    if (what == "data.frame") {
      out <- data.frame(iteration = iteration, values = values)
    } else {
      # ouput in matrix format
      out <- matrix(data = values, nrow = n, byrow = TRUE)
      
      # store each row of the matrix in the list
      if (what == 'list') {
        outlist <- list()
        for (i in 1:nrow(out))
          outlist[[i]] <- out[i, ]
        out <- outlist
      }
    }
  }
  
  return(out)
}

######################################################################

# helper function to calculate one sequence
.markovchainSPHelper <- function(x, t0, mclist, include.t0) {
  # number of transition matrices
  n <- length(mclist@markovchains)
  
  # take care of initial state
  vin <- 0
  if(include.t0) vin <- 1
  
  # a character vector to store a single sequence
  seq <- character(length = n + vin)
  
  if(length(t0) == 0) {
    stateNames <- mclist@markovchains[[1]]@states 
    t0 <- sample(x = stateNames, size = 1,  prob = rep(1 / length(stateNames), length(stateNames)))
  }
  
  if(include.t0) seq[1] <- t0
  
  # calculate one element of sequence in each iteration
  for (i in 1:n) {
    stateNames <- mclist@markovchains[[i]]@states 
    byRow <- mclist@markovchains[[i]]@byrow
    
    # check whether transition matrix follows row-wise or column-wise fashion
    if(byRow) prob <- mclist@markovchains[[i]]@transitionMatrix[which(stateNames == t0), ]
    else prob <- mclist@markovchains[[i]]@transitionMatrix[, which(stateNames == t0)]
    
    # initial state for the next transition matrix
    t0 <- sample(x = stateNames, size = 1,  prob = prob)
    
    # populate the sequence vector
    seq[i+vin] <- t0
  }
  
  return(seq)
}

# Function to generate a list of sequence of states in parallel from non-homogeneous Markov chains.
# 
# Provided any markovchainList object, it returns a list of sequence of states coming 
# from the underlying stationary distribution. 
# 
# @param  n Sample size
# @param object markovchainList object
# @param t0 Initial state
# @param num.cores Number of cores
#   

.markovchainSequenceParallel <- function(n, object,
                                        t0 = character(),
                                        num.cores = NULL, include.t0 = FALSE) {
  # check for the validity of non-uniform markov chain
  verify <- .checkSequence(object@markovchains)
  if (!verify) {
    warning("Warning: some states in the markovchain sequences are not contained in the following states!")
  }
    
  # Calculate the number of cores
  # It's not good to use all cores
  no_cores <- max(1,parallel::detectCores() - 1)
  
  # number of cores specified should be less than or equal to maximum cores available
  if((! is.null(num.cores))  && num.cores <= no_cores + 1 && num.cores >= 1) {
    no_cores <- num.cores
  }
  
  # Initiate cluster
  cl <- parallel::makeCluster(no_cores)
  
  # export the variables to be used in the helper function
  # parallel::clusterExport(cl, "t0")
  
  # export the variables to be used in the helper function
  mclist <- object
  # parallel::clusterExport(cl, "mclist")
   
  # list of n sequence
  listSeq <- tryCatch(parallel::parLapply(cl, 1:n, .markovchainSPHelper, t0, mclist, include.t0), 
                      error=function(e) e, warning=function(w) w)  
  
  # release the resources
  parallel::stopCluster(cl)
  
  return(listSeq)
}

######################################################################

# function to fit a DTMC with Laplacian Smoother
.mcFitLaplacianSmooth <- function(stringchar, byrow, laplacian = 0.01) {
  
  # every element of the matrix store the number of times jth state appears just
  # after the ith state
  origNum <- createSequenceMatrix(stringChar = stringchar, toRowProbs = FALSE)
  
  # add laplacian  to the sequence matrix
  # why? to avoid the cases where sum of row is zero
  newNum <- origNum + laplacian
  
  # store sum of each row  in the vector
  newSumOfRow <- rowSums(newNum)
  
  # helper matrix to convert frequency matrix to transition matrix
  newDen <- matrix(rep(newSumOfRow, length(newSumOfRow)), byrow = FALSE, ncol = length(newSumOfRow))
  
  # transition matrix
  transMatr <- newNum / newDen
  
  # create a markovchain object
  outMc <- new("markovchain", transitionMatrix = transMatr, name = "Laplacian Smooth Fit")

  # transpose the transition matrix
  if (byrow == FALSE) {
    outMc@transitionMatrix <- t(outMc@transitionMatrix)
    outMc@byrow <- FALSE
  }
  
  # wrap markovchain object in a list
  out <- list(estimate = outMc)
  return(out)
}

# function that return a Markov Chain from a given matrix of observations
.matr2Mc <- function(matrData, laplacian = 0) {
  
  # number of columns in the input matrix  
  nCols <- ncol(matrData)
  
  # an empty character vector to store names of possible states
  uniqueVals <- character()
  
  # populate uniqueVals with names of states 
  for(i in 1:nCols) {
    uniqueVals <- union(uniqueVals, unique(as.character(matrData[,i]))) 
  }
  
  # possible states in lexicographical order
  uniqueVals <- sort(uniqueVals)
  
  # create a contingency matrix which store the number of times 
  # jth state appear just after the ith state
  contingencyMatrix <- matrix(rep(0, length(uniqueVals)^2), ncol = length(uniqueVals))
  
  # set the names of rows and columns
  rownames(contingencyMatrix) <- colnames(contingencyMatrix) <- uniqueVals
  
  # fill the contingency matrix
  for (i in 1:nrow(matrData)) {
    for (j in 2:nCols) {
      # state in the ith row and (j-1)th column
      stateBegin <- as.character(matrData[i, j-1])
      
      # index of beginning state 
      whichRow <- which(uniqueVals == stateBegin)
      
      # state in the ith row and jth column
      stateEnd <- as.character(matrData[i, j])
      
      # index of ending state 
      whichCols <- which(uniqueVals == stateEnd)
      
      # update the contingency matrix
      contingencyMatrix[whichRow, whichCols] <- contingencyMatrix[whichRow, whichCols] + 1
    }
  }
  
  # add laplacian correction if needed
  contingencyMatrix <- contingencyMatrix + laplacian
  
  # take care of rows with all entries 0
  sumOfRows <- rowSums(contingencyMatrix) 
  for(i in 1:length(sumOfRows)) {
    if(sumOfRows[i] == 0) {
      contingencyMatrix[i, ] <- 1
      sumOfRows[i] <- length(sumOfRows)
    }
  }
  
  # get a transition matrix and a DTMC
  transitionMatrix <- contingencyMatrix / sumOfRows
  
  # markov chain object to be returned
  outMc <- new("markovchain", transitionMatrix = transitionMatrix)
 
  return(outMc)
}


#' @title markovchainListFit
#' 
#' @description  Given a data frame or a matrix (rows are observations, by cols 
#' the temporal sequence), it fits a non - homogeneous discrete time markov chain 
#' process (storing row)
#' 
#' @param data Either a matrix or a data.frame object.
#' @param laplacian Laplacian correction (default 0).
#' @param byrow Indicates whether distinc stochastic processes trajectiories are shown in distinct rows.
#' @param name Optional name.
#' 
#' @return A list containing two slots:
#' estimate (the estimate)
#' name
#' 
#' @examples
#' 
#' # using holson dataset
#' data(holson)
#' # fitting a single markovchain
#' singleMc <- markovchainFit(data = holson[,2:12])
#' # fitting a markovchainList
#' mclistFit <- markovchainListFit(data = holson[, 2:12], name = "holsonMcList")

markovchainListFit <- function(data, byrow = TRUE, laplacian = 0, name) {
  
  # check the format of input data
  if (! (class(data) %in% c("data.frame", "matrix"))) {
    stop("Error: data must be either a matrix or a data.frame")  
  }
  
  # if input is data frame convert it to matrix
  if(is.data.frame(data)) {
    data <- as.matrix(data)
  } 
  
  # make the entries row wise if it is not
  if(!byrow) {
    data <- t(data) 
  }
  
  # number of columns in the matrix
  nCols <- ncol(data)
  
  # allocate a list of markovchain: a non - homog DTMC process is a 
  # list of DTMC of length n-1, being n the length of the sequence
  markovchains <- list() 
  
  # fit by columns
  for(i in 2:(nCols)) { 
    
    # (i-1)th transition matrix for transition from (i-1)th state to ith state
    estMc <- .matr2Mc(matrData = data[, c(i-1, i)], laplacian = laplacian)
    
    # give name to the markovchain object which is same as the name of (i-1)th column
    if(!is.null(colnames(data))) {
      estMc@name <- colnames(data)[i-1]  
    }
    
    # store one transition matrix at every iteration
    markovchains[[i-1]] <- estMc
  }
  
  # create markovchainList object
  outMcList <- new("markovchainList", markovchains = markovchains)
  
  # wrap the object in a list
  out <- list(estimate = outMcList)
  
  # set the name of markovchainList object as given in the argument
  if(!missing(name)) {
    out$estimate@name <- name 
  }
  
  return(out)
}  

#' Return MultinomialWise Confidence intervals.
#' 
#' @description Return estimated transition matrix assuming a Multinomial Distribution
#' 
#' @param transitionMatrix An estimated transition matrix.
#' @param countsTransitionMatrix Empirical (conts) transition matrix, on which the \code{transitionMatrix} was performed.
#' @param confidencelevel confidence interval level.
#' @return Two matrices containing the confidence intervals.
#' 
#' @seealso \code{markovchainFit}
#' 
#' @references Constructing two-sided simultaneous confidence intervals 
#' for multinomial proportions for small counts in a large number of cells. 
#' Journal of Statistical Software 5(6) (2000)
#'
#' @examples 
#' myletterseq<-sample(x = letters[1:3],size = 120,replace=TRUE)
#' myMcFit<-markovchainFit(data=myletterseq)
#' myMultinomialCI=multinomialConfidenceIntervals(transitionMatrix=myMcFit$estimate,countsTransitionMatrix=createSequenceMatrix(stringChar = myletterseq))
multinomialConfidenceIntervals<-function(transitionMatrix, countsTransitionMatrix, confidencelevel=0.95) {
  
  out<-.multinomialCIRcpp(transMat=transitionMatrix, seqMat=countsTransitionMatrix,confidencelevel=confidencelevel)
  return(out)

  }
