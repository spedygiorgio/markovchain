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
    
    # possible states in the previous markovchain object
    statesNm1 <- states(object[[i - 1]])
    
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

rmarkovchain <- function(n, object, what = "data.frame", useRCpp = TRUE, ...) {
  
  # check the class of the object
  if (class(object) == "markovchain") {
    out <- markovchainSequence(n = n, markovchain = object, useRCpp = useRCpp, ...)
    return(out)
  }
    
  if (class(object) == "markovchainList")
  {
    #######################################################
    if(useRCpp) {
      
      # if include.t0 is not passed as extra argument then set include.t0 as false
      include.t0 <- list(...)$include.t0
      include.t0 <- ifelse(is.null(include.t0), FALSE, include.t0)
      
      # call fast cpp function
      dataList <- .markovchainListRcpp(n, object@markovchains, include.t0)
      
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
      
      #the first iteration may include initial state
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
.markovchainSPHelper <- function(x) {
  # number of transition matrices
  n <- length(mclist@markovchains)
  # a character vector to store a single sequence
  seq <- character(length = n)
  
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
    seq[i] <- t0
  }
  
  return(seq)
}

#' Function to generate a list of sequence of states in parallel from non-homogeneous Markov chains.
#' 
#' Provided any markovchainList object, it returns a list of sequence of states coming 
#' from the underlying stationary distribution. 
#' 
#' @param  n Sample size
#' @param object markovchainList object
#' @param t0 Initial state
#' @param num.cores Number of cores
#'   
#' @export

markovchainSequenceParallel <- function(n, object,
                                        t0 = sample(object@markovchains[[1]]@states, 1),
                                        num.cores = NULL) {
  # check for the validity of non-uniform markov chain
  verify <- .checkSequence(object@markovchains)
  if (!verify) {
    warning("Warning: some states in the markovchain sequences are not contained in the following states!")
  }
    
  # Calculate the number of cores
  # It's not good to use all cores
  no_cores <- parallel::detectCores() - 1
  
  # number of cores specified should be less than or equal to maximum cores available
  if((! is.null(num.cores))  && num.cores <= no_cores + 1 && num.cores >= 1) {
    no_cores <- num.cores
  }
  
  # Initiate cluster
  cl <- parallel::makeCluster(no_cores)
  
  # export the variables to be used in the helper function
  parallel::clusterExport(cl, "t0")
  
  # export the variables to be used in the helper function
  mclist <- object
  parallel::clusterExport(cl, "mclist")
   
  # list of n sequence
  listSeq <- tryCatch(parallel::parLapply(cl, 1:n, .markovchainSPHelper), 
                      error=function(e) e, warning=function(w) w)  
  
  # release the resources
  parallel::stopCluster(cl)
  
  return(listSeq)
}

######################################################################

#core function to get sequence matrix

#createSequenceMatrix <- function(stringchar, toRowProbs = FALSE, sanitize = TRUE) {
#    .Call('markovchain_createSequenceMatrix', PACKAGE = 'markovchain', stringchar, toRowProbs, sanitize)
#}

#functon to fit a Markov chain by MLE

# .mcFitMle<-function(stringchar,byrow)
# {
#   initialMatr<-createSequenceMatrix(stringchar=stringchar,toRowProbs=TRUE)
#   outMc<-new("markovchain", transitionMatrix=initialMatr,name="MLE Fit")
#   if(byrow==FALSE) outMc<-t(outMc)
#   out<-list(estimate=outMc)
#   return(out)
# }


#function to fit a DTMC with Laplacian Smoother


.mcFitLaplacianSmooth <- function(stringchar, byrow, laplacian = 0.01)
{
  origNum <-
    createSequenceMatrix(stringchar = stringchar, toRowProbs = FALSE)
  newNum <- origNum + laplacian
  newSumOfRow <- rowSums(newNum)
  newDen <-
    matrix(rep(newSumOfRow, length(newSumOfRow)),
           byrow = FALSE,
           ncol = length(newSumOfRow))
  transMatr <- newNum / newDen
  outMc <-
    new("markovchain",
        transitionMatrix = transMatr,
        name = "Laplacian Smooth Fit")

  if (byrow == FALSE)
    outMc@transitionMatrix <- t(outMc@transitionMatrix)
  
  out <- list(estimate = outMc)
  return(out)
}


#example

#data=markovchainSequence(10000,markovA,t0="a")
#ciao<-markovchainFit(data=data)


#given a sting of characters, returns the associate one step transition matrix
# .bootstrapCharacterSequences<-function(stringchar, n, size=length(stringchar))
# {
#   contingencyMatrix<-createSequenceMatrix(stringchar=stringchar)
#   samples<-list()
#   itemset<-rownames(contingencyMatrix)
#   for(i in 1:n) #cicle to fill the samples
#   {
#     charseq<-character()
#     char<-sample(x=itemset,size=1)
#     charseq<-c(charseq,char)
#     for(j in 2:size) #cicle to define the element in a list
#     {
#       probsVector<-contingencyMatrix[which(rownames(contingencyMatrix)==char),]
#       char<-sample(x=itemset,size=1, replace=TRUE,prob=probsVector)
#       charseq<-c(charseq,char)
#     }
#     samples[[length(samples)+1]]<-charseq #increase the list
#   }
#   return(samples)
# }


# .fromBoot2Estimate<-function(listMatr)
# {
#   sampleSize<-length(listMatr)
#   matrDim<-nrow(listMatr[[1]])
#   #create the estimate output
#   matrMean<-zeros(matrDim)
#   matrSd<-zeros(matrDim)
#   #create the sample output
#   for(i in 1:matrDim) #move by row
#   {
#     for(j in 1:matrDim) #move by cols
#     {
#       probsEstimated<-numeric()
#       #fill the probs
#       for(k in 1:sampleSize) probsEstimated<-c(probsEstimated,listMatr[[k]][i,j])
#       muCell<-mean(probsEstimated)
#       sdCell<-sd(probsEstimated)
#       matrMean[i,j]<-muCell
#       matrSd[i,j]<-sdCell
#     }
#   }
# 	out<-list(estMu=matrMean, estSigma=matrSd)
#     return(out)
# }
# 
# 
# .mcFitBootStrap<-function(data, nboot=10,byrow=TRUE, parallel=FALSE)
# {
#   #create the list of bootstrap sequence sample
# 	theList<-.bootstrapCharacterSequences(stringchar=data, n=nboot)
# 	if(!parallel)
# 		#convert the list in a probability matrix
# 		pmsBootStrapped<-lapply(X=theList, FUN=createSequenceMatrix, toRowProbs=TRUE,sanitize=TRUE)
# 		 else {
# 		#require(parallel)
# 		type <- if(exists("mcfork", mode = "function")) "FORK" else "PSOCK"
# 		cores <- getOption("mc.cores", detectCores())
# 		cl <- makeCluster(cores, type = type)
# 			clusterExport(cl, varlist = c("createSequenceMatrix","zeros"))
# 			pmsBootStrapped<-parLapply(cl=cl, X=theList, fun="createSequenceMatrix", toRowProbs=TRUE,sanitize=TRUE)
# 		stopCluster(cl)
# 	}
#  
#   estimateList<-.fromBoot2Estimate(listMatr=pmsBootStrapped)
#   #from raw to estimate
#   temp<-estimateList$estMu
#   transMatr<-sweep(temp, 1, rowSums(temp), FUN="/")
#   estimate<-new("markovchain",transitionMatrix=transMatr, byrow=byrow, name="BootStrap Estimate")
#   out<-list(estimate=estimate, standardError=estimateList$estSigma,bootStrapSamples=pmsBootStrapped)
#   return(out)
# }

###############################################
#special function for matrices and data.frames#
###############################################

#function that return a Markov Chain from a given matrix of observations

.matr2Mc<-function(matrData,laplacian=0) {
  #find unique values scanning the matrix
  nCols<-ncol(matrData)
  uniqueVals<-character()
  for(i in 1:nCols) uniqueVals<-union(uniqueVals,unique(as.character(matrData[,i])))
  uniqueVals<-sort(uniqueVals)
  #create a contingency matrix
  contingencyMatrix<-matrix(rep(0,length(uniqueVals)^2),ncol=length(uniqueVals))
  rownames(contingencyMatrix)<-colnames(contingencyMatrix)<-uniqueVals
  #fill the contingency matrix
  for(i in 1:nrow(matrData))
  {
    for( j in 2:nCols)
    {
      stateBegin<-as.character(matrData[i,j-1]);whichRow<-which(uniqueVals==stateBegin)
      stateEnd<-as.character(matrData[i,j]);whichCols<-which(uniqueVals==stateEnd)
      contingencyMatrix[whichRow,whichCols]<-contingencyMatrix[whichRow,whichCols]+1
    }
  }
  #add laplacian correction if needed
  contingencyMatrix<-contingencyMatrix+laplacian
  #get a transition matrix and a DTMC
  transitionMatrix<-contingencyMatrix/rowSums(contingencyMatrix)
  outMc<-new("markovchain",transitionMatrix=transitionMatrix)
 
  return(outMc)
}

#markovchainFit <- function(data, method = "mle", byrow = TRUE, nboot = 10L, laplacian = 0, name = "", parallel = FALSE, confidencelevel = 0.95) {
#    .Call('markovchain_markovchainFit', PACKAGE = 'markovchain', data, method, byrow, nboot, laplacian, name, parallel, confidencelevel)
#}

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
#' #using holson dataset
#' data(holson)
#' #fitting a single markovchain
#' singleMc<-markovchainFit(data=holson[,2:12])
#' #fitting a markovchainList
#' mclistFit<-markovchainListFit(data=holson[,2:12],name="holsonMcList")

markovchainListFit<-function(data,byrow=TRUE, laplacian=0, name) {
  if (!(class(data) %in% c("data.frame","matrix"))) stop("Error: data must be either a matrix or a data.frame")
  if(is.data.frame(data)) data<-as.matrix(data)
  #byrow assumes distinct observations (trajectiories) are per row
  #otherwise transpose
  if(!byrow) data<-t(data)
  nCols<-ncol(data)
  #allocate a the list of markovchain: a non - homog DTMC process is a 
  #list of DTMC of length n-1, being n the length of the sequence
  markovchains<-list() #### allocating #####
  #fit by cols
  for(i in 2:(nCols)) { #if whe have at least two transitions
    # estimate the thransition for period i-1 using cols i-1 and i
    estMc<-.matr2Mc(matrData = data[,c(i-1,i)],laplacian = laplacian)
    if(!is.null(colnames(data))) estMc@name<-colnames(data)[i-1]
    markovchains[[i-1]]<-estMc  #### save this in mc #####
  }
  #return the fitted list
  outMcList<-new("markovchainList",markovchains=markovchains)
  out<-list(estimate=outMcList)
  if(!missing(name)) out$estimate@name<-name
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
#' myMultinomialCI=multinomialConfidenceIntervals(transitionMatrix=myMcFit$estimate,countsTransitionMatrix=createSequenceMatrix(stringchar = myletterseq))
multinomialConfidenceIntervals<-function(transitionMatrix, countsTransitionMatrix, confidencelevel=0.95) {
  
  out<-.multinomialCIRcpp(transMat=transitionMatrix, seqMat=countsTransitionMatrix,confidencelevel=confidencelevel)
  return(out)

  }
