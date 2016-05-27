#sampler for univariate markov chains
# markovchainSequence<-function(n,markovchain, t0=sample(markovchain@states,1),include.t0=FALSE)
# {
#   if(!(t0%in%markovchain@states)) stop("Error! Initial state not defined")
#   chain=character()
#   state=t0
#   for(i in 1:n) {
#     rowProbs<-markovchain@transitionMatrix[which(markovchain@states==state),]
#     outstate<-sample(size=1, x=markovchain@states, prob=rowProbs)
#     chain=c(chain, outstate)
#     state=outstate
#   }
#   if(include.t0) out<-c(t0, chain) else out<-chain
#   return(out)
# }

markovchainSequence<-function (n, markovchain, t0 = sample(markovchain@states, 1),
                               include.t0 = FALSE, useRCpp = TRUE)
{
  if (!(t0 %in% markovchain@states))
    stop("Error! Initial state not defined")
  
  # call to cpp implmentation of markovchainSequence
  if (useRCpp) {
    return(.markovchainSequenceRcpp(n, markovchain, t0, include.t0))
  }
  
  chain <- rep(NA,n)# CHANGED
  state <- t0
  for (i in 1:n) {
    rowProbs <- markovchain@transitionMatrix[which(markovchain@states == state), ]
    outstate <- sample(size = 1, x = markovchain@states,
                       prob = rowProbs)
    chain[i] <- outstate #CHANGED
    state <- outstate
  }
  if (include.t0)
    out <- c(t0, chain)
  else out <- chain
  return(out)
}


################
#random sampler#
################


#check if the subsequent states are included in the previous ones

# TODO: too strong contidion; should be changed by checking that
# all states that can be reached in one step at t-1 are named  
# in object[[t]]

.checkSequence <- function(object)
{
  out <- TRUE
  if (length(object) == 1)
    return(out) #if the size of the list is one do
  for (i in 2:length(object))
  {
    statesNm1 <- states(object[[i - 1]]) #evalutate mc n-1
    statesN <- states(object[[i]]) #evaluate mc n
    intersection <-
      intersect(statesNm1, statesN) #check the ibntersection
    if (setequal(intersection, statesNm1) == FALSE) {
      #the states at n-1
      out <- FALSE
      break
    }
  }
  return(out)
}


rmarkovchain <- function(n, object, what = "data.frame", useRCpp = TRUE, ...)
{
  if (class(object) == "markovchain")
    out <- markovchainSequence(n = n, markovchain = object, useRCpp = useRCpp, ...)
  if (class(object) == "markovchainList")
  {
    #######################################################
    if(useRCpp) {
      include.t0 <- list(...)$include.t0
      include.t0 <- ifelse(is.null(include.t0), FALSE, include.t0)
      
      dataList <- .markovchainListRcpp(n, object@markovchains, include.t0)
      
      if (what == "data.frame")
        out <- data.frame(iteration = dataList[[1]], values = dataList[[2]])
      else {
        out <- matrix(data = dataList[[2]], nrow = n, byrow = TRUE)
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
    
    verify <- .checkSequence(object = object)
    if (!verify)
      warning(
        "Warning: some states in the markovchain sequences are not contained in the following states!"
      )
    iteration <- numeric()
    values <- character()
    for (i in 1:n)
      #number of replicates
    {
      #the first iteration may include initial state
      sampledValues <-
        markovchainSequence(n = 1, markovchain = object[[1]], ...)
      outIter <- rep(i, length(sampledValues))
      
      if (length(object) > 1)
      {
        for (j in 2:length(object))
        {
          pos2take <- length(sampledValues)
          newVals <-
            markovchainSequence(n = 1,
                                markovchain = object[[j]],
                                t0 = sampledValues[pos2take]) #the initial state is in the ending position of the mclist
          outIter <- c(outIter, i)
          sampledValues <- c(sampledValues, newVals)
        }
      }
      iteration <- c(iteration, outIter)
      values <- c(values, sampledValues)
    }
    #defining the output
    if (what == "data.frame")
      out <- data.frame(iteration = iteration, values = values)
    else
      #matrix
    {
      out <- matrix(data = values,
                    nrow = n,
                    byrow = TRUE)
      if (what == 'list')
        #or list?
      {
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

# helper function
.markovchainSPHelper <- function(x) {
  # number of transition matrices
  n <- length(mclist@markovchains)
  # a character vector to store a single sequence
  seq <- character(length = n)
  
  for (i in 1:n) {
    stateNames <- mclist@markovchains[[i]]@states 
    byRow <- mclist@markovchains[[i]]@byrow
    
    if(byRow) prob <- mclist@markovchains[[i]]@transitionMatrix[which(stateNames == t0), ]
    else prob <- mclist@markovchains[[i]]@transitionMatrix[, which(stateNames == t0)]
    
    t0 <- sample(x = stateNames, size = 1,  prob = prob)
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
  verify <- .checkSequence(object = object)
  if (!verify) {
    warning("Warning: some states in the markovchain sequences are not contained in the following states!")
  }
    
  # Calculate the number of cores
  no_cores <- parallel::detectCores() - 1
  
  if((! is.null(num.cores))  && num.cores <= no_cores + 1 && num.cores >= 1) {
    no_cores <- num.cores
  }
  
  # Initiate cluster
  cl <- parallel::makeCluster(no_cores)
  
  # export the variables to be used in the helper function
  parallel::clusterExport(cl, "t0")
  
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
