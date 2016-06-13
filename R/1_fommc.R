# This function will calculate the first parameter i.e. all transition matrix
# reference : W.-K. Ching et al. / Linear Algebra and its Applications
# page no : 496
# link to pdf copy : http://hkumath.hku.hk/~imr/IMRPreprintSeries/2007/IMR2007-15.pdf

.allTransMat <- function(seqMat) {
  # number of sequence
  s <- nrow(seqMat)
  
  # size of each sequence
  n <- ncol(seqMat)
  
  # possible states
  sname <- sort(unique(as.vector(seqMat)))
  
  # number of distinct states
  m <- length(sname)
  
  # list to store all transition matrices
  lstTransMat <- list()
  
  # populate list with all transition matrices
  for (i in 1:s) {
    for (j in 1:s) {
      
      # temporary matrix
      tmat <- matrix(0, nrow = m, ncol = m)
      colnames(tmat) <- sname 
      rownames(tmat) <- sname
      
      # populate temporary matrix 
      for (k in 1:(n - 1)) {
        dat <- tmat[seqMat[i, k + 1], seqMat[j, k]]
        dat <- dat + 1
        tmat[seqMat[i, k + 1], seqMat[j, k]] <- dat
      }
      
      # calculate columns sums to get transition matrix
      csum <- colSums(tmat)
      for (l in 1:length(csum)) {
        if (csum[l] == 0) {
          tmat[, l] <- 1 / m
          csum[l] <- 1
        }
      }
      
      # convert frequency to probabilty
      # take care of rows with all entries equal to 1
      tmat <- t(t(tmat) / csum)
      
      # every time add a one step transition matrix
      lstTransMat[[s * (i - 1) + j]] <- tmat
    }
  }
  
  return(lstTransMat)
}


# calculate frequency probabilty for each sequence
# useful in finding the value of lambdas
# same reference as mentioned above

.allSeq2freqProb <- function(seqMat) {
  
  # number of sequence
  s <- nrow(seqMat)
  
  # possible states
  sname <- sort(unique(as.vector(seqMat)))
  
  # number of distinct states
  m <- length(sname)
  
  # to store state probability distribution for all sequences
  freqProbList <- list()
  
  # populate the list
  for(i in 1:s) {
    
    # create temporary matrix
    fmat <- matrix(0, ncol = 1, nrow = m)
    rownames(fmat) <- sname
    
    # populate temporary matrix
    for(j in seqMat[i, ]) {
      fmat[j, 1] <- fmat[j, 1] + 1
    }
    
    # normalize to get probabilty distribution
    fmat <- fmat/colSums(fmat)
    
    # everytime add a probabilty distribution for ith sequence
    freqProbList[[i]] <- fmat
  }
  
  return(freqProbList)
}

# objective function to pass to solnp
.fn2 <- function(params, ...) {
  tmat <- 0
  hdata <- list(...)
  
  s <- hdata$s
  w <- hdata$w
  lstTransMat <- hdata$lstTransMat
  freqProbList <- hdata$freqProbList
  
  for (i in 1:s) {
    tmat <- tmat + (params[i]*(lstTransMat[[s*(w-1) + i]] %*% freqProbList[[i]]) )
  }
  return(sum((tmat-freqProbList[[w]])^2))
}

# equality constraint function to pass to solnp
.eqn2 <- function(params, ...){
  return(sum(params))
}

#' @name fitFirstOrderMultivarMC
#' @title Function to fit First Order Multivariate Markovchain
#' @description Given a number of sequence(in the form of matrix) this function will fit the 
#'              First Order Multivariate Markovchain
#' 
#' @param seqMat Matrix whose each row is a sequence
#' 
#' @return A list of two parameters
#'         1. lambdas
#'         2. transition matrices
#'         
#' @references W.-K. Ching et al. / Linear Algebra and its Applications
#' 
#' @author Giorgio Spedicato, Deepak Yadav
#' 
#' @seealso \code{\link{markovchainFit}}
#' 
#' @examples 
#' smat <- matrix(c('a','a','b','b','a','c','a',
#'                  'b','a','b','c','c','b','a',
#'                  'a','b','b','a','c','a','b',
#'                  'c','c','c','c','a','a','a'),
#'                nrow = 4, byrow = TRUE)
#'                
#' fommcFit <- fitFirstOrderMultivarMC(smat)
#' 
#' @export
#' 
fitFirstOrderMultivarMC <- function(seqMat) {
  
  # list of transition matrices
  lstTransMat <- .allTransMat(seqMat)
  
  # list of freq probability
  freqProbList <- .allSeq2freqProb(seqMat)
  
  # number of sequences
  s <- nrow(seqMat)
  
  # which row's lambdas
  w <- 1
  
  # an s x s matrix to store all lambda values
  lambdas <- matrix(nrow = s, ncol = s)
  
  # populate lambdas matrix
  for(i in 1:s) {
    w <- i 
    params <- rep(1/s, s)
    model <- Rsolnp::solnp(params, fun = .fn2, eqfun = .eqn2, eqB = 1, LB = rep(0, s), control = list(trace = 0),
                           lstTransMat = lstTransMat, freqProbList = freqProbList, s = s, w = w)
    
    lambda <- model$pars 
    lambdas[i, ] <- lambda
  }
  
  # Aim is to return a List comprises of
  # 1. Lambda = the matrix consists of all lambda's values
  # 2. P: a list of transition matrices
  # constraint : store everything in a row-wise fashion
  
  lambdas <- t(lambdas)
  
  for(i in 1:length(lstTransMat)) {
    lstTransMat[[i]] <- t(lstTransMat[[i]])
  }
  
  return(list(Lambda = lambdas, P = lstTransMat))
}