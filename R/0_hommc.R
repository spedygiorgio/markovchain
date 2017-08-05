#' An S4 class for representing High Order Multivariate Markovchain (HOMMC)
#' 
#' @slot order an integer equal to order of Multivariate Markovchain
#' @slot states a vector of states present in the HOMMC model
#' @slot P array of transition matrices
#' @slot Lambda a vector which stores the weightage of each transition matrices in P
#' @slot byrow if FALSE each column sum of transition matrix is 1 else row sum = 1
#' @slot name a name given to hommc
#' 
#' @author Giorgio Spedicato, Deepak Yadav
#' 
#' @examples 
#' statesName <- c("a", "b")
#' 
#' P <- array(0, dim = c(2, 2, 4), dimnames = list(statesName, statesName))
#' P[,,1] <- matrix(c(0, 1, 1/3, 2/3), byrow = FALSE, nrow = 2)
#' P[,,2] <- matrix(c(1/4, 3/4, 0, 1), byrow = FALSE, nrow = 2)
#' P[,,3] <- matrix(c(1, 0, 1/3, 2/3), byrow = FALSE, nrow = 2)
#' P[,,4] <- matrix(c(3/4, 1/4, 0, 1), byrow = FALSE, nrow = 2)
#' 
#' Lambda <- c(0.8, 0.2, 0.3, 0.7)
#' 
#' ob <- new("hommc", order = 1, states = statesName, P = P, 
#'           Lambda = Lambda, byrow = FALSE, name = "FOMMC")

hommc <- setClass("hommc",
                      slots = list(order = "numeric", states  =  "character",
                      P = "array", Lambda = "numeric", byrow = "logical",
                      name = "character")
)

# internal method to show hommc object in informative way
.showHommc <- function(object) {
  
  # whether data in transition matrices are stored in column-wise or row-wise fashion
  if(object@byrow == TRUE) {
    direction <- "(by rows)" 
  } else {
    direction <- "(by cols)" 
  }
  
  # display order and unique states
  cat("Order of multivariate markov chain =", object@order, "\n")
  cat("states =", object@states, "\n")
  
  cat("\n")
  cat("List of Lambda's and the corresponding transition matrix", direction,":\n")
  
  # display transition matrices and the corresponding lambdas
  n <- object@order
  s <- sqrt((dim(object@P))[3]/n)
  
  for(i in 1:s) {
    for(j in 1:s) {
      t <- n * s * (i-1) + (j-1) * n
      for(k in 1:n) {
        cat("Lambda", k, "(", i, ",", j, ") : ", object@Lambda[t+k],"\n", sep = "")
        cat("P", k, "(", i, ",", j, ") : \n", sep = "")
        print(object@P[, , t+k])
        cat("\n")
      }  
    }
  }
}

#' @title Function to display the details of hommc object
#' @description This is a convenience function to display the slots of hommc object
#'              in proper format
#' 
#' @param object An object of class hommc
#' 
#' @rdname hommc-show
#' @export                             
setMethod("show", "hommc",
          function(object){
            .showHommc(object)
          }
)

# all transition matrices
# n*s*s n = order s = number of categorical sequences
# verified using two examples from research paper
.allTransMat <- function(data, order = 2) {
  n <- order # order
  uelement <- sort(unique(as.character(data))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(data) # number of categorical sequence
  lseq <- ncol(data) # length of each categorical sequence
  
  # store all transition matrices
  allTmat <- array(dim = c(length(uelement), length(uelement), n*s*s), dimnames = list(uelement, uelement))
  
  t <- 1 # help
  
  for(i in 1:s) {
    for(j in 1:s) {
      x <- data[j, ] # jth sequence
      y <- data[i, ] # ith sequence
      
      # jumps
      for(h in 1:n) {
        # column wise
        allTmat[ , , t] <- t(createSequenceMatrix(matrix(c(x[1:(lseq-h)], y[-(1:h)]),
                            ncol = 2, byrow = FALSE), toRowProbs = TRUE, 
                            possibleStates = uelement, sanitize = TRUE))
        t <- t + 1
      }
    }
  }
  return(allTmat)
}

# distribution of each categorical sequence based on the frequency
# verified using two examples from research paper
.allFreqProbMat <- function(data) {
  
  uelement <- sort(unique(as.character(data))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(data) # number of categorical sequence
  
  # frequency based probability for all sequences
  freqMat <- array(0, dim = c(m, 1, s), dimnames = list(uelement)) 
  
  for(i in 1:s) {
    idata <- data[i, ] # ith categorical sequence
    
    # populate frequency matrix
    for(j in idata) {
      freqMat[j, 1, i] <- freqMat[j, 1, i] + 1
    }
    
    # normalization
    freqMat[, , i] <- freqMat[, , i] / sum(freqMat[, , i])
  }
  
  return(freqMat)
}

# objective function to pass to solnp
.fn3 <- function(params, ...) {
  
  hdata <- list(...)
  
  # calculate error
  error <- 0
  
  # number of categorical sequence
  s <- hdata$s
  
  # order
  n <- hdata$n
  
  # number of uniq states || dimension of t-matrix
  m <- hdata$m
  
  # array of transition matrices
  allTmat <- hdata$allTmat
  
  # all frequency matrix
  freqMat <- hdata$freqMat
  
  # norm
  Norm <- hdata$Norm
  
  for(i in 1:s) {
    helper <- matrix(0, nrow = m*n, ncol = 1)
    for(j in 1:s) {
      helper2 <- matrix(0, nrow = m, ncol = 1)
      y <- n * (j - 1 + s * (i - 1))
      
      for(k in 1:n) {
        helper2 <- helper2 + params[y + k] * (allTmat[ , , y + k] %*% matrix(freqMat[ , , j]))
      }
      
      helper[1:m, ] <- helper[1:m, ] + helper2
      
      if(i == j && n>= 2) {
        for(k in 2:n) {
          p <- (k - 1) * m
          helper[(p + 1):(p + m)] <- freqMat[ , , j]
        }
      }
    }
    error <- error + sum(abs((helper - freqMat[ , , i]) ^ Norm))
  }
  
  return(error ^ (1 / Norm))
}

# equality constraint function to pass to solnp
.eqn3 <- function(params, ...) {
  
  hdata <- list(...)
  
  # number of categorical sequence
  s <- hdata$s
  
  # order
  n <- hdata$n
  
  toReturn <- numeric()
  
  for(i in 1:s) {
    toReturn[i] <- sum(params[((i - 1) * n * s + 1):(i * n * s)])
  }
  return(toReturn)
}

#' Function to fit Higher Order Multivariate Markov chain
#'
#' @description Given a matrix of categorical sequences it fits 
#'              Higher Order Multivariate Markov chain.
#'
#' @param seqMat a matrix or a data frame where each column 
#'               is a categorical sequence
#' @param order Multivariate Markov chain order. Default is 2.
#' @param Norm Norm to be used. Default is 2.
#' 
#' @return an hommc object
#' 
#' @examples 
#' data <- matrix(c('2', '1', '3', '3', '4', '3', '2', '1', '3', '3', '2', '1', 
#'                c('2', '4', '4', '4', '4', '2', '3', '3', '1', '4', '3', '3')), 
#'                ncol = 2, byrow = FALSE)
#'                
#' fitHighOrderMultivarMC(data, order = 2, Norm = 2)                
#' 
#' @references W.-K. Ching et al. / Linear Algebra and its Applications
#' 
#' @author Giorgio Spedicato, Deepak Yadav
#' 
#' @export
fitHighOrderMultivarMC <- function(seqMat, order = 2, Norm = 2) {
  
  message("This function is experimental")
  
  if(class(seqMat) == "data.frame") {
    seqMat <- as.matrix(seqMat)
  }
  
  seqMat <- t(seqMat)
  
  # array of transition matrices
  allTmat <- .allTransMat(seqMat, order = order)
  
  # array of freq probability
  freqMat <- .allFreqProbMat(seqMat)
  
  n <- order # order
  uelement <- sort(unique(as.character(seqMat))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(seqMat) # number of categorical sequence
  
  lmbda <- rep(1 / (n * s), n * s * s)
  
  fit <- Rsolnp::solnp(pars = lmbda, fun =  .fn3, eqfun = .eqn3, eqB = rep(1, s),
                       LB = rep(0, n * s * s), control = list(trace = 0), 
                       allTmat = allTmat, freqMat = freqMat, n = n, m = m, s = s, Norm = Norm)
  
  
  return(new("hommc", order = order, Lambda = fit$pars, P = allTmat, states = uelement, byrow = FALSE))
}


#' Simulate a higher order multivariate markovchain
#' 
#' @description 
#' This function provides a prediction of states for a higher order 
#' multivariate markovchain object
#' 
#' @usage predictHommc(hommc,t,init)
#' 
#' @param hommc a hommc-class object
#' @param t no of iterations to predict
#' @param init matrix of previous states size of which depends on hommc
#' 
#' @details 
#' The user is required to provide a matrix of giving n previous coressponding
#' every categorical sequence. Dimensions of the init are s X n, where s is 
#' number of categorical sequences and n is order of the homc.
#' 
#' @return 
#' The function returns a matrix of size s X t displaying t predicted states 
#' in each row coressponding to every categorical sequence. 
#' 
#' @author Vandit Jain
#' 
#' 
#' @export
predictHommc <- function(hommc, t, init) {
  ## order of markovchain 
  n <- hommc@order
  
  ## number of categorical sequences
  s <- sqrt((dim(hommc@P))[3]/n)
  
  ## list of states
  states <- hommc@states
  
  ## size of set of all possible states
  m <- length(states)
  
  ## if initial states not provided take statndard example
  if(missing(init)) {
    init <- matrix(rep(states[1],s*n),nrow = s,byrow = TRUE)
  }
  
  if(dim(init)!= c(s,n)){
    stop("Please provide sufficient number of previous states")
  }
  
  if(class(hommc)!= "hommc") {
    stop("Please provide a valid hommc-class object")
  }
  
  if(t <=0)
    stop("T should be a positive integer")
  
  for(i in 1:s)
  {
    for(j in 1:n)
    {
      if(!(init[i,j] %in% states))
        stop("invalid states in provided state matrix init")
    }
  }
  
  ## initialize result matrix
  result <- matrix(NA,nrow = s,ncol = t)
  
  ## runs loop according to hommc class structure
  for(i in 1:t)
  {
    for(j in 1:s)
    {
      ## initialises probability according
      rowProbs <- rep(0,m)
      
      ## probability for current sequence depends all sequence
      for(k in 1:s)
      {
        ## gets index of coressponding in the 3-D array P 
        index <- n * s * (j-1) + n * (k-1)
        
        ## iterates for all order 1 to n
        for(h in 1:n)
        {
          prev <- init[j,n-h+1]
          label <- which(prev == states)
          rowProbs <- rowProbs + hommc@Lambda[h + index] * hommc@P[label, ,h + index]
        }
      }
      
      ## uses sample function from base package
      curr <- sample(size = 1, x = states, prob = rowProbs)
      
      ## changes init for next t iteration
      for(temp in 2:n)
      {
        if(temp <= n)
        init[j,temp-1] = init[j,temp]
      }
      init[j,n] = curr;
      result[j,i] = curr;
    }
  }
  
  ## returns result
  return(result)
}




