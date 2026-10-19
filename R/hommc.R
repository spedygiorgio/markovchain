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
#'           
#'@export
# S4 class 'hommc' = High Order Multivariate Markov Chain (Ching, Ng and Fung 2008): s categorical
# sequences observed together, each depending on the last `order` values of ALL sequences.
# Slots: P = 3D array of transition matrices (one per lag and per pair of sequences),
# Lambda = their weights.
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
      # t is the index of transition matrix for transition from i sequence to j sequence
      # order of transition matrices in P is P1{1,1},P2{1,1}..Pn{1,1},P1{1,2}....Pn{s,s}
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
# Builds the n*s*s transition matrices P_h^{(i,j)}: from sequence j at lag h to sequence i now.
# At each step: pair the values x[t-h] of sequence j with the values y[t] of sequence i, count the
# transitions, turn the counts into probabilities (createSequenceMatrix) and transpose, so
# the matrices are stored by column.
.allTransMat <- function(data, order = 2) {
  n <- order # order
  uelement <- sort(unique(as.character(data))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(data) # number of categorical sequence
  lseq <- ncol(data) # length of each categorical sequence
  
  # store all transition matrices
  # 3rd dimension has n*s*s entries ordered as P1{1,1}, P2{1,1}, ..., Pn{1,1}, P1{1,2}, ...
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
# Marginal distribution of every sequence from the frequencies of its states
# (one column per sequence, stacked on the third dimension).
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

# Objective function minimised by Rsolnp::solnp. For each sequence i it measures how far the
# model's prediction  sum_j sum_k lambda_{ijk} P_k^{(ij)} x_j  is from its own distribution x_i
# (norm `Norm`), and adds up the errors. The weights `params` are the unknowns.
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
  
  # `helper` stacks the predictions of the n lags; the lag-1 block gets the weighted sum, the blocks for
  # lags 2..n get x_j itself (the fixed part of the Ching et al. formulation).
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

# Equality constraints for solnp: for each sequence i, its weights must sum to 1.
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

#' Fit a higher order multivariate Markov chain
#'
#' @description Fits the higher order multivariate Markov chain of Ching, Ng
#'   and Fung (2008) to several categorical sequences observed together: the
#'   next state of each sequence depends on the last `order` states of all the
#'   sequences, through a weighted sum of lag-specific transition matrices.
#'
#' @param seqMat A matrix or data frame with one categorical sequence per
#'   column, all of the same length; numeric or factor values are taken as
#'   state labels (compared as character strings).
#' @param order Order of the chain (default 2).
#' @param Norm Norm used to estimate the weights (default 2): the weights
#'   minimize this norm of the difference between the state frequencies of
#'   each sequence and their image under the model, as in Ching et al.
#'   (2008). The minimization needs the suggested package \pkg{Rsolnp}.
#'
#' @return An [`hommc`][hommc-class] object. If \pkg{Rsolnp} is not
#'   available, `NULL` with a warning.
#'
#' @examples
#' data <- matrix(c("2", "1", "3", "3", "4", "3", "2", "1", "3", "3", "2", "1",
#'                  "2", "4", "4", "4", "4", "2", "3", "3", "1", "4", "3", "3"),
#'                ncol = 2)
#' if (requireNamespace("Rsolnp", quietly = TRUE)) {
#'   fitHighOrderMultivarMC(data, order = 2, Norm = 2)
#' }
#'
#' @references Ching, W.-K., Ng, M. K. and Fung, E. S. (2008). Higher-order
#'   multivariate Markov chains and their applications. Linear Algebra and its
#'   Applications 428(2-3), 492-507. \doi{10.1016/j.laa.2007.05.021}
#' @family higher order chains
#' @seealso [predictHommc()] to simulate from the fit.
#' @author Giorgio Spedicato, Deepak Yadav
#'
#' @export
# Fitting: estimate the transition matrices by counting, and the weights by minimising .fn3 with
# the constraints of .eqn3 (non-negative weights summing to 1 for every sequence).
fitHighOrderMultivarMC <- function(seqMat, order = 2, Norm = 2) {
  # Rsolnp is an optional dependency: without it we warn and return NULL.
  if (requireNamespace("Rsolnp", quietly = TRUE)) {
  message("This function is experimental")
  
  # Step 1: validation. Data frame -> matrix; no missing values; valid order.
  if(is.data.frame(seqMat) == TRUE) {
    seqMat <- as.matrix(seqMat)
  }
  if (!is.matrix(seqMat) || anyNA(seqMat) || ncol(seqMat) < 1L) {
    stop("seqMat must be a matrix or data frame of states without missing values")
  }
  if (length(order) != 1L || is.na(order) || order < 1 || order != floor(order) ||
      order >= nrow(seqMat)) {
    stop("order must be a positive integer smaller than the length of the sequences")
  }
  # States are matched by name: numeric or factor codes are used as labels (not as positions).
  # states are matched by name: numeric or factor codes are used as labels
  seqMat <- matrix(as.character(seqMat), nrow = nrow(seqMat), ncol = ncol(seqMat))
  
  # Transpose so that each ROW is a sequence (internal convention).
  seqMat <- t(seqMat)
  
  # Step 2: transition matrices and marginal distributions from the data.
  # array of transition matrices
  allTmat <- .allTransMat(seqMat, order = order)
  
  # array of freq probability
  freqMat <- .allFreqProbMat(seqMat)
  
  n <- order # order
  uelement <- sort(unique(as.character(seqMat))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(seqMat) # number of categorical sequence
  
  # Step 3: start from equal weights and run the constrained optimiser.
  lmbda <- rep(1 / (n * s), n * s * s)
  
  fit <- Rsolnp::solnp(pars = lmbda, fun =  .fn3, eqfun = .eqn3, eqB = rep(1, s),
                       LB = rep(0, n * s * s), control = list(trace = 0), 
                       allTmat = allTmat, freqMat = freqMat, n = n, m = m, s = s, Norm = Norm)
  return(new("hommc", order = order, Lambda = fit$pars, P = allTmat, states = uelement, byrow = FALSE))
  } else {
    warning("package Rsolnp is required by fitHighOrderMultivarMC() and is not available: NULL is returned",
            call. = FALSE)
    return(NULL)
  }
}


#' Simulate from a higher order multivariate Markov chain
#'
#' @description Simulates the next `t` states of each sequence of a higher
#'   order multivariate Markov chain, given their last `order` states.
#'
#' @param hommc An [`hommc`][hommc-class] object, for instance fitted by
#'   [fitHighOrderMultivarMC()].
#' @param t Number of steps to simulate, a positive integer.
#' @param init A matrix with one row per sequence and `order` columns, the
#'   last states of each sequence, the most recent in the last column. If
#'   missing, every sequence starts from the first state, repeated.
#'
#' @details At each step the next state of sequence \eqn{j} is drawn from
#'   \eqn{\sum_k \sum_h \lambda_{jkh} P_h^{(jk)} x^{(k)}_{t-h+1}}, where
#'   \eqn{x^{(k)}_{t-h+1}} is the state of sequence \eqn{k} \eqn{h - 1}
#'   steps before the current one (Ching et al., 2008). The matrices are read
#'   by column (`P[to, from]`), as returned by [fitHighOrderMultivarMC()], or
#'   by row when the `byrow` slot of `hommc` is `TRUE`; all the sequences are
#'   drawn from the same past before it is updated.
#'
#' @return A character matrix with one row per sequence and `t` columns, the
#'   simulated states.
#'
#' @references Ching, W.-K., Ng, M. K. and Fung, E. S. (2008). Higher-order
#'   multivariate Markov chains and their applications. Linear Algebra and its
#'   Applications 428(2-3), 492-507. \doi{10.1016/j.laa.2007.05.021}
#' @author Vandit Jain
#' @family higher order chains
#' @examples
#' data <- matrix(c("2", "1", "3", "3", "4", "3", "2", "1", "3", "3", "2", "1",
#'                  "2", "4", "4", "4", "4", "2", "3", "3", "1", "4", "3", "3"),
#'                ncol = 2)
#' if (requireNamespace("Rsolnp", quietly = TRUE)) {
#'   fit <- fitHighOrderMultivarMC(data, order = 2, Norm = 2)
#'   # the last two states of each sequence
#'   init <- t(data[11:12, ])
#'   predictHommc(fit, 5, init)
#' }
#' @export
# Predicts (simulates) t future values of every sequence, starting from the last `order` values
# of each sequence given in `init`.
predictHommc <- function(hommc, t, init) {
  if(!is(hommc, "hommc")) {
    stop("Please provide a valid hommc-class object")
  }

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
  init <- as.matrix(init)
  
  if(!all(dim(init) == c(s,n))){
    stop("Please provide sufficient number of previous states")
  }
  
  if(length(t) != 1L || is.na(t) || t <= 0 || t != floor(t))
    stop("T should be a positive integer")
  
  if(!all(init %in% states))
    stop("invalid states in provided state matrix init")
  
  ## initialize result matrix
  result <- matrix(NA,nrow = s,ncol = t)
  
  # Model: x_{t+1}^{(j)} = sum_k sum_h Lambda_{jkh} P_h^{(jk)} x_{t-h+1}^{(k)} gives the PROBABILITIES of
  # the next state of sequence j; we then draw one state at random.
  ## The model (Ching et al.) is
  ##   x_{t+1}^{(j)} = sum_k sum_h Lambda_{jkh} P_h^{(jk)} x_{t-h+1}^{(k)},
  ## where P_h^{(jk)} maps the state of sequence k (the "from" state) h steps
  ## before to the next state of sequence j. The matrices are stored by column
  ## (P[to, from]) unless byrow = TRUE, and init[k, n] is the most recent state
  ## of sequence k. All s sequences are drawn from the same past before the
  ## past is updated.
  # Loop over the t time steps to predict.
  for(i in 1:t)
  {
    current <- character(s)
    for(j in 1:s)
    {
      ## initialises probability according
      probs <- rep(0,m)
      
      ## probability for current sequence depends all sequence
      # Every sequence k contributes, at every lag h, with its weight times the column of the
      # transition matrix that corresponds to the state k had h steps ago.
      for(k in 1:s)
      {
        ## gets index of coressponding in the 3-D array P 
        # index is the index of transition matrix for transition from sequence k to sequence j
        # order of transition matrices in P is P1{1,1},P2{1,1}..Pn{1,1},P1{1,2}....Pn{s,s}
        index <- n * s * (j-1) + n * (k-1)
        
        ## iterates for all order 1 to n
        for(h in 1:n)
        {
          prev <- init[k,n-h+1]
          label <- which(prev == states)
          column <- if (isTRUE(hommc@byrow)) hommc@P[label, , h + index] else hommc@P[, label, h + index]
          probs <- probs + hommc@Lambda[h + index] * column
        }
      }
      
      ## uses sample function from base package
      current[j] <- sample(size = 1, x = states, prob = probs)
    }
    
    # Slide the window: forget the oldest column of states and append the newly drawn ones.
    ## changes init for next t iteration
    if (n > 1L) init[, 1:(n - 1L)] <- init[, 2:n, drop = FALSE]
    init[, n] <- current
    result[, i] <- current
  }
  
  ## returns result
  return(result)
}
