#' @name verifyMarkovProperty
#' 
#' @rdname verifyMarkovProperty
#' 
#' @title Various functions to perform statistical inference of DTMC
#' @description These functions verify the Markov property, assess 
#'              the order and stationarity of the Markov chain.
#' 
#' @param sequence An empirical sequence.
#' @param ... Parameters for chi-square test.
#' @param hypothetic A transition matrix for a hypothetic markov chain sequence.
#' @param nblocks Number of blocks.
#' 
#' @return Verification result
#' 
#' @references Monika, Anderson and Goodman.
#' 
#' @author Tae Seung Kang, Giorgio Alfredo Spedicato
#' 
#' @seealso \code{markovchain}
#' 
#' @examples 
#' sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b",
#'               "a", "b", "a", "a", "b", "b", "b", "a")
#' mcFit <- markovchainFit(data = sequence, byrow = FALSE)
#' verifyMarkovProperty(sequence)
#' assessOrder(sequence)
#' assessStationarity(sequence, 1)
#' divergenceTest(sequence, mcFit$estimate@transitionMatrix)
#' 
NULL

#' @rdname verifyMarkovProperty
#' 
#' @export

# check if the sequence holds the Markov property
verifyMarkovProperty <- function(sequence, ...) {
  n <- length(sequence)
  u <- unique(sequence)
  stateNames <- u
  nelements <- length(stateNames)
  mat <- matlab::zeros(nrow = nelements, ncol = 3)
  
  # SSO: state sequence occurrences
  # TSO: two state occurences
  dimnames(mat) <- list(stateNames, c("SSO", "TSO", "TSO-SSO"))
  
  # numeric vector initialized with zero for all states
  SSO <- numeric()
  for(i in 1:nelements) {
    sname <- stateNames[i]
    SSO[sname] <- 0
  }
  
  # numeric vector initialized with zero for all states
  TSO <- SSO
  
  # store the output to be returned
  out <- list()
  
  for(present in stateNames) {
    for(future in stateNames) {
      
      for(i in 1:nelements) {
        TSO[i] <- SSO[i] <- 0  
      }
      
      # populate TSO and SSO vector
      for(i in 1:(n-1)) {
        # let the ith state as past state
        past <- sequence[i]
        
        # if next state in the sequence is present state
        if(sequence[i+1] == present) {
          TSO[past] <- TSO[past] + 1
          
          # if next to next state in the sequence is future state
          if((i < n - 1) && (sequence[i+2] == future)) {
            SSO[past] <- SSO[past] + 1
          }
        }
      }
      
      # populate the matrix
      # first column corresponds to SSO, second to TSO and
      # third to their difference
      
      for(i in 1:(length(SSO))) {
        mat[i, 1] <- SSO[i]
        mat[i, 2] <- TSO[i]
        mat[i, 3] <- TSO[i] - SSO[i]
      }
      
      # chi-squared test
      
      # between SSO and TSO-SSO
      table <- as.data.frame(mat[, c(1, 3)])
      
      # an object of class htest
      res <- chisq.test(table, ...)
      
      # extract all information from htest object
      # and stored the result in the form of list
      res <- c(res)
      
      # SSO and TSO
      table <- as.data.frame(mat[ , c(1, 2)])
      
      # stored the table in the list
      res[["table"]] <- table
      
      # store the result corresponding to present state and future state
      out[[paste0(present, future)]] <- res
    }
  }
  
  return(out)
}

#the out should be a function with following slots:
#stateTransitionSequenceTable (as in page 25)
#statistic: the chi - square statistic 
#p-value: the p-value of the statistic with appropriate degrees of freedom (see p 25-27 on the calculation)

#also the following should run: data(blanden); myMc<-as(blanden,"markovchain");sequenza<-rmarkovchain(n = 100,myMc)
#verifyMarkovProperty(sequenza)
#http://stats.stackexchange.com/questions/37386/check-memoryless-property-of-a-markov-chain 

#' @rdname verifyMarkovProperty
#' @export

# check if sequence is of first order or of second order
assessOrder <- function(sequence) {
  
  # length of sequence
  n <- length(sequence)
  
  # unique states
  states <- unique(sequence)
  
  # number of unique states
  nelements <- length(states)
  
  TStat <- 0
  for(present in states) {
    # going to be a transition matrix 
    mat <- matlab::zeros(nelements)
    dimnames(mat) <- list(states, states)
    
    # populate transition matrix
    for(i in 1:(n - 2)) {
      if(present == sequence[i + 1]) {
        past <- sequence[i]
        future <- sequence[i+2]
        mat[past, future] <- mat[past, future] + 1 
      }
    }
    
    # chi-squared test
    res <- chisq.test(mat)
    TStat <- TStat + res$statistic
  }
  
  k <- nelements
  df <- k*(k-1)^2
  pvalue <- 1-pchisq(q = TStat, df)
  
  # returning the output
  cat("The assessOrder test statistic is: ", TStat, " the Chi-Square d.f. are: ", df, " the p-value is: ", pvalue, "\n")
  out <- list(statistic = TStat[[1]], p.value = pvalue[[1]])
  
  return(out)
}

#' @rdname verifyMarkovProperty
#' @export

# check if sequence is stationary
assessStationarity <- function(sequence, nblocks) {
  
  # length of sequence
  n <- length(sequence)
  
  # size of each block
  blocksize <- n / nblocks
  
  # vector of unique states
  states <- unique(sequence)
  
  # number of states
  nstates <- length(states) 
  
  # sum of the statistics
  TStat <- 0 
  
  # chi-squared test for each state
  for(i in states) {
    
    # init matrix
    mat <- matlab::zeros(nblocks, nstates)
    dimnames(mat) <- list(1:nblocks, states)
    
    # compute the transition matrix from sequence
    for(j in 1:(n-1)) {
      if(sequence[j] == i) {
        # row index
        b <- ceiling(j / blocksize) 
        
        # next state
        future <- sequence[j+1] 
        
        # update transition matrix
        mat[b, future] <- mat[b, future] + 1 
      }
    }
    
    # vector to store row sum of matrix
    rowsums <- rowSums(mat)
    
    # store the indices with zero row sum
    indices <- which(rowsums == 0) 
    
    # update rows with zero sum
    for(k in indices) mat[k, ] <- 1/nstates 
    
    # update row sum after checking zero sum row
    rowsums <- rowSums(mat)
    
    # row-wise normalize. 
    mat <- mat/rowsums 
    
    # Some columns may still be all zeros. This causes NaN for chi-squared test.
    # chi-squared test
    res <- chisq.test(mat)
    TStat <- TStat + res$statistic
  }
  
  k <- nstates
  
  # degree of freedom
  df <- k * (nblocks - 1) * (k-1)
  pvalue <- 1-pchisq(TStat, df)
  
  # returning the output
  cat("The assessStationarity test statistic is: ", TStat, " the Chi-Square d.f. are: ", df," the p-value is: ", pvalue,"\n")
  out <- list(statistic = TStat[[1]], p.value = pvalue[[1]])
  
  return(out)
}

# sequence to transition frequencey matrix
.seq2mat <- function(sequence) {
  
  # basic requirement to create transition matrix
  n <- length(sequence)
  states <- unique(sequence)
  nstates <- length(states)
  
  # create transition matrix
  mat <- matlab::zeros(nstates)
  dimnames(mat) <- list(states, states)
  
  # populate transition matrix
  for(i in 1:(n-1)) {
    from <- sequence[i]
    to <- sequence[i+1]
    mat[from, to] <- mat[from, to] + 1
  }
  
  return (mat)
}


#' @rdname verifyMarkovProperty
#' @export




# divergence test for the hypothesized one and an empirical transition matrix from sequence
divergenceTest <- function(sequence, hypothetic) {
  # length of sequence
  n <- length(sequence)
  # frequency matrix
  empirical <- .seq2mat(sequence)
  # unique states as per the sequence given
  M <- nrow(empirical)
  v <- numeric()
  out <- 2 * n / .phi2(1)
  sum <- 0
  c <- 0
  for (i in 1:M) {
    sum2 <- 0
    sum3 <- 0
    for (j in 1:M) {
      if (hypothetic[i, j] > 0) {
        c <- c + 1
      }
      sum2 <-
        sum2 + hypothetic[i, j] * .phi(empirical[i, j] / hypothetic[i, j])
      for (k in 2:n) {
        if (sequence[k - 1] == rownames(empirical)[i] &&
            sequence[k] == rownames(empirical)[j]) {
          sum3 <- sum3 + 1
        }
      }
    }
    v[i] <- sum3
    sum <- sum + ((v[i] / n) * sum2)
  }
  TStat <- out * sum
  pvalue <- 1 - pchisq(TStat, c - M)
  cat(
    "The Divergence test statistic is: ",
    TStat,
    " the Chi-Square d.f. are: ",
    c - M,
    " the p-value is: ",
    pvalue,
    "\n"
  )
  out <- list(statistic = TStat, p.value = pvalue)
  return(out)
}

.seq2mat <- function(sequence) {
  n <- length(sequence)
  states <- unique(sequence)
  nstates <- length(states)
  
  mat <- matlab::zeros(nstates)
  dimnames(mat) <- list(states, states)
  
  for(i in 1:(n-1)) {
    from <- sequence[i]
    to <- sequence[i+1]
    mat[from, to] <- mat[from, to] + 1
  }
  for (i in 1:nstates){
    mat[i,] = mat[i,]/rowSums(mat)[i]
  }
  return (mat)
}


# phi function for divergence test
.phi <- function(x) {
  out <- x*log(x) - x + 1
  return(out)
}

# another phi function for divergence test
.phi2 <- function(x) {
  out <- 1/x
  return(out)
}
