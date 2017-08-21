# function to simulate a ctmc
rctmc <- function(n, ctmc, initDist = numeric(), T = 0, include.T0 = TRUE, 
                  out.type = "list") {
  if (identical(initDist, numeric()))
    state <- sample(ctmc@states, 1) # sample state randomly
  else if (length(initDist) != dim(ctmc) | round(sum(initDist), 5) != 1)
    stop("Error! Provide a valid initial state probability distribution")
  else 
    state <- sample(ctmc@states, 1, prob = initDist) # if valid probability distribution,
  # sample accordingly
  
  # obtain transition probability matrix from the generator matrix
  trans <- generatorToTransitionMatrix(ctmc@generator)
  
  
  states <- c()
  time <- c()
  if (include.T0 == TRUE){
    states <- c(states, state)
    time <- c(time, 0)
  }
  
  t <- 0
  i <- 1
  while (i <= n){
    idx <- which(ctmc@states == state)
    t <- t + rexp(1, -ctmc@generator[idx, idx])
    state <- ctmc@states[sample(1:dim(ctmc), 1, prob = trans[idx, ])]
    
    if(T > 0 & t > T)
      break
    
    states <- c(states, state)
    time <- c(time, t)
    i <- i + 1
  }
  
  out <- list(states, time)
  if (out.type == "list")
    return(out)
  else if(out.type == "df"){
    df <- data.frame(matrix(unlist(out), nrow = length(states)))
    names(df) <- c("states", "time")
    return(df)
  }
  else
    stop("Not a valid output type")
}

#' @title Return the generator matrix for a corresponding transition matrix
#' 
#' @description Calculate the generator matrix for a 
#'              corresponding transition matrix
#' 
#' @param P transition matrix between time 0 and t
#' @param t time of observation
#' @param method "logarithm" returns the Matrix logarithm of the transition matrix
#'
#' @return A matrix that represent the generator of P
#' @export
#'
#' @examples
#' mymatr <- matrix(c(.4, .6, .1, .9), nrow = 2, byrow = TRUE)
#' Q <- transition2Generator(P = mymatr)
#' expm::expm(Q)
#'  
#' @seealso \code{\link{rctmc}}
transition2Generator<-function(P, t = 1,method = "logarithm") {
  if (method == "logarithm") {
    Q = logm(P)/t
  } #else 
  return(Q)
}


#' Returns a generator matrix corresponding to frequency matrix
#' 
#' @description The function provides interface to calculate generator matrix corresponding to 
#' a frequency matrix and time taken
#' 
#' @param P relative frequency matrix
#' @param t (default value = 1)
#' @param method one among "QO"(Quasi optimaisation), "WA"(weighted adjustment), "DA"(diagonal adjustment)
#' @param logmethod method for computation of matrx algorithm (by default : Eigen)
#' 
#' @return returns a generator matix with same dimnames
#' 
#' @references E. Kreinin and M. Sidelnikova: Regularization Algorithms for
#' Transition Matrices. Algo Research Quarterly 4(1):23-40, 2001
#' 
#' @examples
#' sample <- matrix(c(150,2,1,1,1,200,2,1,2,1,175,1,1,1,1,150),nrow = 4,byrow = TRUE)
#' sample_rel = rbind((sample/rowSums(sample))[1:dim(sample)[1]-1,],c(rep(0,dim(sample)[1]-1),1)) 
#' freq2Generator(sample_rel,1)
#' 
#' data(tm_abs)
#' tm_rel=rbind((tm_abs/rowSums(tm_abs))[1:7,],c(rep(0,7),1))
#' ## Derive quasi optimization generator matrix estimate
#' freq2Generator(tm_rel,1)
#' 
#' @export
#' 
freq2Generator <- function(P,t = 1,method = "QO",logmethod = "Eigen"){
  if(method == "QO"){
    out <- ctmcd::gmQO(P, t, logmethod)
  } else if(method == "WA") {
    out <- ctmcd::gmWA(P, t, logmethod)
  } else if(method == "DA") {
    out <- ctmcd::gmDA(P, t, logmethod)
  }
  return(out)
}

#' @title Returns expected hitting time from state i to state j
#' 
#' @description Returns expected hitting time from state i to state j
#' 
#' @usage ExpectedTime(C,i,j,useRCpp)
#' 
#' @param C A CTMC S4 object
#' @param i Initial state i
#' @param j Final state j
#' @param useRCpp logical whether to use Rcpp
#' 
#' @details According to the theorem, holding times for all states except j should be greater than 0.
#' 
#' @return A numerical value that returns expected hitting times from i to j
#' 
#' @references Markovchains, J. R. Norris, Cambridge University Press
#' 
#' @author Vandit Jain
#' 
#' @examples
#' states <- c("a","b","c","d")
#' byRow <- TRUE
#' gen <- matrix(data = c(-1, 1/2, 1/2, 0, 1/4, -1/2, 0, 1/4, 1/6, 0, -1/3, 1/6, 0, 0, 0, 0),
#' nrow = 4,byrow = byRow, dimnames = list(states,states))
#' ctmc <- new("ctmc",states = states, byrow = byRow, generator = gen, name = "testctmc")
#' ExpectedTime(ctmc,1,4,TRUE)
#' 
#' @export
ExpectedTime <- function(C,i,j,useRCpp = TRUE){
  # take generator from ctmc-class object
  Q <- C@generator
  
  # in case where generator is written column wise
  if(C@byrow==FALSE){
    Q <- t(Q)
  }
  NoofStates <- dim(C)
  
  Exceptj <- c(1:NoofStates)
  # create vector with all values from 1:NoofStates except j 
  Exceptj <- which(Exceptj!=j)
  
  # build matrix with vlaues from Q such that row!=j or column!=j
  Q_Exceptj <- Q[Exceptj,Exceptj]
  
  # check for positivity of holding times except for state j
  if(!all(diag(Q_Exceptj)!=0)){
    stop("Holding times for all states except j should be greater than 0")
  }
  
  # get b for solving the system of linear equation Ax = b where A is Q_Exceptj
  b <- rep(-1,dim(Q_Exceptj)[1])
  
  # use solve function from base packge to solve Ax = b
  if(useRCpp == TRUE){
    out <- .ExpectedTimeRCpp(Q_Exceptj,b)
  } else {
    out <- solve(Q_Exceptj,b)
  }
  
  # out will be of size NoofStates-1, hence the adjustment for different cases of i>=<j
  if(i<j){
    return(out[[i]])
  } else if(i == j) {
    return(0)
  } else {
    return(out[[i-1]])
  }
}

#' Calculating probability from a ctmc object
#' 
#' @description 
#' This function returns the probability of every state at time t under different conditions
#' 
#' @usage probabilityatT(C,t,x0,useRCpp)
#' 
#' @param C A CTMC S4 object
#' @param t final time t
#' @param x0 initial state
#' @param useRCpp logical whether to use RCpp implementation
#' 
#' @details The initial state is not mandatory, In case it is not provided, 
#' function returns a matrix of transition function at time \code{t} else it returns
#' vector of probaabilities of transition to different states if initial state was \code{x0}
#' 
#' @return returns a vector or a matrix in case \code{x0} is provided or not respectively.
#' 
#' @references INTRODUCTION TO STOCHASTIC PROCESSES WITH R, ROBERT P. DOBROW, Wiley
#' 
#' @author Vandit Jain
#' 
#' @examples
#' states <- c("a","b","c","d")
#' byRow <- TRUE
#' gen <- matrix(data = c(-1, 1/2, 1/2, 0, 1/4, -1/2, 0, 1/4, 1/6, 0, -1/3, 1/6, 0, 0, 0, 0),
#' nrow = 4,byrow = byRow, dimnames = list(states,states))
#' ctmc <- new("ctmc",states = states, byrow = byRow, generator = gen, name = "testctmc")
#' probabilityatT(ctmc,1,useRCpp = TRUE)
#' 
#' @export
probabilityatT <- function(C, t, x0, useRCpp = TRUE){
  
  if(class(C) != "ctmc"){
    stop("Provided object is not a ctmc object")
  }
  if(t < 0){
    stop("Time provided should be greater than equal to 0")
  }
  # take generator from ctmc-class object
  Q <- C@generator
  
  # in case where generator is written column wise
  if(C@byrow==FALSE){
    Q <- t(Q)
  }
  NoofStates <- dim(C)
  
  
  # calculate transition functoin at time t using Kolmogorov backward equation
  if(useRCpp == TRUE){
    P <- .probabilityatTRCpp(t*Q);
  }
  else{
  P <- expm(t*Q)
  }
  
  
  # return the row vector according to state at time t=0 if x0 is provided
  # if x0 not provided returns the whole matrix
  # hence returns probability of being at state j at time t if state at t =0 is x0
  # where j is from 1:noofStates
  if(missing(x0)){
    P <- matrix(P, nrow = NoofStates, dimnames = list(C@states,C@states))
    return(P)
  } else {
    if(x0 > NoofStates || x0 < 1){
      stop("Initial state provided is not correct")
    }
    return(P[x0,])
  }
}




#' Calculating full conditional probability using lower rate transition matrix
#' 
#' This function calculates full conditional probability at given 
#' time s using lower rate transition matrix
#' 
#' @usage impreciseProbabilityatT(C,i,t,s,error,useRCpp)
#' 
#' @param C a ictmc class object
#' @param i initial state at time t
#' @param t initial time t. Default value = 0
#' @param s final time
#' @param error error rate. Default value = 0.001
#' @param useRCpp logical whether to use RCpp implementation; by default TRUE
#' 
#' @references Imprecise Continuous-Time Markov Chains, Thomas Krak et al., 2016
#' 
#' @author Vandit Jain
#' 
#' @examples
#' states <- c("n","y")
#' Q <- matrix(c(-1,1,1,-1),nrow = 2,byrow = TRUE,dimnames = list(states,states))
#' range <- matrix(c(1/52,3/52,1/2,2),nrow = 2,byrow = 2)
#' name <- "testictmc"
#' ictmc <- new("ictmc",states = states,Q = Q,range = range,name = name)
#' impreciseProbabilityatT(ictmc,2,0,1,10^-3,TRUE)
#'
impreciseProbabilityatT <- function(C, i, t=0, s, error = 10^-3, useRCpp = TRUE){
  ##  input validity checking
  if(s <= t){
    stop("Please provide time points such that initial time is greater than or equal to end point")
  }
  
  if(!class(C) == 'ictmc'){
    stop("Please provide a valid ictmc-class object")
  }
  noOfstates <-length(C@states)
  
  if(i <= 0 || i > noOfstates){
    stop("Please provide a valid initial state")
  }
  ### validity checking ends
  if(useRCpp == TRUE) {
    Qgx <- .impreciseProbabilityatTRCpp(C,i,t,s,error)
  } else {
  ## extract values from ictmc object 
  Q <- C@Q
  range <- C@range
  
  ### calculate ||QI_i||
  
  #initialise Q norm value
  QNorm <- -1
  
  for(i in 1:noOfstates){
    sum <- 0
    for(j in 1:noOfstates){
      sum <- sum + abs(Q[i,j])
    }
    QNorm <- max(sum*range[i,2],QNorm)
  }
  
  ### calculate no of iterations
  # The 1 is for norm of I_s i.e. ||I_s|| which equals 1
  n <- max((s-t)*QNorm, (s-t)*(s-t)*QNorm*QNorm*1/(2*error))
  
  ### calculate delta
  delta <- (s-t)/n
  
  ### build I_i vector
  Ii <- rep(0, noOfstates)
  Ii[i] <- 1
  
  
  ### calculate value of lower operator _QI_i(x) for all x belongs to no ofStates
  values <- Q%*%Ii
  Qgx <- rep(0, noOfstates)
  for(i in 1:noOfstates){
    Qgx[i] <- min(values[i]*range[i,1], values[i]*range[i,2])
  }
  Qgx <- delta*Qgx
  Qgx <- Ii + Qgx
  for(iter in 1:n-1){
    temp <- Qgx
    values <- Q%*%Qgx
    for(i in 1:noOfstates){
      Qgx[i] <- min(values[i]*range[i,1], values[i]*range[i,2])
    }
    Qgx <- delta*Qgx
    Qgx <- temp + Qgx
  }
  }
  return(Qgx)
}



#' Check if CTMC is irreducible
#' 
#' @description 
#' This function verifies whether a CTMC object is irreducible
#' 
#' @usage is.CTMCirreducible(ctmc)
#' 
#' @param ctmc a ctmc-class object
#' 
#' @references 
#' Continuous-Time Markov Chains, Karl Sigman, Columbia University
#' 
#' @author Vandit Jain
#' 
#' @return a boolean value as described above.
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
#' is.CTMCirreducible(molecularCTMC)
#' 
#' @export
is.CTMCirreducible <- function(ctmc) {
  
  if(!class(ctmc) == 'ctmc') {
    stop("please provide a valid ctmc class object")
  }
  
  ## gets the embeded chain matrix
  embeddedChainMatrix <- generatorToTransitionMatrix(ctmc@generator)
  
  
  ## forms a markovchain object related to embedded transition matrix
  markovchainObject <- new("markovchain",states = ctmc@states,
                           transitionMatrix = embeddedChainMatrix)
  
  ## returns result using is.irreducible function on embedded chain transition matrix
  return(is.irreducible(markovchainObject))
  
}



#' checks if ctmc object is time reversible
#' 
#' @description 
#' The function returns checks if provided function is time reversible
#' 
#' @usage is.TimeReversible(ctmc)
#' 
#' @param ctmc a ctmc-class object
#' 
#' @return Returns a boolean value stating whether ctmc object is time reversible
#' 
#' @author Vandit Jain
#' 
#' @return a boolean value as described above
#' 
#' @references 
#' INTRODUCTION TO STOCHASTIC PROCESSES WITH R, ROBERT P. DOBROW, Wiley
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
#' is.TimeReversible(molecularCTMC)
#' 
#' @export
is.TimeReversible <- function(ctmc) {
  
  if(!class(ctmc) == "ctmc") {
    stop("please provide a valid ctmc-class object")
  }
  
  ## get steady state probabilities
  Pi <- steadyStates(ctmc)
  
  ## initialise boolean result
  check <- TRUE
  
  ## no of states
  m <- length(ctmc@states)
  
  ## checks for byrow
  if(ctmc@byrow == FALSE)
    gen <- t(ctmc@generator)
  else
    gen <- ctmc@generator
  
  ## iterates for every state
  for( i in 1:m)
  {
    for(j in 1:m)
    {
      if(Pi[i]*gen[i,j] != Pi[j]*gen[j,i]) {
        check <- FALSE
        break
      }
    }
  }
  
  return(check)
  
}



# `generator/nextki` <- function(k) {
#   if(k >= 0) return(-1-k)
#   return(-k)
# }
# 
# `generator/nextk` <- function(k, Kmin, Kmax) {
#   if(is.null(k)) {
#     k <- rep(0, length(Kmin))
#     return(list(ans = TRUE, k = k))
#   }
#   
#   if(length(Kmin) == 0) {
#     return(list(ans = FALSE, k = k))
#   }
#   
#   i <- 1
#   kl <- k
#   kl[i] <- `generator/nextki`(kl[i])
#   while (kl[i] > Kmax[i] || kl[i] < Kmin[i]) {
#     kl[i] <- 0
#     i <- i+1
#     if(i > length(kl)) {
#       k <- kl
#       return(list(ans = FALSE, k = k))
#     }
#     kl[i] <- `generator/nextki`(kl[i])
#   }
#   k <- kl
#   return(list(ans = TRUE, k = k))
# }
# 
# `generator/generator` <- function(Po, Di, odi) {
#   P <- Po
#   N <- nrow(P)
#   
#   if(Di > 22) return(NULL) # bad idea
#   options(digits = Di)
#   odigs <- odi
#   
#   rSum <- rowSums(P)
#   if(! all(abs(1-rSum) < 0.001)) {
#     stop("Sum of each rows of Po should be equal to 1")
#   }
#   
#   P <- P/rSum
#   d <- det(P)
#   
#   if(d <= 0) {
#     cat("Matrix has non-positive determinant")
#     return(NULL)
#   }
#   
#   diagP <- 1
#   for(i in 1:nrow(P)) diagP <- diagP * P[i, i]
#   
#   if(d >= diagP) {
#     cat("Determinant exceeds product of diagonal elements\n")
#     return(NULL)
#   }
#   
#   E <- eigen(P)[[1]]
#   B <- eigen(P)[[2]]
#   
#   print("Eigenvalues")
#   print(E)
#   
#   # risky 
#   if(length(unique(E)) != length(E)) {
#     warning("Matrix does not have distinct eigenvalues")
#   }
#   
#   L <- abs(log(d))
#   addigs <- 2 + round(log10(1/Matrix::rcond(B))) + round(L/log(10)) # problem
#   
#   if(options()$digits < odigs + addigs) {
#     if(odigs + addigs > 100) {
#       print("Eigenvector matrix is singular")
#       return(NULL)
#     }
#     
#     cat('Going to', odigs + addigs, "digits")
#     return(`generator/generator`(Po, odigs + addigs, odigs))
#   }
#   
#   Bi <- solve(B)
#   
#   posevs <- NULL
#   negevs <- NULL
#   bestj <- NULL
#   bestQ <- NULL
#   marks <- rep(TRUE, length(E))
#   
#   for(i in 1:length(E)) { 
#     if(marks[i] && !(Re(E[i]) > 0 && Im(E[i]) == 0)) { # invalid comparison of complex number
#       cj <- Conj(E[i])
#       best <- Inf
#       if(i+1 <= length(E)) {
#         for(j in (i+1):length(E)) {
#          if(marks[j]) {
#            score <- abs(cj-E[j])
#            if(score < best) {
#              best <- score
#              bestj <- j
#            }
#          }
#         }
#       }
#         
#       if(best > 10^(3-options()$digits)) {
#         cat("Unpaired non-positive eigenvalue", E[i])
#         return(NULL)
#       }
#       marks[bestj] <- FALSE
#       if(Im(E[i]) >= 0) {
#         posevs <- c(posevs, i)
#         negevs <- c(negevs, bestj)
#         if(Im(E[bestj]) == 0) {
#           E[bestj] <- complex(real = E[bestj], imaginary = 0)
#         }
#       } else {
#         posevs <- c(posevs, bestj)  
#         negevs <- c(negevs, i)
#         if(Im(E[i]) == 0) {
#           E[i] <- complex(real = E[i], imaginary = 0)
#         }
#       }
#     }
#   }
#   
#   npairs <- length(posevs)
#   # display conjugate pairs
#   
#   Kmax <- rep(0, npairs)
#   Kmin <- Kmax
#   
#   for(i in 1:npairs) {
#     a <- Arg(E[posevs[i]])
#     Kmax[i] <- trunc((L-a)/2*pi)
#     Kmin[i] <- trunc((-L-a)/2*pi)
#   }
#   
#   # display K-max
#   # display K-min
#   
#   best <- -0.001
#   DD <- diag(log(E))
#   DK <- matlab::zeros(N)
#   res <- list(); p <- 1
#   k <- NULL
#   while(TRUE) {
#     
#     dlist <- `generator/nextk`(k, Kmin, Kmax)
#     k <- dlist$k
#     
#     if(dlist$ans == FALSE) {break}
#     
#     # display value of k
#     for(i in 1:npairs) {
#       ke <- complex(real = 0, imaginary = 2*pi*k[i])
#       DK[posevs[i], posevs[i]] <- ke
#       DK[negevs[i], negevs[i]] <- -ke
#     }
#     
#     Q <- B %*% (DD + DK) %*% Bi
#     # Q <- fnormal(Re(Q), options()$digits, 5*(10^(-1-odigs))) # define fnormal of maple
#     qmin <- Q[1,2]
#     for(i in 1:N) {
#       for(j in 1:N) {
#         if(i != j) {
#           if(Q[i, j] < qmin) qmin <- Q[i, j]
#         }
#       }
#     }
#     
#     if(EnvAllGenerators == TRUE) {
#       if(qmin > -.001) {
#         cat("Possible generator with qmin =", qmin)
#         res[[p]] <- round(Q, odigs)
#         p <- p + 1  
#       } else {
#         cat("qmin =", qmin)  
#       }
#       
#     } else {
#       if(qmin >= 0) {
#         cat("Found a generator")
#         return(round(Q, odigs))
#       } else {
#         if(qmin > best) {
#           best <- qmin
#           bestQ <- Q
#         }
#         if(qmin > -.001) {
#           cat("Approximate generator with qmin = ", qmin)
#         } else {
#           cat("qmin =", qmin)
#         }
#       }
#     }
#   }
#   
#   if(EnvAllGenerators == TRUE) {
#     return(res)
#   }
#   
#   warning("No completely valid generator found")
#   
#   if(! is.null(bestQ)) {
#     return(round(bestQ, odigs))
#   } else return(NULL)
#   
# }
# 
# generator <- function(Po, digits = 10) {
#   odigs <- digits
#   options(digits = 15)
#   if(is.matrix(Po)) {
#     P <- Po
#   } else {
#     stop("Po must be matrix")
#   }
#   
#   if(nrow(P) != ncol(P)) {
#     print(P)
#     stop('Po must be square matrix')
#   }
#   
#   if(! all(P >= 0)) {
#     print(P)
#     stop('Po must be non negative square matrix')
#   }
#   
#   `generator/generator`(P, options()$digits, odigs)
# }
