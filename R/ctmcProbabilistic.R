#' Simulate a continuous time Markov chain
#'
#' @description Simulates the jumps of a continuous time Markov chain: in each
#'   state \eqn{i} the chain waits an exponential time with rate
#'   \eqn{-Q_{ii}}, then jumps to \eqn{j \ne i} with probability
#'   \eqn{-Q_{ij} / Q_{ii}} (the embedded jump chain, see
#'   [generatorToTransitionMatrix()]).
#'
#' @param n Number of jumps to simulate: a positive integer, or `Inf` to
#'   simulate up to time `T`.
#' @param ctmc A [`ctmc`][ctmc-class] object.
#' @param initDist Distribution of the initial state, a probability vector in
#'   the order of `states(ctmc)`. The default, an empty vector, draws the initial
#'   state uniformly at random.
#' @param T Time horizon: the jumps after time `T` are not returned. The
#'   default `0` means no horizon, which requires a finite `n`.
#' @param include.T0 Logical: include the initial state, at time 0, in the
#'   output?
#' @param out.type `"list"` (default) or `"df"`, the form of the output.
#'
#' @details The simulation stops after `n` jumps, at the first jump after
#'   `T` (when `T > 0`), or when an absorbing state (zero exit rate) is
#'   reached, whichever comes first. The generator is used in the orientation
#'   given by the `byrow` slot of `ctmc`.
#'
#' @return With `out.type = "list"`, an unnamed list of two vectors of equal
#'   length: the states visited (character) and the times at which they were
#'   entered (numeric). With `out.type = "df"`, a data frame with the same
#'   content in the columns `states` and `time`.
#'
#' @references Anderson, D. F. (2013). Introduction to Stochastic Processes
#'   with Applications in the Biosciences. University of Wisconsin at Madison.
#'
#' @author Sai Bhargav Yalamanchi
#' @family continuous time Markov chains
#' @examples
#' energyStates <- c("sigma", "sigma_star")
#' gen <- matrix(c(-3, 3, 1, -1), nrow = 2, byrow = TRUE,
#'               dimnames = list(energyStates, energyStates))
#' molecularCTMC <- new("ctmc", states = energyStates, byrow = TRUE,
#'                      generator = gen, name = "Molecular Transition Model")
#'
#' # all the jumps up to time 1
#' rctmc(n = Inf, ctmc = molecularCTMC, T = 1)
#' # five jumps from an initial state drawn from c(0.8, 0.2)
#' rctmc(n = 5, ctmc = molecularCTMC, initDist = c(0.8, 0.2),
#'       include.T0 = FALSE, out.type = "df")
#' @export
# rctmc(): simulates a continuous time chain. Algorithm: in the current state i the chain waits an
# exponential time with rate -Q[i,i] (the total rate of leaving i), then jumps to a state j drawn
# with probabilities given by the embedded jump chain. It stops after n jumps or at time T.
rctmc <- function(n, ctmc, initDist = numeric(), T = 0, include.T0 = TRUE, 
                  out.type = "list") {
  # Step 1: validate the arguments.
  if (!is(ctmc, "ctmc")) stop("ctmc must be a ctmc object")
  if (length(n) != 1L || is.na(n) || n <= 0 ||
      (!is.finite(n) && !is.infinite(n)) ||
      (is.finite(n) && (n != floor(n) || n > .Machine$integer.max))) {
    stop("n must be a positive integer or Inf")
  }
  if (length(T) != 1L || is.na(T) || !is.finite(T) || T < 0) {
    stop("T must be a finite non-negative number")
  }
  if (is.infinite(n) && T <= 0) {
    stop("a finite positive T is required when n is infinite")
  }
  if (!out.type %in% c("list", "df")) {
    stop("out.type must be 'list' or 'df'")
  }
  # Step 2: choose the initial state (random, or according to the distribution initDist).
  if (identical(initDist, numeric()))
    state <- sample(ctmc@states, 1) # sample state randomly
  else if (!is.numeric(initDist) || length(initDist) != dim(ctmc) ||
           anyNA(initDist) || any(!is.finite(initDist)) || any(initDist < 0) ||
           abs(sum(initDist) - 1) > 1e-8)
    stop("Error! Provide a valid initial state probability distribution")
  else 
    state <- sample(ctmc@states, 1, prob = initDist) # if valid probability distribution,
  # sample accordingly
  
  # Step 3: probabilities of the jump chain: P[i,j] = Q[i,j] / (-Q[i,i]), for j != i.
  # The generator is taken by rows: a generator stored by columns (byrow = FALSE)
  # is transposed first (before version 1.2.1 it was read by rows, which
  # simulated the wrong chain).
  generator <- if (ctmc@byrow) ctmc@generator else t(ctmc@generator)
  trans <- generatorToTransitionMatrix(generator)
  
  
  states <- c()
  time <- c()
  if (include.T0) {
    states <- c(states, state)
    time <- c(time, 0)
  }
  
  # Step 4: simulation loop. `t` = current time; a state with zero exit rate is absorbing.
  t <- 0
  i <- 1
  while (i <= n){
    idx <- which(ctmc@states == state)
    if (generator[idx, idx] == 0) {
      # absorbing state; stay here forever
      t <- Inf
    } else {
      # Holding time in the current state: exponential with rate -Q[idx, idx].
      t <- t + rexp(1, -generator[idx, idx])
    }
    
    # Stop if the time limit T has been exceeded (or the chain is stuck in an absorbing state).
    if((T > 0 & t > T) | (is.infinite(t)))
      break
    
    # Next state: drawn with the jump probabilities of the current state.
    state <- ctmc@states[sample(1:dim(ctmc), 1, prob = trans[idx, ])]
    states <- c(states, state)
    time <- c(time, t)
    i <- i + 1
  }
  
  # Step 5: return the visited states and the jump times, as a list or a data frame.
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


#' Generator matrix from a transition matrix
#'
#' @description Computes the generator \eqn{Q} of a continuous time Markov
#'   chain whose transition matrix over a period of length `t` is `P`, as the
#'   matrix logarithm \eqn{Q = \log(P) / t}, so that
#'   \eqn{\exp(t Q) = P}.
#'
#' @param P A transition matrix (rows summing to one) over a period of length
#'   `t`.
#' @param t Length of the period, a positive number.
#' @param method Method used: only `"logarithm"`, the principal matrix
#'   logarithm computed by [expm::logm()], is available.
#'
#' @details Not every transition matrix is the transition matrix of a
#'   continuous time chain (the embedding problem): the matrix returned is a
#'   valid generator only if its off-diagonal entries are non-negative, which
#'   is not checked. [freq2Generator()] gives regularized estimates that are
#'   always valid generators.
#'
#' @return A square matrix with the dimensions of `P`, the generator
#'   \eqn{\log(P) / t}.
#'
#' @family continuous time Markov chains
#' @examples
#' P <- matrix(c(0.4, 0.6, 0.1, 0.9), nrow = 2, byrow = TRUE)
#' Q <- transition2Generator(P)
#' Q
#' expm::expm(Q) # gives back P
#' @export
# Generator from a transition matrix over a period t: Q = log(P) / t (matrix logarithm).
# It only gives a valid generator if the matrix logarithm has non-negative off-diagonal entries.
transition2Generator<-function(P, t = 1,method = "logarithm") {
  if (method == "logarithm") {
    Q = logm(P)/t
  } #else 
  return(Q)
}


#' Generator matrix estimated from a matrix of transition frequencies
#'
#' @description Estimates the generator of a continuous time Markov chain from
#'   the empirical transition matrix observed over a period of length `t`,
#'   with the regularization algorithms of Kreinin and Sidelnikova (2001),
#'   which return a valid generator (non-negative off-diagonal rates, rows
#'   summing to zero) even when the matrix logarithm of `P` is not one. The
#'   computation is done by the suggested package \pkg{ctmcd}.
#'
#' @param P Matrix of relative transition frequencies over the period: a
#'   square matrix whose rows sum to one.
#' @param t Length of the observation period, a positive number (default 1).
#' @param method Regularization algorithm: `"QO"` (quasi-optimization, the
#'   default), `"WA"` (weighted adjustment) or `"DA"` (diagonal adjustment).
#' @param logmethod Method of the matrix logarithm used by \pkg{ctmcd}
#'   (default `"Eigen"`, see `ctmcd::gmQO()`).
#'
#' @return The estimated generator, a matrix with the dimensions and dimnames
#'   of `P`; `NULL`, with a warning, if \pkg{ctmcd} is not installed.
#'
#' @references Kreinin, A. and Sidelnikova, M. (2001). Regularization
#'   algorithms for transition matrices. Algo Research Quarterly 4(1), 23-40.
#'
#' @family continuous time Markov chains
#' @examples
#' counts <- matrix(c(150, 2, 1, 1,
#'                    1, 200, 2, 1,
#'                    2, 1, 175, 1,
#'                    1, 1, 1, 150), nrow = 4, byrow = TRUE)
#' # relative frequencies, with the last state made absorbing
#' rel <- rbind((counts / rowSums(counts))[1:3, ], c(0, 0, 0, 1))
#' if (requireNamespace("ctmcd", quietly = TRUE)) freq2Generator(rel, 1)
#'
#' # rating transitions (quasi-optimization estimate)
#' data(tm_abs)
#' tmRel <- rbind((tm_abs / rowSums(tm_abs))[1:7, ], c(rep(0, 7), 1))
#' if (requireNamespace("ctmcd", quietly = TRUE)) freq2Generator(tmRel, 1)
#' @export
#' 
# Estimates a generator from the frequencies P of a discrete time process over a period t,
# using the methods of the optional package ctmcd (QO, WA, DA).
freq2Generator <- function(P,t = 1,method = "QO",logmethod = "Eigen"){
  if(requireNamespace('ctmcd', quietly = TRUE)) {
  if(method == "QO"){
    out <- ctmcd::gmQO(P, t, logmethod)
  } else if(method == "WA") {
    out <- ctmcd::gmWA(P, t, logmethod)
  } else if(method == "DA") {
    out <- ctmcd::gmDA(P, t, logmethod)
  }
  } else {
    warning("package ctmcd is not installed")
    out = NULL
  }
  return(out)
}

#' Expected hitting time of a state in a continuous time Markov chain
#'
#' @description Expected time for a continuous time Markov chain started in
#'   state `i` to reach state `j` for the first time.
#'
#' @param C A [`ctmc`][ctmc-class] object.
#' @param i,j Indices (positions in `states(C)`) of the starting and target
#'   states.
#' @param useRCpp Logical: solve the linear system in C++ (default) or with
#'   [solve()]?
#'
#' @details The expected hitting times \eqn{m_k} of `j` from the states
#'   \eqn{k \ne j} solve \eqn{\sum_{l \ne j} Q_{kl} m_l = -1} (Norris 1997,
#'   Section 3.3). Every state other than `j` must have a positive exit
#'   rate, otherwise an error is given; the times are finite only if `j` is
#'   reachable from every state, otherwise the system is singular.
#'
#' @return A number, the expected hitting time of `j` from `i` (0 if
#'   `i == j`).
#'
#' @references Norris, J. R. (1997). Markov Chains. Cambridge University
#'   Press.
#'
#' @author Vandit Jain
#' @family continuous time Markov chains
#' @examples
#' states <- c("a", "b", "c", "d")
#' gen <- matrix(c(-1, 1/2, 1/2, 0,
#'                 1/4, -1/2, 0, 1/4,
#'                 1/6, 0, -1/3, 1/6,
#'                 0, 0, 0, 0), nrow = 4, byrow = TRUE,
#'               dimnames = list(states, states))
#' ctmc <- new("ctmc", states = states, byrow = TRUE, generator = gen,
#'             name = "testctmc")
#' ExpectedTime(ctmc, 1, 4) # from "a" to the absorbing state "d"
#' @export
# Expected time to reach state j starting from state i. Let Q' be the generator restricted to the states
# other than j; the vector m of expected hitting times solves the linear system  Q' m = -1.
ExpectedTime <- function(C,i,j,useRCpp = TRUE){
  # take generator from ctmc-class object
  Q <- C@generator
  
  # in case where generator is written column wise
  if (!C@byrow) {
    Q <- t(Q)
  }
  NoofStates <- dim(C)
  
  Exceptj <- c(1:NoofStates)
  # create vector with all values from 1:NoofStates except j 
  Exceptj <- which(Exceptj!=j)
  
  # build matrix with vlaues from Q such that row!=j or column!=j
  Q_Exceptj <- Q[Exceptj,Exceptj,drop = FALSE]
  
  # check for positivity of holding times except for state j
  if(!all(diag(Q_Exceptj)!=0)){
    stop("Holding times for all states except j should be greater than 0")
  }
  
  # Right-hand side b = (-1, ..., -1) and the system  Q' x = b.
  # get b for solving the system of linear equation Ax = b where A is Q_Exceptj
  b <- rep(-1,dim(Q_Exceptj)[1])
  
  # use solve function from base packge to solve Ax = b
  if (useRCpp) {
    out <- .ExpectedTimeRCpp(Q_Exceptj,b)
  } else {
    out <- solve(Q_Exceptj,b)
  }
  
  # The solution omits state j: adjust the index for states i > j.
  # out will be of size NoofStates-1, hence the adjustment for different cases of i>=<j
  if(i<j){
    return(out[[i]])
  } else if(i == j) {
    return(0)
  } else {
    return(out[[i-1]])
  }
}

#' Transition probabilities of a continuous time Markov chain at time t
#'
#' @description Computes the transition matrix \eqn{P(t) = \exp(t Q)} of a
#'   continuous time Markov chain with generator \eqn{Q}, or one of its rows.
#'
#' @param C A [`ctmc`][ctmc-class] object.
#' @param t Time, a non-negative number.
#' @param x0 Optional index (position in `states(C)`) of the initial state.
#' @param useRCpp Logical: compute the matrix exponential in C++ (default) or
#'   with [expm::expm()]?
#'
#' @return Without `x0`, the matrix \eqn{P(t)}, with the states as dimnames,
#'   whose entry \eqn{(i, j)} is the probability of being in state `j` at time
#'   `t` when starting from `i`. With `x0`, row `x0` of that matrix: the
#'   distribution of the state at time `t` when starting from `x0`.
#'
#' @references Dobrow, R. P. (2016). Introduction to Stochastic Processes with
#'   R. Wiley.
#'
#' @author Vandit Jain
#' @family continuous time Markov chains
#' @examples
#' states <- c("a", "b", "c", "d")
#' gen <- matrix(c(-1, 1/2, 1/2, 0,
#'                 1/4, -1/2, 0, 1/4,
#'                 1/6, 0, -1/3, 1/6,
#'                 0, 0, 0, 0), nrow = 4, byrow = TRUE,
#'               dimnames = list(states, states))
#' ctmc <- new("ctmc", states = states, byrow = TRUE, generator = gen,
#'             name = "testctmc")
#' probabilityatT(ctmc, 1)
#' probabilityatT(ctmc, 1, x0 = 1)
#' @export
# Transition probabilities of a continuous chain at time t: P(t) = exp(t Q) (matrix exponential).
probabilityatT <- function(C, t, x0, useRCpp = TRUE){
  
  if(!is(C,"ctmc")){
    stop("Provided object is not a ctmc object")
  }
  if(t < 0){
    stop("Time provided should be greater than equal to 0")
  }
  # take generator from ctmc-class object
  Q <- C@generator
  
  # in case where generator is written column wise
  if (!C@byrow) {
    Q <- t(Q)
  }
  NoofStates <- dim(C)
  
  
  # Matrix exponential, in C++ or with expm::expm().
  # calculate transition functoin at time t using Kolmogorov backward equation
  if (useRCpp) {
    P <- .probabilityatTRCpp(t*Q);
  }
  else{
  P <- expm(t*Q)
  }
  
  
  # With x0 only the row of the starting state is returned, otherwise the whole matrix.
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




#' Lower transition probability of an imprecise continuous time Markov chain
#'
#' @description For an imprecise continuous time Markov chain, whose
#'   transition rates are only known to lie in intervals, computes the lower
#'   probability of being in state `i` at time `s`, as a function of the state
#'   at time `t` (Krak et al. 2017).
#'
#' @param C An `ictmc` object: its slot `Q` is a generator and row \eqn{k} of
#'   the rates may be scaled by any factor in `range[k, ]`.
#' @param i Index (position in the states of `C`) of the state whose lower
#'   probability is computed.
#' @param t Initial time (default 0).
#' @param s Final time, larger than `t`.
#' @param error Maximum approximation error (default 0.001).
#' @param useRCpp Logical: use the C++ implementation (default)? It requires
#'   `t` and `s` to be integer-valued.
#'
#' @details The lower transition operator is approximated by
#'   \eqn{n} steps \eqn{g \leftarrow g + \delta \underline{Q} g}, with
#'   \eqn{\delta = (s - t) / n} and \eqn{n} chosen so that the error is below
#'   `error`, starting from the indicator of state `i`; \eqn{\underline{Q}}
#'   takes, in each row, the smaller of the values given by the two extremes of
#'   the interval of rates.
#'
#' @return A numeric vector with one element per state \eqn{x}: the lower
#'   probability that the chain is in state `i` at time `s` given that it is in
#'   \eqn{x} at time `t`. With degenerate intervals (`range` equal to 1) it is
#'   column `i` of \eqn{\exp((s - t) Q)}, up to `error`.
#'
#' @references Krak, T., De Bock, J. and Siebes, A. (2017). Imprecise
#'   continuous-time Markov chains. International Journal of Approximate
#'   Reasoning 88, 452-528.
#'
#' @author Vandit Jain
#' @family continuous time Markov chains
#' @examples
#' states <- c("n", "y")
#' Q <- matrix(c(-1, 1, 1, -1), nrow = 2, byrow = TRUE,
#'             dimnames = list(states, states))
#' range <- matrix(c(1/52, 3/52, 1/2, 2), nrow = 2, byrow = TRUE)
#' ictmc <- new("ictmc", states = states, Q = Q, range = range,
#'              name = "testictmc")
#' impreciseProbabilityatT(ictmc, i = 2, t = 0, s = 1, error = 1e-3)
#' @export
# Imprecise CTMC: the rates are only known to lie within intervals. This function computes the LOWER
# expectation of the indicator of state i at time s (Krak et al. 2017) by repeated small Euler steps:
#   g <- g + delta * lower(Q g),  delta = (s - t)/n ,
# where the number of steps n is chosen so that the error is below `error`.
impreciseProbabilityatT <- function(C, i, t=0, s, error = 10^-3, useRCpp = TRUE){
  if (length(error) != 1L || is.na(error) || !is.finite(error) || error <= 0) {
    stop("error must be a finite positive number")
  }
  if (length(t) != 1L || length(s) != 1L || anyNA(c(t, s)) ||
      any(!is.finite(c(t, s)))) {
    stop("t and s must be finite scalar time points")
  }
  ##  input validity checking
  if(s <= t){
    stop("Please provide time points such that initial time is greater than or equal to end point")
  }
  
  if(!is(C,'ictmc')){
    stop("Please provide a valid ictmc-class object")
  }
  noOfstates <-length(C@states)
  if (length(i) != 1L || is.na(i) || !is.finite(i) || i != floor(i)) {
    stop("Please provide a valid initial state")
  }
  if (useRCpp && (t != floor(t) || s != floor(s))) {
    stop("t and s must be integer-valued when useRCpp is TRUE")
  }
  
  if(i <= 0 || i > noOfstates){
    stop("Please provide a valid initial state")
  }
  ### validity checking ends
  if (useRCpp) {
    Qgx <- .impreciseProbabilityatTRCpp(C,i,t,s,error)
  } else {
  ## extract values from ictmc object 
  Q <- C@Q
  range <- C@range
  
  # Norm of the operator Q: the largest row sum of |Q| times the upper range.
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
  
  # n = number of Euler steps needed to reach the required precision.
  ### calculate no of iterations
  # The 1 is for norm of I_s i.e. ||I_s|| which equals 1
  n <- max((s-t)*QNorm, (s-t)*(s-t)*QNorm*QNorm*1/(2*error))
  
  ### calculate delta
  delta <- (s-t)/n
  
  ### build I_i vector
  Ii <- rep(0, noOfstates)
  Ii[i] <- 1
  
  
  # Lower operator: for each state take the minimum of the two extremes of the interval.
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



#' Check whether a continuous time Markov chain is irreducible
#'
#' @description A continuous time Markov chain is irreducible if every state
#'   can be reached from every other state, which is the case if and only if
#'   its embedded jump chain is irreducible.
#'
#' @param ctmc A [`ctmc`][ctmc-class] object.
#'
#' @return `TRUE` if the chain is irreducible, `FALSE` otherwise.
#'
#' @references Sigman, K. (2009). Continuous-time Markov chains. Lecture
#'   notes, Columbia University.
#'
#' @author Vandit Jain
#' @family continuous time Markov chains
#' @examples
#' energyStates <- c("sigma", "sigma_star")
#' gen <- matrix(c(-3, 3, 1, -1), nrow = 2, byrow = TRUE,
#'               dimnames = list(energyStates, energyStates))
#' molecularCTMC <- new("ctmc", states = energyStates, byrow = TRUE,
#'                      generator = gen, name = "Molecular Transition Model")
#' is.CTMCirreducible(molecularCTMC)
#' @export
# A CTMC is irreducible iff its embedded jump chain is irreducible.
is.CTMCirreducible <- function(ctmc) {
  
  if(!is(ctmc, 'ctmc') ) {
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



#' Check whether a continuous time Markov chain is time reversible
#'
#' @description A continuous time Markov chain with stationary distribution
#'   \eqn{\pi} is time reversible if it satisfies the detailed balance
#'   equations \eqn{\pi_i Q_{ij} = \pi_j Q_{ji}} for all states \eqn{i, j}.
#'
#' @param ctmc A [`ctmc`][ctmc-class] object.
#'
#' @details The equations are checked up to a relative tolerance of
#'   \eqn{\sqrt{\epsilon}} (about \eqn{1.5 \times 10^{-8}}) of the largest
#'   flow \eqn{\pi_i |Q_{ij}|}, because \eqn{\pi} is computed numerically.
#'   When the chain has several stationary distributions they must all
#'   satisfy them.
#'
#' @return `TRUE` if the chain is time reversible, `FALSE` otherwise.
#'
#' @references Dobrow, R. P. (2016). Introduction to Stochastic Processes with
#'   R. Wiley.
#'
#' @author Vandit Jain
#' @family continuous time Markov chains
#' @examples
#' energyStates <- c("sigma", "sigma_star")
#' gen <- matrix(c(-3, 3, 1, -1), nrow = 2, byrow = TRUE,
#'               dimnames = list(energyStates, energyStates))
#' molecularCTMC <- new("ctmc", states = energyStates, byrow = TRUE,
#'                      generator = gen, name = "Molecular Transition Model")
#' is.TimeReversible(molecularCTMC)
#'
#' # a birth-death chain is reversible, a cycle is not
#' s <- c("a", "b", "c")
#' birthDeath <- matrix(c(-0.3, 0.3, 0,
#'                        0.1, -0.5, 0.4,
#'                        0, 0.7, -0.7), nrow = 3, byrow = TRUE,
#'                      dimnames = list(s, s))
#' is.TimeReversible(new("ctmc", states = s, byrow = TRUE, generator = birthDeath))
#' cycle <- matrix(c(-1, 1, 0, 0, -1, 1, 1, 0, -1), nrow = 3, byrow = TRUE,
#'                 dimnames = list(s, s))
#' is.TimeReversible(new("ctmc", states = s, byrow = TRUE, generator = cycle))
#' @export
# Time reversibility (detailed balance in continuous time): pi_i Q_ij = pi_j Q_ji for all i, j.
is.TimeReversible <- function(ctmc) {
  
  if(!is(ctmc,"ctmc")) {
    stop("please provide a valid ctmc-class object")
  }
  
  ## get steady state probabilities
  Pi <- steadyStates(ctmc)
  
  ## initialise boolean result
  check <- TRUE
  
  ## checks for byrow
  if (!ctmc@byrow)
    gen <- t(ctmc@generator)
  else
    gen <- ctmc@generator
  
  # Detailed balance pi_i Q_ij = pi_j Q_ji, checked on the flow matrix
  # F = diag(pi) Q against its transpose, up to a relative tolerance: pi is
  # computed numerically, so an exact comparison (as before version 1.2.1)
  # wrongly rejected reversible chains. With several stationary
  # distributions (one per row, or per column if byrow = FALSE) each must
  # satisfy it.
  if (is.null(Pi)) return(FALSE)
  Pi <- as.matrix(Pi)
  if (!ctmc@byrow) Pi <- t(Pi)
  for (k in seq_len(nrow(Pi))) {
    flows <- Pi[k, ] * gen
    scale <- max(abs(flows))
    if (any(abs(flows - t(flows)) > sqrt(.Machine$double.eps) * scale)) check <- FALSE
  }
  
  return(check)
  
}



# Old, disabled code kept for reference (a search for generators of a transition matrix).
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
