#' Build a birth-death Markov chain
#'
#' Constructs a \code{markovchain} object for a birth-death process: a chain
#' on linearly ordered states \eqn{1,2,\ldots,n} that, from any state, can
#' only move to itself or to an immediately adjacent state.
#'
#' @param p A numeric vector of length \eqn{n-1}: \code{p[i]} is the "birth"
#'   probability of moving from state \eqn{i} up to state \eqn{i+1}.
#' @param q A numeric vector of length \eqn{n-1}, the same length as
#'   \code{p}: \code{q[i]} is the "death" probability of moving from state
#'   \eqn{i+1} down to state \eqn{i}.
#' @param states An optional character vector of \eqn{n=\code{length(p)}+1}
#'   state names. Defaults to \code{as.character(1:n)}.
#'
#' @return A new, row-stochastic \code{markovchain} object on \eqn{n}
#'   states, with transition matrix
#'   \deqn{P_{ii}=1-p_i-q_{i-1},\quad P_{i,i+1}=p_i,\quad P_{i,i-1}=q_{i-1}}
#'   (boundary terms \eqn{q_0} and \eqn{p_n} are understood to not exist,
#'   i.e. \eqn{P_{11}=1-p_1} and \eqn{P_{nn}=1-q_{n-1}}).
#'
#' @details
#' \strong{Why \code{p} and \code{q} have length \eqn{n-1}, not \eqn{n}.}
#' Every birth-death transition is a move between two adjacent states, and
#' there are exactly \eqn{n-1} adjacent pairs among \eqn{n} linearly ordered
#' states: \code{p[i]}/\code{q[i]} unambiguously describe the pair
#' \eqn{(i,i+1)}. This sidesteps a common source of confusion in this
#' construction, namely what to do with a stray "birth probability of the
#' top state" or "death probability of the bottom state" -- quantities
#' that do not correspond to any actual transition, since there is no state
#' \eqn{n+1} to be born into or state \eqn{0} to die into. Some
#' implementations accept two length-\eqn{n} vectors and quietly renormalize
#' every row so that any such leftover probability mass is redistributed
#' among the transitions that do exist; \code{birthDeath()} instead makes
#' the \eqn{n-1} genuine transition probabilities the only inputs, so there
#' is no leftover mass to (silently) dispose of in the first place.
#'
#' Every row's diagonal entry is determined by the requirement that the row
#' sums to \eqn{1}, so \code{p} and \code{q} alone fully determine \eqn{P}:
#' no separate "staying" probability is accepted or needed. \code{p+q} is
#' allowed to reach \eqn{1} for an interior state (no staying probability
#' there), but each element of \code{p} and \code{q} must itself lie in
#' \eqn{[0,1]} and \code{p[i]+q[i]} for the shared index \eqn{i} need not be
#' checked against 1 the way it would for a single state's own two
#' probabilities, since \code{p[i]} leaves state \eqn{i} while \code{q[i]}
#' leaves state \eqn{i+1}: the actual per-state constraint,
#' \eqn{p_i+q_{i-1}\le 1}, is checked directly on the assembled diagonal.
#'
#' The two boundary states \eqn{1} and \eqn{n} are reflecting only in the
#' weak sense that no birth/death carries them outside \eqn{\{1,\ldots,n\}}
#' -- they still generally have a positive probability of staying put
#' (\eqn{1-p_1} and \eqn{1-q_{n-1}} respectively) rather than being forced
#' to bounce back, unlike \code{\link{gamblersRuin}}'s absorbing ends or
#' \code{\link{toBoundedChain}}'s explicit reflecting condition, which can
#' be applied afterwards to force deterministic bouncing or absorption at
#' the ends of any chain, including one built here.
#'
#' @seealso \code{\link{gamblersRuin}}, \code{\link{toBoundedChain}},
#'   \code{\link{urnModel}}
#'
#' @examples
#' # A simple 4-state birth-death chain with constant birth/death rates.
#' bd <- birthDeath(p = c(0.3, 0.4, 0.5), q = c(0.2, 0.3, 0.1))
#' bd
#' rowSums(bd@transitionMatrix)
#'
#' @export
birthDeath <- function(p, q, states = NULL) {
  if (!is.numeric(p) || length(p) < 1L || anyNA(p) ||
      any(!is.finite(p)) || any(p < 0) || any(p > 1)) {
    stop("p must be a non-empty numeric vector with entries in [0, 1].")
  }
  if (!is.numeric(q) || length(q) != length(p) || anyNA(q) ||
      any(!is.finite(q)) || any(q < 0) || any(q > 1)) {
    stop("q must be a numeric vector of entries in [0, 1], the same length as p.")
  }

  m <- length(p) # number of adjacent pairs
  n <- m + 1L    # number of states

  if (!is.null(states)) {
    if (!is.character(states) || length(states) != n || anyNA(states) ||
        anyDuplicated(states)) {
      stop("states, if supplied, must be a character vector of ", n,
           " unique, non-missing names.")
    }
  } else {
    states <- as.character(seq_len(n))
  }

  P <- matrix(0, n, n, dimnames = list(states, states))
  # Interior diagonal: state i (2 <= i <= n-1) has birth p[i] and death
  # q[i-1] (the pair (i-1, i) is described by index i-1 of p/q).
  diagVals <- numeric(n)
  diagVals[1] <- 1 - p[1]
  if (n > 2L) {
    diagVals[2:(n - 1)] <- 1 - p[2:m] - q[1:(m - 1)]
  }
  diagVals[n] <- 1 - q[m]
  if (any(diagVals < -sqrt(.Machine$double.eps))) {
    stop(paste0(
      "Some implied diagonal entry is negative: for interior state i, ",
      "p[i] + q[i-1] must not exceed 1."
    ))
  }
  diag(P) <- pmax(diagVals, 0)

  for (i in seq_len(m)) {
    P[i, i + 1L] <- p[i]
    P[i + 1L, i] <- q[i]
  }

  new("markovchain", states = states, byrow = TRUE, transitionMatrix = P,
      name = "Birth-Death Chain")
}

#' Build a gambler's ruin Markov chain
#'
#' Constructs the classic gambler's ruin chain: a gambler with a fortune
#' between \code{0} and \code{upperBound} wins each round (and gains one
#' unit) with probability \code{prob}, otherwise loses one unit; play stops
#' as soon as the fortune reaches \code{0} (ruin) or \code{upperBound}
#' (the gambler's target).
#'
#' @param upperBound A single positive integer: the fortune at which the
#'   gambler stops (having won). The chain has \code{upperBound + 1} states,
#'   \eqn{0,1,\ldots,\code{upperBound}}.
#' @param prob A single number in \eqn{[0,1]}: the probability of winning
#'   an individual round (moving up by one unit) while the fortune is
#'   strictly between \code{0} and \code{upperBound}.
#' @param states An optional character vector of \code{upperBound + 1}
#'   state names, in increasing order of fortune. Defaults to
#'   \code{as.character(0:upperBound)}.
#'
#' @return A new, row-stochastic \code{markovchain} object with
#'   \code{upperBound + 1} states. States \code{"0"} and
#'   \code{as.character(upperBound)} are absorbing; every interior state
#'   \eqn{i} has \eqn{P_{i,i+1}=\code{prob}} and \eqn{P_{i,i-1}=1-\code{prob}}.
#'
#' @details
#' This is the special case of \code{\link{birthDeath}} with constant birth
#' probability \code{prob} and constant death probability \code{1-prob} at
#' every interior state, together forced to be \emph{absorbing} rather than
#' merely reflecting at the two ends -- which is why it is provided as its
#' own constructor rather than expressed purely in terms of
#' \code{birthDeath()}, which cannot produce absorbing boundaries by itself
#' (see \code{\link{toBoundedChain}} for turning any chain's ends
#' absorbing or reflecting after construction).
#'
#' With \code{prob != 0.5}, the classical ruin probability of reaching
#' \code{0} before \code{upperBound}, starting from fortune \eqn{i}, is
#' \deqn{P(\text{ruin}\mid X_0=i) =
#'   \frac{\left(\frac{1-\code{prob}}{\code{prob}}\right)^{i} -
#'         \left(\frac{1-\code{prob}}{\code{prob}}\right)^{\code{upperBound}}}
#'        {1-\left(\frac{1-\code{prob}}{\code{prob}}\right)^{\code{upperBound}}},}
#' and \eqn{i/\code{upperBound}} when \code{prob = 0.5}; this is a standard
#' textbook result (see Norris (1998), Section 1.3) and is not itself
#' computed by this function, but can be read off from
#' \code{\link{absorptionProbabilities}} applied to the returned chain.
#'
#' @references
#' Norris, J. R. (1998). \emph{Markov Chains}. Cambridge University Press.
#'
#' @seealso \code{\link{birthDeath}}, \code{\link{absorptionProbabilities}},
#'   \code{\link{toBoundedChain}}
#'
#' @examples
#' ruin <- gamblersRuin(upperBound = 5, prob = 0.4)
#' ruin
#' absorbingStates(ruin)
#'
#' @export
gamblersRuin <- function(upperBound, prob, states = NULL) {
  if (length(upperBound) != 1L || !is.numeric(upperBound) ||
      is.na(upperBound) || upperBound != as.integer(upperBound) ||
      upperBound < 2L) {
    stop("upperBound must be a single integer of at least 2.")
  }
  upperBound <- as.integer(upperBound)
  if (length(prob) != 1L || !is.numeric(prob) || is.na(prob) ||
      prob < 0 || prob > 1) {
    stop("prob must be a single number in [0, 1].")
  }

  n <- upperBound + 1L
  if (!is.null(states)) {
    if (!is.character(states) || length(states) != n || anyNA(states) ||
        anyDuplicated(states)) {
      stop("states, if supplied, must be a character vector of ", n,
           " unique, non-missing names.")
    }
  } else {
    states <- as.character(0:upperBound)
  }

  P <- matrix(0, n, n, dimnames = list(states, states))
  P[1, 1] <- 1
  P[n, n] <- 1
  if (n > 2L) {
    for (i in 2:(n - 1L)) {
      P[i, i - 1L] <- 1 - prob
      P[i, i + 1L] <- prob
    }
  }

  new("markovchain", states = states, byrow = TRUE, transitionMatrix = P,
      name = paste0("Gambler's Ruin (upperBound = ", upperBound, ")"))
}

#' Build an Ehrenfest urn model Markov chain
#'
#' Constructs the Ehrenfest diffusion model: \code{balls} balls are split
#' between two urns, A and B; at each step, one of the \code{balls} balls is
#' chosen uniformly at random and moved to the other urn. The chain tracks
#' the number of balls in urn A.
#'
#' @param balls A single positive integer, the total number of balls. The
#'   chain has \code{balls + 1} states, \eqn{0,1,\ldots,\code{balls}}
#'   (the possible counts of balls in urn A).
#' @param states An optional character vector of \code{balls + 1} state
#'   names, in increasing order of ball count. Defaults to
#'   \code{as.character(0:balls)}.
#'
#' @return A new, row-stochastic \code{markovchain} object with
#'   \code{balls + 1} states. From state \eqn{i} (\eqn{0<i<\code{balls}}),
#'   \deqn{P_{i,i-1} = i/\code{balls}, \qquad P_{i,i+1} = 1 - i/\code{balls},}
#'   the probability that the ball moved was one of the \eqn{i} currently in
#'   urn A (decreasing A's count) versus one of the \eqn{\code{balls}-i}
#'   currently in urn B (increasing it). States \eqn{0} and \code{balls}
#'   (all balls in one urn) are \emph{reflecting}: the next ball moved must
#'   come from the only non-empty urn, so \eqn{P_{0,1}=P_{\code{balls},
#'   \code{balls}-1}=1} exactly.
#'
#' @details
#' The Ehrenfest model is the classical example of a chain whose
#' equilibrium behaviour matches thermodynamic intuition despite every
#' individual transition being fully reversible: its stationary
#' distribution is \eqn{\mathrm{Binomial}(\code{balls}, 1/2)} (each ball is,
#' at equilibrium, independently in urn A or B with probability \eqn{1/2}),
#' sharply concentrated around \eqn{\code{balls}/2} for large
#' \code{balls} even though the chain only ever moves one ball at a time
#' and is reflecting, not absorbing, at the boundaries. It is irreducible
#' and reversible for every \code{balls}, but periodic with period \eqn{2}
#' (the parity of the ball count in urn A alternates every step): pass the
#' result through \code{\link{lazyChain}} first if an aperiodic chain is
#' needed, e.g. for \code{\link{mixingTime}}.
#'
#' @references
#' Ehrenfest, P. and Ehrenfest, T. (1907). Uber zwei bekannte Einwande gegen
#' das Boltzmannsche H-Theorem. \emph{Physikalische Zeitschrift}, 8,
#' 311-314.
#'
#' @seealso \code{\link{birthDeath}}, \code{\link{lazyChain}}
#'
#' @examples
#' ehrenfest <- urnModel(balls = 4)
#' ehrenfest
#' steadyStates(ehrenfest) # approximately Binomial(4, 0.5): 1/16 6/16 ...
#' dbinom(0:4, 4, 0.5)
#'
#' @export
urnModel <- function(balls, states = NULL) {
  if (length(balls) != 1L || !is.numeric(balls) || is.na(balls) ||
      balls != as.integer(balls) || balls < 1L) {
    stop("balls must be a single positive integer.")
  }
  balls <- as.integer(balls)

  n <- balls + 1L
  if (!is.null(states)) {
    if (!is.character(states) || length(states) != n || anyNA(states) ||
        anyDuplicated(states)) {
      stop("states, if supplied, must be a character vector of ", n,
           " unique, non-missing names.")
    }
  } else {
    states <- as.character(0:balls)
  }

  P <- matrix(0, n, n, dimnames = list(states, states))
  P[1, 2] <- 1
  P[n, n - 1L] <- 1
  if (n > 2L) {
    for (i in 2:(n - 1L)) {
      k <- i - 1L # ball count in urn A at this state (0-indexed)
      P[i, i - 1L] <- k / balls
      P[i, i + 1L] <- 1 - (k / balls)
    }
  }

  new("markovchain", states = states, byrow = TRUE, transitionMatrix = P,
      name = paste0("Ehrenfest Urn Model (balls = ", balls, ")"))
}

# Shared argument checking for tauchen() and rouwenhorst(): both discretize
# the same AR(1) process y_t = (1-rho)*alpha + rho*y_{t-1} + eps_t.
.checkAR1Args <- function(alpha, sigma, rho, size) {
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha)) {
    stop("alpha must be a single finite number.")
  }
  if (length(sigma) != 1L || !is.numeric(sigma) || !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be a single finite positive number.")
  }
  if (length(rho) != 1L || !is.numeric(rho) || !is.finite(rho) ||
      rho <= -1 || rho >= 1) {
    stop("rho must be a single finite number in (-1, 1).")
  }
  if (length(size) != 1L || !is.numeric(size) || is.na(size) ||
      size != as.integer(size) || size < 2L) {
    stop("size must be a single integer of at least 2.")
  }
  invisible(NULL)
}

#' Discretize an AR(1) process into a Markov chain (Tauchen's method)
#'
#' Approximates the stationary first-order autoregressive process
#' \deqn{y_t = (1-\rho)\alpha + \rho y_{t-1} + \varepsilon_t, \qquad
#'   \varepsilon_t \overset{\mathrm{iid}}{\sim} \mathcal N(0,\sigma^2)}
#' by a finite-state Markov chain on an evenly spaced grid, following
#' Tauchen (1986).
#'
#' @param alpha A single finite number: the unconditional mean of the
#'   process.
#' @param sigma A single finite positive number: the standard deviation of
#'   the innovation \eqn{\varepsilon_t}.
#' @param rho A single number in \eqn{(-1,1)}: the autocorrelation
#'   (persistence) of the process.
#' @param size A single integer of at least \code{2}: the number of grid
#'   points (states) of the discretized chain.
#' @param k A single positive number, the half-width of the grid in units
#'   of the process's unconditional standard deviation
#'   \eqn{\sigma_y=\sigma/\sqrt{1-\rho^2}}. The default, \code{3}, follows
#'   Tauchen (1986)'s own recommendation and covers the great majority of
#'   the stationary distribution's mass for typical \code{rho}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{chain}}{The discretized \code{markovchain} object, with
#'       state names equal to the grid values of \eqn{y} formatted to 4
#'       significant digits.}
#'     \item{\code{states}}{The numeric grid of \eqn{y}-values themselves,
#'       in the same order as \code{chain}'s states. Returning the actual
#'       levels alongside the chain, rather than only generic state labels
#'       \code{"1"}, \code{"2"}, \ldots, is deliberate: the whole point of
#'       discretizing an AR(1) process is usually to do further numeric
#'       work with the levels (e.g. plugging them into a pricing formula),
#'       and re-deriving the grid from \code{alpha}, \code{sigma}, \code{rho}
#'       and \code{k} a second time by hand is both extra work and a place
#'       for an off-by-one or rounding mismatch to creep in.}
#'   }
#'
#' @details
#' The grid is \eqn{n=\code{size}} evenly spaced points
#' \eqn{y_1<\cdots<y_n} spanning
#' \eqn{[\alpha-k\sigma_y,\ \alpha+k\sigma_y]}, with half-spacing
#' \eqn{w=(y_n-y_1)/(2(n-1))}. Writing \eqn{\Phi} for the standard normal
#' CDF, the transition probabilities from grid point \eqn{y_i} are
#' \deqn{P_{i1} = \Phi\!\left(\frac{y_1-(1-\rho)\alpha-\rho y_i+w}{\sigma}\right),}
#' \deqn{P_{in} = 1-\Phi\!\left(\frac{y_n-(1-\rho)\alpha-\rho y_i-w}{\sigma}\right),}
#' \deqn{P_{ij} = \Phi\!\left(\frac{y_j-(1-\rho)\alpha-\rho y_i+w}{\sigma}\right) -
#'   \Phi\!\left(\frac{y_j-(1-\rho)\alpha-\rho y_i-w}{\sigma}\right),
#'   \quad 1<j<n,}
#' i.e. the probability that \eqn{y_t} (a normal draw centred at the AR(1)
#' conditional mean) lands in the half-open bin around \eqn{y_j}, with the
#' two end bins extended to \eqn{\pm\infty} so that rows sum to exactly
#' \eqn{1}.
#'
#' Tauchen's method is simple and fast (\eqn{O(n^2)} normal CDF
#' evaluations) but the grid width is fixed by \code{k} regardless of
#' \code{size}: for a coarse grid (small \code{size}) it under-resolves the
#' bulk of the distribution, and for \code{rho} close to \eqn{\pm1} the true
#' unconditional variance is large and sensitive to \code{k}. See
#' \code{\link{rouwenhorst}} for an alternative that tends to match the
#' persistence of near-unit-root processes more accurately and needs no
#' arbitrary grid-width parameter.
#'
#' @references
#' Tauchen, G. (1986). Finite state markov-chain approximations to
#' univariate and vector autoregressions. \emph{Economics Letters}, 20(2),
#' 177-181.
#'
#' @seealso \code{\link{rouwenhorst}}
#'
#' @examples
#' out <- tauchen(alpha = 0, sigma = 1, rho = 0.9, size = 5)
#' out$states
#' out$chain
#' # The chain's own stationary variance should be close to the AR(1)'s
#' # theoretical unconditional variance sigma^2 / (1 - rho^2).
#' pi <- as.numeric(steadyStates(out$chain))
#' sum(pi * (out$states - sum(pi * out$states))^2)
#' 1 / (1 - 0.9^2)
#'
#' @export
tauchen <- function(alpha, sigma, rho, size, k = 3) {
  .checkAR1Args(alpha, sigma, rho, size)
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k <= 0) {
    stop("k must be a single finite positive number.")
  }
  n <- as.integer(size)

  yStd <- sigma / sqrt(1 - rho^2)
  yMax <- alpha + k * yStd
  yMin <- alpha - k * yStd
  y <- seq(yMin, yMax, length.out = n)
  w <- 0.5 * (yMax - yMin) / (n - 1)

  P <- matrix(0, n, n)
  condMean <- (1 - rho) * alpha + rho * y # length n, one per row
  for (i in seq_len(n)) {
    P[i, 1] <- stats::pnorm((y[1] - condMean[i] + w) / sigma)
    P[i, n] <- 1 - stats::pnorm((y[n] - condMean[i] - w) / sigma)
    if (n > 2L) {
      for (j in 2:(n - 1L)) {
        P[i, j] <- stats::pnorm((y[j] - condMean[i] + w) / sigma) -
          stats::pnorm((y[j] - condMean[i] - w) / sigma)
      }
    }
  }

  states <- format(signif(y, 4))
  states <- trimws(states)
  if (anyDuplicated(states)) {
    states <- as.character(seq_len(n))
  }
  dimnames(P) <- list(states, states)

  list(
    chain = new("markovchain", states = states, byrow = TRUE,
                transitionMatrix = P,
                name = paste0("Tauchen AR(1) Approximation (rho = ", rho, ")")),
    states = stats::setNames(y, states)
  )
}

#' Discretize an AR(1) process into a Markov chain (Rouwenhorst's method)
#'
#' Approximates the stationary first-order autoregressive process
#' \deqn{y_t = (1-\rho)\alpha + \rho y_{t-1} + \varepsilon_t, \qquad
#'   \varepsilon_t \overset{\mathrm{iid}}{\sim} \mathcal N(0,\sigma^2)}
#' by a finite-state Markov chain, following Rouwenhorst (1995). Unlike
#' \code{\link{tauchen}}, no grid-width parameter is needed and the method
#' remains accurate for \code{rho} close to \eqn{\pm 1}.
#'
#' @param alpha A single finite number: the unconditional mean of the
#'   process.
#' @param sigma A single finite positive number: the standard deviation of
#'   the innovation \eqn{\varepsilon_t}.
#' @param rho A single number in \eqn{(-1,1)}: the autocorrelation
#'   (persistence) of the process.
#' @param size A single integer of at least \code{2}: the number of grid
#'   points (states) of the discretized chain.
#'
#' @return A named list with two elements, \code{chain} and \code{states},
#'   in the same form as returned by \code{\link{tauchen}}.
#'
#' @details
#' The grid is \eqn{n=\code{size}} evenly spaced points spanning
#' \eqn{[\alpha-\psi,\ \alpha+\psi]} with
#' \eqn{\psi=\sigma_y\sqrt{n-1}}, where
#' \eqn{\sigma_y=\sigma/\sqrt{1-\rho^2}} is the process's unconditional
#' standard deviation: this particular width (rather than a fixed multiple
#' of \eqn{\sigma_y} as in \code{\link{tauchen}}) is what the method needs
#' in order to match the AR(1)'s variance and first-order autocorrelation
#' exactly at every \code{size}, including for \code{rho} near \eqn{\pm1}.
#' The transition matrix is built recursively. Let
#' \eqn{\theta=(1+\rho)/2} and, for two states,
#' \deqn{\Theta_2 = \begin{pmatrix}\theta & 1-\theta\\ 1-\theta & \theta\end{pmatrix}.}
#' For \eqn{m} states (\eqn{2<m\le n}), form the \eqn{m\times m} matrix
#' \deqn{\Theta_m = \theta\begin{pmatrix}\Theta_{m-1} & 0\\ 0 & 0\end{pmatrix}
#'   + (1-\theta)\begin{pmatrix}0 & \Theta_{m-1}\\ 0 & 0\end{pmatrix}
#'   + (1-\theta)\begin{pmatrix}0 & 0\\ \Theta_{m-1} & 0\end{pmatrix}
#'   + \theta\begin{pmatrix}0 & 0\\ 0 & \Theta_{m-1}\end{pmatrix},}
#' then divide every interior row (all but the first and last) by \eqn{2}
#' to restore row-stochasticity, since those rows receive contributions
#' from two of the four corner blocks above.
#'
#' Rouwenhorst's method reproduces the AR(1)'s unconditional variance and
#' lag-1 autocorrelation \eqn{\rho} exactly, for every \code{size}
#' (Kopecky and Suen (2010) find it outperforms \code{\link{tauchen}} and
#' several other methods across the persistence range typically seen in
#' quarterly macroeconomic and actuarial time series, e.g. discretized
#' short-rate or inflation processes for reserving and ALM work).
#'
#' @references
#' Rouwenhorst, K. G. (1995). Asset pricing implications of equilibrium
#' business cycle models. In T. F. Cooley (Ed.), \emph{Frontiers of
#' Business Cycle Research}, 294-330. Princeton University Press.
#'
#' Kopecky, K. A. and Suen, R. M. H. (2010). Finite state Markov-chain
#' approximations to highly persistent processes. \emph{Review of Economic
#' Dynamics}, 13(3), 701-714.
#'
#' @seealso \code{\link{tauchen}}
#'
#' @examples
#' out <- rouwenhorst(alpha = 0, sigma = 1, rho = 0.9, size = 5)
#' out$states
#' pi <- as.numeric(steadyStates(out$chain))
#' sum(pi * (out$states - sum(pi * out$states))^2) # matches 1/(1-rho^2) closely
#' 1 / (1 - 0.9^2)
#'
#' @export
rouwenhorst <- function(alpha, sigma, rho, size) {
  .checkAR1Args(alpha, sigma, rho, size)
  n <- as.integer(size)
  theta <- (1 + rho) / 2

  buildTheta <- function(m) {
    if (m == 2L) {
      return(matrix(c(theta, 1 - theta, 1 - theta, theta), 2, 2, byrow = TRUE))
    }
    inner <- buildTheta(m - 1L)
    out <- matrix(0, m, m)
    out[1:(m - 1L), 1:(m - 1L)] <- out[1:(m - 1L), 1:(m - 1L)] + theta * inner
    out[1:(m - 1L), 2:m] <- out[1:(m - 1L), 2:m] + (1 - theta) * inner
    out[2:m, 1:(m - 1L)] <- out[2:m, 1:(m - 1L)] + (1 - theta) * inner
    out[2:m, 2:m] <- out[2:m, 2:m] + theta * inner
    if (m > 2L) {
      out[2:(m - 1L), ] <- out[2:(m - 1L), ] / 2
    }
    out
  }

  P <- buildTheta(n)

  yStd <- sigma / sqrt(1 - rho^2)
  psi <- yStd * sqrt(n - 1)
  y <- seq(alpha - psi, alpha + psi, length.out = n)

  states <- format(signif(y, 4))
  states <- trimws(states)
  if (anyDuplicated(states)) {
    states <- as.character(seq_len(n))
  }
  dimnames(P) <- list(states, states)

  list(
    chain = new("markovchain", states = states, byrow = TRUE,
                transitionMatrix = P,
                name = paste0("Rouwenhorst AR(1) Approximation (rho = ", rho, ")")),
    states = stats::setNames(y, states)
  )
}

#' Build a population-genetics Markov chain (Moran or Wright-Fisher)
#'
#' Constructs the Markov chain tracking the number of copies of a mutant
#' allele in a population of \code{n} individuals under mutation and
#' viability selection, for either the Moran or the Wright-Fisher model of
#' population genetics.
#'
#' @param model Either \code{"moran"} or \code{"wright-fisher"}, selecting
#'   which reproduction scheme generates the transition matrix.
#' @param n A single integer of at least \code{2}: the (constant)
#'   population size. The chain has \code{n + 1} states,
#'   \eqn{0,1,\ldots,n}, the possible mutant-allele counts.
#' @param s A single finite number greater than \code{-1}: the selection
#'   coefficient. The mutant allele's fitness relative to the wild-type is
#'   \eqn{1+s} (\code{s = 0} is neutral drift, \code{s > 0} favours the
#'   mutant, \code{-1 < s < 0} disfavours it).
#' @param u A single number in \eqn{[0,1]}: the backward mutation rate,
#'   mutant to wild-type. Defaults to \code{1e-9}, matching the effectively
#'   mutation-free chain conventionally used for the neutral/absorbing case.
#' @param v A single number in \eqn{[0,1]}: the forward mutation rate,
#'   wild-type to mutant. Defaults to \code{1e-9}.
#' @param states An optional character vector of \code{n + 1} state names.
#'   Defaults to \code{as.character(0:n)}.
#'
#' @return A new, row-stochastic \code{markovchain} object on \code{n + 1}
#'   states. States \code{"0"} and \code{as.character(n)} (loss and
#'   fixation of the mutant allele) are absorbing, following the standard
#'   textbook treatment of both models; every interior state has a
#'   transition row determined by \code{model}, \code{s}, \code{u} and
#'   \code{v} as detailed below.
#'
#' @details
#' \strong{Moran model.} At each step one individual is chosen to reproduce
#' (with probability proportional to its type's relative fitness, mutant
#' fitness \eqn{r=1+s} against wild-type fitness \eqn{1}) and one
#' individual, chosen uniformly at random, dies and is replaced by the
#' offspring, which mutates with probability \code{u} or \code{v}
#' depending on the parent's type. From interior state \eqn{i}
#' (\eqn{0<i<n}), writing \eqn{r_i=(1+s)i} and \eqn{m_i=n-i},
#' \deqn{P_{i,i-1}=\frac{i}{n}\cdot\frac{r_i v + m_i(1-u)}{r_i+m_i},\qquad
#'   P_{i,i+1}=\frac{m_i}{n}\cdot\frac{r_i(1-u) + m_i v}{r_i+m_i},}
#' with \eqn{P_{ii}} the remainder.
#'
#' \strong{Wright-Fisher model.} Generations are discrete and
#' non-overlapping: from interior state \eqn{i}, the mutant-allele
#' frequency \eqn{k=i/n} is first updated for mutation,
#' \eqn{p=k(1-u)+(1-k)v}, then for viability selection,
#' \eqn{p'=\min\bigl(p(1+s)/(1+ps),\,1\bigr)}, and the next generation's
#' mutant count is \eqn{\mathrm{Binomial}(n,p')}-distributed:
#' \eqn{P_{ij}=\binom{n}{j}p'^{\,j}(1-p')^{n-j}}.
#'
#' Both constructions follow \pkg{PyDTMC}'s
#' \code{population_genetics_model}, with one deliberate difference: this
#' function uses mutant relative fitness \eqn{1+s} in \emph{both} models
#' (PyDTMC's Moran implementation instead uses \eqn{1-s}, so that a
#' positive \code{s} there favours the wild-type rather than the mutant,
#' the opposite of its own Wright-Fisher convention and of the usual
#' textbook one). \code{s = 0} (neutral drift) and the mutation-rate terms
#' are unaffected by this choice and match PyDTMC exactly; away from
#' \code{s = 0} the two packages' Moran chains differ by construction, by
#' design, to keep \code{s}'s meaning consistent between \code{"moran"} and
#' \code{"wright-fisher"} within this package. As \code{u, v} shrink to
#' \code{0}, both reduce to the classical drift-only chains with absorbing
#' loss/fixation states, going back to Wright (1931) and Moran (1958); see
#' Ewens (2004) for a modern textbook treatment of both.
#'
#' @references
#' Wright, S. (1931). Evolution in Mendelian populations. \emph{Genetics},
#' 16(2), 97-159.
#'
#' Moran, P. A. P. (1958). Random processes in genetics.
#' \emph{Mathematical Proceedings of the Cambridge Philosophical Society},
#' 54(1), 60-71.
#'
#' Ewens, W. J. (2004). \emph{Mathematical Population Genetics I:
#' Theoretical Introduction} (2nd ed.). Springer.
#'
#' @seealso \code{\link{birthDeath}}, \code{\link{gamblersRuin}}
#'
#' @examples
#' # Neutral drift (s = 0): absorption ("fixation") probability of the
#' # mutant allele from state i equals i/n, the classical result.
#' neutral <- populationGeneticsModel(model = "moran", n = 8, s = 0)
#' ap <- absorptionProbabilities(neutral)
#' ap["4", "8"] # close to 4/8 = 0.5
#'
#' # Positive selection increases the fixation probability from any
#' # interior starting count.
#' favoured <- populationGeneticsModel(model = "wright-fisher", n = 8, s = 0.5)
#' absorptionProbabilities(favoured)["4", "8"]
#'
#' @export
populationGeneticsModel <- function(model = c("moran", "wright-fisher"),
                                     n, s = 0, u = 1e-9, v = 1e-9,
                                     states = NULL) {
  model <- match.arg(model)
  if (length(n) != 1L || !is.numeric(n) || is.na(n) ||
      n != as.integer(n) || n < 2L) {
    stop("n must be a single integer of at least 2.")
  }
  n <- as.integer(n)
  if (length(s) != 1L || !is.numeric(s) || !is.finite(s) || s <= -1) {
    stop("s must be a single finite number greater than -1.")
  }
  if (length(u) != 1L || !is.numeric(u) || is.na(u) || u < 0 || u > 1) {
    stop("u must be a single number in [0, 1].")
  }
  if (length(v) != 1L || !is.numeric(v) || is.na(v) || v < 0 || v > 1) {
    stop("v must be a single number in [0, 1].")
  }

  size <- n + 1L
  if (!is.null(states)) {
    if (!is.character(states) || length(states) != size || anyNA(states) ||
        anyDuplicated(states)) {
      stop("states, if supplied, must be a character vector of ", size,
           " unique, non-missing names.")
    }
  } else {
    states <- as.character(0:n)
  }

  P <- matrix(0, size, size, dimnames = list(states, states))
  P[1L, 1L] <- 1
  P[size, size] <- 1

  ui <- 1 - u
  vi <- 1 - v

  if (model == "moran") {
    r <- 1 + s
    for (i in seq_len(n - 1L)) {
      nmi <- n - i
      ri <- r * i
      pm1 <- (i / n) * ((ri * v + nmi * vi) / (ri + nmi))
      pp1 <- (nmi / n) * ((ri * ui + nmi * v) / (ri + nmi))
      P[i + 1L, i] <- pm1
      P[i + 1L, i + 1L] <- 1 - pm1 - pp1
      P[i + 1L, i + 2L] <- pp1
    }
  } else {
    for (i in seq_len(n - 1L)) {
      k <- i / n
      pm <- k * ui + (1 - k) * v
      ps <- min((pm * (1 + s)) / (1 + pm * s), 1)
      P[i + 1L, ] <- stats::dbinom(0:n, n, ps)
    }
  }

  modelLabel <- if (model == "moran") "Moran" else "Wright-Fisher"
  new("markovchain", states = states, byrow = TRUE, transitionMatrix = P,
      name = paste0(modelLabel, " Population Genetics Model (n = ", n, ")"))
}
