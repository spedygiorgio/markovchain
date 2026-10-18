# Markov chains obtained by discretizing a continuous transition kernel with a
# quadrature rule (the Nystrom method), and the Tauchen-Hussey approximation
# of an AR(1) process, which is the same construction with Gauss-Hermite
# nodes.
#
# The idea, for a junior reader: a continuous-state process moves from x to a
# point y with density f(x, y). For any function h,
#   E[h(Y) | X = x] = integral of h(y) f(x, y) dy
#                   ~ sum_j w_j h(y_j) f(x, y_j)          (a quadrature rule)
# so, on the quadrature nodes y_1, ..., y_n, the numbers w_j f(x_i, y_j),
# divided by their row sum, form a transition matrix whose expectations
# approximate those of the continuous process.

# Internal helper: nodes and weights of a quadrature rule for the integral of a
# function with respect to the Lebesgue measure (plain dx), so that
#   integral g(x) dx  ~  sum_k weights[k] * g(nodes[k]).
# Every rule returned here has positive weights, which is what makes the
# resulting matrix a valid transition matrix.
.quadratureRule <- function(type, n, interval = c(0, 1), mean = 0, sd = 1) {
  # Gauss rules from the eigen-decomposition of the Jacobi (tridiagonal)
  # matrix of the orthogonal polynomials (Golub and Welsch 1969): the nodes
  # are the eigenvalues and the weights mu0 times the squared first
  # components of the normalized eigenvectors.
  golubWelsch <- function(offDiagonal, mu0) {
    J <- matrix(0, n, n)
    if (n > 1L) {
      k <- seq_len(n - 1L)
      J[cbind(k, k + 1L)] <- offDiagonal
      J[cbind(k + 1L, k)] <- offDiagonal
    }
    e <- eigen(J, symmetric = TRUE)
    o <- order(e$values)
    list(nodes = e$values[o], weights = mu0 * e$vectors[1L, o]^2)
  }
  a <- interval[1L]
  b <- interval[2L]
  # map a rule on [-1, 1] to [a, b]
  toInterval <- function(t, w) list(nodes = (b - a) / 2 * t + (a + b) / 2,
                                    weights = (b - a) / 2 * w)
  switch(type,
    "gauss-legendre" = {
      k <- seq_len(n - 1L)
      r <- golubWelsch(k / sqrt(4 * k^2 - 1), 2)
      # the eigenvalues of a symmetric matrix are only accurate to rounding:
      # symmetrize the nodes around the centre, as the exact rule is
      r$nodes <- (r$nodes - rev(r$nodes)) / 2
      r$weights <- (r$weights + rev(r$weights)) / 2
      toInterval(r$nodes, r$weights)
    },
    "gauss-chebyshev" = {
      # Gauss-Chebyshev (first kind) integrates g(t) / sqrt(1 - t^2) with
      # equal weights pi / n; for a plain integral the factor sqrt(1 - t^2)
      # moves into the weights.
      t <- -cos((2 * seq_len(n) - 1) * pi / (2 * n))
      toInterval(t, pi / n * sqrt(1 - t^2))
    },
    "fejer" = {
      # Fejer's first rule: the same Chebyshev nodes, with the weights that
      # integrate exactly every polynomial of degree n - 1 (Fejer 1933).
      theta <- (2 * seq_len(n) - 1) * pi / (2 * n)
      m <- seq_len(n %/% 2L)
      w <- vapply(theta, function(th)
        2 / n * (1 - 2 * sum(cos(2 * m * th) / (4 * m^2 - 1))), numeric(1))
      toInterval(-cos(theta), w)
    },
    "gauss-hermite" = {
      # Gauss-Hermite for the weight exp(-t^2), turned into a rule for the
      # normal distribution N(mean, sd^2): nodes mean + sqrt(2) sd t,
      # probability weights w / sqrt(pi). Dividing by the normal density at
      # each node gives weights for the plain integral.
      k <- seq_len(n - 1L)
      r <- golubWelsch(sqrt(k / 2), sqrt(pi))
      r$nodes <- (r$nodes - rev(r$nodes)) / 2
      r$weights <- (r$weights + rev(r$weights)) / 2
      z <- mean + sqrt(2) * sd * r$nodes
      list(nodes = z, weights = (r$weights / sqrt(pi)) / stats::dnorm(z, mean, sd))
    },
    "trapezoid" = {
      x <- seq(a, b, length.out = n)
      w <- rep((b - a) / (n - 1), n)
      w[c(1L, n)] <- w[c(1L, n)] / 2
      list(nodes = x, weights = w)
    },
    "simpson" = {
      x <- seq(a, b, length.out = n)
      # 1, 4, 2, 4, ..., 2, 4, 1 times h / 3
      w <- rep(c(2, 4), length.out = n)
      w[c(1L, n)] <- 1
      list(nodes = x, weights = (b - a) / (n - 1) / 3 * w)
    },
    "midpoint" = {
      h <- (b - a) / n
      list(nodes = a + (seq_len(n) - 0.5) * h, weights = rep(h, n))
    }
  )
}

# Internal helper: state names from numeric nodes, as tauchen() does.
.nodeStateNames <- function(x) {
  states <- trimws(format(signif(x, 4)))
  if (anyDuplicated(states)) {
    states <- as.character(seq_along(x))
  }
  states
}

#' Markov chain approximating a continuous transition kernel by quadrature
#'
#' Builds a finite-state Markov chain from a transition density (or any
#' non-negative transition kernel) of a continuous-state process, by
#' evaluating the kernel at the nodes of a quadrature rule (the Nystrom
#' method). It is the counterpart of PyDTMC's \code{fit_function}.
#'
#' @param f A function of two arguments, \code{f(x, y)}, returning the
#'   (unnormalized) density of moving from \code{x} to \code{y}. It should be
#'   vectorized: it is called once with two vectors of length
#'   \code{size^2}. A function that is not vectorized also works (it is then
#'   called once per pair), more slowly. Values must be finite and
#'   non-negative.
#' @param size A single integer of at least \code{2}: the number of nodes,
#'   i.e. of states of the chain.
#' @param quadrature The quadrature rule, one of \code{"gauss-legendre"}
#'   (the default), \code{"gauss-chebyshev"}, \code{"fejer"},
#'   \code{"gauss-hermite"}, \code{"trapezoid"}, \code{"simpson"} and
#'   \code{"midpoint"}; see Details.
#' @param interval A numeric vector \code{c(a, b)} with \code{a < b}, the
#'   interval on which the process lives, for every rule except
#'   \code{"gauss-hermite"}.
#' @param mean,sd Mean and standard deviation of the normal distribution
#'   that defines the nodes of \code{"gauss-hermite"}, which covers the
#'   whole real line; ignored by the other rules.
#' @param states Optional character vector of \code{size} distinct state
#'   names. By default the nodes, formatted to 4 significant digits, are
#'   used (or \code{"1"}, \ldots, \code{"size"} if that formatting is not
#'   unique).
#' @param name Optional name of the chain.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{chain}}{The \code{markovchain} object (row-stochastic).}
#'     \item{\code{states}}{The nodes, a numeric vector named by the states,
#'       in the order of the states of \code{chain}.}
#'     \item{\code{weights}}{The quadrature weights for the plain integral
#'       \eqn{\int g(y)\,dy \approx \sum_j w_j g(y_j)}, named by the states.}
#'   }
#'
#' @details
#' With nodes \eqn{y_1, \ldots, y_n} and weights \eqn{w_1, \ldots, w_n} of
#' a rule for \eqn{\int g(y)\,dy}, the transition matrix is
#' \deqn{P_{ij} = \frac{w_j f(y_i, y_j)}{\sum_{k} w_k f(y_i, y_k)}.}
#' If \eqn{f(x, \cdot)} is the transition density of the continuous
#' process, \eqn{\sum_j P_{ij} h(y_j)} approximates
#' \eqn{E[h(Y_{t+1}) \mid Y_t = y_i]} for every smooth \eqn{h}, and the
#' normalization corrects the mass lost outside the nodes. The kernel does
#' not need to integrate to one: only its values at the nodes matter.
#'
#' The rules, all with positive weights:
#' \itemize{
#'   \item \code{"gauss-legendre"}: Gauss-Legendre nodes, exact for
#'     polynomials of degree \eqn{2n - 1}. A good default for smooth kernels
#'     on a bounded interval.
#'   \item \code{"gauss-chebyshev"}: Gauss-Chebyshev nodes of the first
#'     kind, \eqn{t_k = -\cos((2k - 1)\pi / (2n))} on \eqn{[-1, 1]}, with the
#'     weights \eqn{(\pi/n)\sqrt{1 - t_k^2}} of a plain integral. It is exact
#'     for \eqn{g(t) = q(t)/\sqrt{1 - t^2}} with \eqn{q} a polynomial of
#'     degree \eqn{2n - 1}, so it suits kernels with that kind of
#'     behaviour at the end points.
#'   \item \code{"fejer"}: Fejer's first rule, the same Chebyshev nodes with
#'     the weights that integrate exactly every polynomial of degree
#'     \eqn{n - 1}. This is the rule that PyDTMC calls
#'     \code{"gauss-chebyshev"}.
#'   \item \code{"gauss-hermite"}: Gauss-Hermite nodes for the normal
#'     distribution \eqn{N(\code{mean}, \code{sd}^2)}, for processes on the
#'     real line. The weights of the plain integral are the probability
#'     weights divided by the normal density at the nodes. With an AR(1)
#'     density this is the method of Tauchen and Hussey (1991); see
#'     \code{\link{tauchenHussey}}.
#'   \item \code{"trapezoid"} and \code{"simpson"}: equally spaced nodes
#'     including the end points; Simpson's rule needs an odd \code{size}.
#'   \item \code{"midpoint"}: the centres of \code{size} equal cells.
#' }
#' Newton-Cotes rules of high order are not offered, because their weights
#' become negative and would give negative probabilities.
#'
#' A row whose kernel values are all zero at the nodes (the process cannot
#' reach any node from there) stops the function with an error: PyDTMC
#' replaces such a row with a uniform distribution, which invents
#' transitions; here the interval, the number of nodes or the kernel should
#' be changed instead.
#'
#' For the same kernel, interval and number of nodes, the
#' \code{"gauss-legendre"}, \code{"fejer"}, \code{"trapezoid"} and
#' \code{"simpson"} rules give the same matrices as PyDTMC 10.0.0
#' (\code{"gauss-legendre"}, \code{"gauss-chebyshev"},
#' \code{"trapezoid-rule"} and \code{"simpson-rule"}); PyDTMC's midpoint
#' rule (\code{"newton-cotes"}) ignores the lower end of the interval.
#'
#' @references
#' Tauchen, G. and Hussey, R. (1991). Quadrature-based methods for
#' obtaining approximate solutions to nonlinear asset pricing models.
#' \emph{Econometrica}, 59(2), 371-396. \doi{10.2307/2938261}
#'
#' Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss quadrature
#' rules. \emph{Mathematics of Computation}, 23(106), 221-230.
#' \doi{10.1090/S0025-5718-69-99647-1}
#'
#' Fejer, L. (1933). Mechanische Quadraturen mit positiven Cotesschen
#' Zahlen. \emph{Mathematische Zeitschrift}, 37, 287-309.
#' \doi{10.1007/BF01474575}
#'
#' @seealso \code{\link{tauchenHussey}}, \code{\link{tauchen}},
#'   \code{\link{rouwenhorst}}
#'
#' @examples
#' # A process on [0, 1] that moves towards 0.6 times its current value,
#' # with a Gaussian kernel truncated to the interval.
#' kernel <- function(x, y) dnorm(y, mean = 0.6 * x + 0.2, sd = 0.15)
#' out <- fitFunction(kernel, size = 7, interval = c(0, 1))
#' out$chain
#' # expected next value from each node, against the continuous one
#' cbind(chain = as.numeric(out$chain@transitionMatrix %*% out$states),
#'       target = 0.6 * out$states + 0.2)
#'
#' @export
fitFunction <- function(f, size, quadrature = c("gauss-legendre", "gauss-chebyshev",
                                                "fejer", "gauss-hermite", "trapezoid",
                                                "simpson", "midpoint"),
                        interval = c(0, 1), mean = 0, sd = 1, states = NULL, name = NULL) {
  # Step 1: check the arguments.
  if (!is.function(f)) {
    stop("f must be a function of two arguments, f(x, y).")
  }
  quadrature <- match.arg(quadrature)
  if (length(size) != 1L || !is.numeric(size) || is.na(size) ||
      size != as.integer(size) || size < 2L) {
    stop("size must be a single integer of at least 2.")
  }
  n <- as.integer(size)
  if (quadrature == "simpson" && n %% 2L == 0L) {
    stop("Simpson's rule needs an odd size.")
  }
  if (quadrature == "gauss-hermite") {
    if (length(mean) != 1L || !is.numeric(mean) || !is.finite(mean)) {
      stop("mean must be a single finite number.")
    }
    if (length(sd) != 1L || !is.numeric(sd) || !is.finite(sd) || sd <= 0) {
      stop("sd must be a single finite positive number.")
    }
    if (n > 150L) {
      # beyond this the outer nodes sit so far in the tails that the normal
      # density at them underflows
      stop("gauss-hermite supports at most 150 nodes.")
    }
  } else {
    if (length(interval) != 2L || !is.numeric(interval) || any(!is.finite(interval)) ||
        interval[1L] >= interval[2L]) {
      stop("interval must be a numeric vector c(a, b) of finite numbers with a < b.")
    }
  }
  if (!is.null(states)) {
    states <- as.character(states)
    if (length(states) != n || anyNA(states) || anyDuplicated(states) || any(!nzchar(states))) {
      stop("states must be size distinct, non-empty names.")
    }
  }

  # Step 2: nodes and weights of the rule.
  rule <- .quadratureRule(quadrature, n, interval, mean, sd)
  x <- rule$nodes

  # Step 3: kernel values K[i, j] = f(x_i, x_j), with i the departure node.
  # In a column-major matrix element (i, j) has index (j - 1) n + i, so the
  # departure nodes cycle fastest.
  from <- rep(x, times = n)
  to <- rep(x, each = n)
  values <- tryCatch(f(from, to), error = function(e) NULL)
  if (!is.numeric(values) || length(values) != n * n) {
    # not vectorized: one call per pair
    values <- tryCatch(mapply(function(u, v) as.numeric(f(u, v))[1L], from, to),
                       error = function(e) stop("f could not be evaluated: ", conditionMessage(e)))
  }
  K <- matrix(as.numeric(values), n, n)
  if (any(!is.finite(K)) || any(K < 0)) {
    stop("f must return finite, non-negative values at the nodes.")
  }

  # Step 4: weight column j by w_j and normalize the rows.
  P <- K * matrix(rule$weights, n, n, byrow = TRUE)
  rowTotals <- rowSums(P)
  if (any(rowTotals <= 0)) {
    stop("f is zero at every node for the departure node(s) ",
         paste(format(signif(x[rowTotals <= 0], 4)), collapse = ", "),
         ": change interval, size or f.")
  }
  P <- P / rowTotals

  if (is.null(states)) {
    states <- .nodeStateNames(x)
  }
  dimnames(P) <- list(states, states)
  if (is.null(name)) {
    name <- paste0("Quadrature approximation (", quadrature, ", ", n, " nodes)")
  }
  list(chain = new("markovchain", states = states, byrow = TRUE, transitionMatrix = P,
                   name = name),
       states = stats::setNames(x, states),
       weights = stats::setNames(rule$weights, states))
}

#' Discretize an AR(1) process into a Markov chain (Tauchen-Hussey method)
#'
#' Approximates the stationary first-order autoregressive process
#' \deqn{y_t = (1-\rho)\alpha + \rho y_{t-1} + \varepsilon_t, \qquad
#'   \varepsilon_t \overset{\mathrm{iid}}{\sim} \mathcal N(0,\sigma^2)}
#' by a finite-state Markov chain on Gauss-Hermite nodes, following
#' Tauchen and Hussey (1991), with the choice of the standard deviation of
#' the nodes proposed by Floden (2008) as an option.
#'
#' @inheritParams tauchen
#' @param sigmaWeight How the standard deviation \eqn{\hat\sigma} of the
#'   normal distribution defining the nodes is chosen:
#'   \code{"conditional"} (the default, Tauchen and Hussey's original
#'   choice) uses \eqn{\hat\sigma = \sigma}; \code{"unconditional"} uses
#'   \eqn{\hat\sigma = \sigma_y = \sigma/\sqrt{1-\rho^2}}; \code{"floden"}
#'   uses Floden's compromise \eqn{\hat\sigma = w\sigma + (1-w)\sigma_y}
#'   with \eqn{w = 1/2 + \rho/4}.
#'
#' @return A named list with two elements, \code{chain} and \code{states},
#'   in the same form as returned by \code{\link{tauchen}}.
#'
#' @details
#' The nodes \eqn{y_1, \ldots, y_n} and probability weights
#' \eqn{\omega_1, \ldots, \omega_n} are those of the \eqn{n}-point
#' Gauss-Hermite rule for \eqn{N(\alpha, \hat\sigma^2)}, and
#' \deqn{P_{ij} \propto \omega_j \frac{\phi(y_j;\ (1-\rho)\alpha + \rho y_i,\
#'   \sigma)}{\phi(y_j;\ \alpha,\ \hat\sigma)},}
#' normalized by rows, with \eqn{\phi(\cdot; m, s)} the normal density. It
#' is \code{\link{fitFunction}} with the AR(1) transition density and
#' \code{quadrature = "gauss-hermite"}.
#'
#' The method is accurate for moderate persistence, but with
#' \eqn{\hat\sigma = \sigma} the nodes are too concentrated when \eqn{\rho}
#' is close to one and the chain understates the unconditional variance;
#' with \eqn{\hat\sigma = \sigma_y} it overstates the persistence and
#' understates the conditional variance. Floden (2008) shows that his
#' weighted choice is a good compromise; for very persistent processes
#' \code{\link{rouwenhorst}} is usually preferable (Kopecky and Suen 2010).
#' The vignette reproduces the conditional and unconditional standard
#' deviations of Floden's Table 1.
#'
#' @references
#' Tauchen, G. and Hussey, R. (1991). Quadrature-based methods for
#' obtaining approximate solutions to nonlinear asset pricing models.
#' \emph{Econometrica}, 59(2), 371-396. \doi{10.2307/2938261}
#'
#' Floden, M. (2008). A note on the accuracy of Markov-chain approximations
#' to highly persistent AR(1) processes. \emph{Economics Letters}, 99(3),
#' 516-520. \doi{10.1016/j.econlet.2007.09.040}
#'
#' Kopecky, K. A. and Suen, R. M. H. (2010). Finite state Markov-chain
#' approximations to highly persistent processes. \emph{Review of Economic
#' Dynamics}, 13(3), 701-714.
#'
#' @seealso \code{\link{fitFunction}}, \code{\link{tauchen}},
#'   \code{\link{rouwenhorst}}
#'
#' @examples
#' out <- tauchenHussey(alpha = 0, sigma = sqrt(0.03), rho = 0.95, size = 5,
#'                      sigmaWeight = "floden")
#' out$chain
#' # unconditional standard deviation of the chain (Floden 2008, Table 1: 0.4792)
#' pi <- as.numeric(steadyStates(out$chain))
#' sqrt(sum(pi * (out$states - sum(pi * out$states))^2))
#'
#' @export
tauchenHussey <- function(alpha, sigma, rho, size,
                          sigmaWeight = c("conditional", "floden", "unconditional")) {
  .checkAR1Args(alpha, sigma, rho, size)
  sigmaWeight <- match.arg(sigmaWeight)
  sigmaY <- sigma / sqrt(1 - rho^2)
  w <- 1 / 2 + rho / 4
  sigmaHat <- switch(sigmaWeight,
                     conditional = sigma,
                     unconditional = sigmaY,
                     floden = w * sigma + (1 - w) * sigmaY)
  out <- fitFunction(function(x, y) stats::dnorm(y, (1 - rho) * alpha + rho * x, sigma),
                     size = size, quadrature = "gauss-hermite", mean = alpha, sd = sigmaHat,
                     name = paste0("Tauchen-Hussey AR(1) Approximation (rho = ", rho, ")"))
  list(chain = out$chain, states = out$states)
}
