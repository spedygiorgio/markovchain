context("fitFunction() and tauchenHussey()")

test_that("each quadrature rule integrates what it should exactly", {
  rule <- markovchain:::.quadratureRule
  # Gauss-Legendre: polynomials of degree 2n - 1
  r <- rule("gauss-legendre", 6, c(-1, 2))
  expect_equal(sum(r$weights * r$nodes^11), (2^12 - 1) / 12, tolerance = 1e-12)
  expect_equal(r$nodes, rev(1 - r$nodes), tolerance = 1e-14)
  # Fejer's first rule: polynomials of degree n - 1
  r <- rule("fejer", 6, c(0, 1))
  expect_equal(sum(r$weights * r$nodes^5), 1 / 6, tolerance = 1e-14)
  # Gauss-Chebyshev: q(t) / sqrt(1 - t^2), q of degree 2n - 1
  r <- rule("gauss-chebyshev", 6, c(-1, 1))
  expect_equal(sum(r$weights * r$nodes^10 / sqrt(1 - r$nodes^2)), pi * 945 / 3840,
               tolerance = 1e-13)
  # Gauss-Hermite: moments of N(mean, sd^2) up to degree 2n - 1
  r <- rule("gauss-hermite", 6, mean = 1, sd = 2)
  p <- r$weights * dnorm(r$nodes, 1, 2)
  expect_equal(sum(p), 1, tolerance = 1e-14)
  expect_equal(sum(p * r$nodes^4), 1 + 6 * 4 + 3 * 16, tolerance = 1e-12)
  # Newton-Cotes rules
  r <- rule("simpson", 5, c(0, 1))
  expect_equal(r$weights, c(1, 4, 2, 4, 1) / 12)
  expect_equal(sum(r$weights * r$nodes^3), 1 / 4)
  r <- rule("trapezoid", 4, c(0, 3))
  expect_equal(r$weights, c(0.5, 1, 1, 0.5))
  r <- rule("midpoint", 4, c(1, 3))
  expect_equal(r$nodes, c(1.25, 1.75, 2.25, 2.75))
  for (q in c("gauss-legendre", "gauss-chebyshev", "fejer", "gauss-hermite", "trapezoid",
              "simpson", "midpoint")) {
    expect_true(all(rule(q, 9, c(0, 1))$weights > 0), info = q)
  }
})

test_that("fitFunction() gives the matrices of PyDTMC 10.0.0's fit_function", {
  # PyDTMC: MarkovChain.fit_function(rule, 5 states, f, (0, 2)) with
  # f = exp(-(y - 0.6 x)^2 / 0.08); matrices by rows.
  pydtmc <- list(
  "gauss-legendre" = c(0.791140802418, 0.208830908905, 2.8288674867e-05, 1.92881035561e-12, 2.13055632504e-19,
    0.332066922115, 0.666174076404, 0.00175899914307, 2.33777455145e-09, 1.96257939158e-15,
    0.0207906022208, 0.81299758386, 0.166194713319, 1.71003206455e-05, 2.79826561064e-10,
    7.73844953788e-05, 0.0589842707576, 0.933499768528, 0.00743620432052, 2.37189833054e-06,
    4.82853711507e-07, 0.00279717056876, 0.862893155345, 0.133984388205, 0.000324803027985),
  "gauss-chebyshev" = c(0.664971443553, 0.335009793001, 1.87634459628e-05, 1.36861026117e-13, 5.98879683735e-21,
    0.21495996972, 0.783959733519, 0.00108029656684, 1.93867195515e-10, 6.14109440548e-17,
    0.00886909221109, 0.795810152939, 0.195314512113, 6.24268818941e-06, 4.86527317626e-11,
    1.00747941818e-05, 0.0222413236042, 0.972213091341, 0.00553444904338, 1.06121681069e-06,
    5.0985184664e-08, 0.000814798273689, 0.876284393296, 0.122730399053, 0.000170358391711),
  "simpson-rule" = c(0.850517191312, 0.149476469521, 6.33916517782e-06, 2.07592081784e-12, 1.64043490343e-22,
    0.1178347639, 0.880577285202, 0.00158792878761, 2.21112881971e-08, 7.42961944791e-17,
    0.00291427282627, 0.92603746292, 0.0710062222117, 4.2042035774e-05, 6.00676399263e-12,
    1.70431791283e-05, 0.230278432925, 0.75080200421, 0.0189024048498, 1.14836037663e-07,
    6.04182869392e-09, 0.0034711659049, 0.481229050205, 0.515166697715, 0.000133080133065),
  "trapezoid-rule" = c(0.919217915692, 0.0807752330944, 6.85121272269e-06, 1.12180190289e-12, 1.77294141514e-22,
    0.210527738395, 0.786635192367, 0.00283704948628, 1.97524030391e-08, 1.32740197183e-16,
    0.00542735302979, 0.86229610776, 0.132237391016, 3.91481827268e-05, 1.11866083576e-11,
    1.94688057979e-05, 0.131526109546, 0.857657969948, 0.0107963205211, 1.31179781602e-07,
    8.15712585808e-09, 0.00234322575125, 0.649711557204, 0.347765536234, 0.000179672653697))
  ours <- c("gauss-legendre" = "gauss-legendre", "gauss-chebyshev" = "fejer",
            "simpson-rule" = "simpson", "trapezoid-rule" = "trapezoid")
  kernel <- function(x, y) exp(-(y - 0.6 * x)^2 / 0.08)
  for (q in names(ours)) {
    P <- fitFunction(kernel, 5, ours[[q]], interval = c(0, 2))$chain@transitionMatrix
    expect_equal(as.numeric(t(P)), pydtmc[[q]], tolerance = 1e-9, info = q)
  }
})

test_that("fitFunction() returns nodes, weights and a valid chain", {
  kernel <- function(x, y) dnorm(y, mean = 0.6 * x + 0.2, sd = 0.15)
  out <- fitFunction(kernel, size = 7, interval = c(0, 1))
  expect_s4_class(out$chain, "markovchain")
  expect_equal(unname(rowSums(out$chain@transitionMatrix)), rep(1, 7))
  expect_identical(names(out$states), states(out$chain))
  expect_identical(names(out$weights), states(out$chain))
  expect_equal(sum(out$weights), 1, tolerance = 1e-14) # length of [0, 1]
  # a kernel that is not vectorized gives the same chain
  scalar <- function(x, y) { stopifnot(length(x) == 1L); kernel(x, y) }
  expect_equal(fitFunction(scalar, size = 7)$chain@transitionMatrix,
               out$chain@transitionMatrix)
  # user-supplied names and chain name
  named <- fitFunction(kernel, 3, states = c("low", "mid", "high"), name = "test")
  expect_identical(states(named$chain), c("low", "mid", "high"))
  expect_identical(named$chain@name, "test")
})

test_that("fitFunction() converges to the continuous expectation", {
  # On [0, 1] a kernel proportional to exp(-(y - m(x))^2 / (2 s^2)) truncated
  # to [0, 1]; with many nodes the chain's conditional mean must approach the
  # truncated normal mean.
  m <- function(x) 0.3 + 0.4 * x
  s <- 0.2
  kernel <- function(x, y) dnorm(y, m(x), s)
  truncMean <- function(x) {
    a <- (0 - m(x)) / s; b <- (1 - m(x)) / s
    m(x) + s * (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
  }
  err <- sapply(c(5, 10, 20), function(n) {
    out <- fitFunction(kernel, n)
    max(abs(out$chain@transitionMatrix %*% out$states - truncMean(out$states)))
  })
  expect_lt(err[3], 1e-10)
  expect_true(all(diff(err) < 0))
})

test_that("tauchenHussey() reproduces Floden (2008), Table 1", {
  moments <- function(out) {
    P <- out$chain@transitionMatrix
    z <- out$states
    p <- as.numeric(steadyStates(out$chain))
    m <- sum(p * z)
    condVar <- P %*% z^2 - (P %*% z)^2
    c(sqrt(sum(p * condVar)), sqrt(sum(p * (z - m)^2)))
  }
  # conditional and unconditional standard deviations implied by the chain
  hsz <- list( # rho = 0.95, sigma_eps^2 = 0.030
    conditional   = list(`5` = c(0.1576, 0.3275), `9` = c(0.1670, 0.4303)),
    unconditional = list(`5` = c(0.0101, 0.5622), `9` = c(0.0561, 0.5556)),
    floden        = list(`5` = c(0.1410, 0.4792), `9` = c(0.1692, 0.5407)))
  for (s in names(hsz)) for (n in names(hsz[[s]])) {
    out <- tauchenHussey(0, sqrt(0.03), 0.95, as.integer(n), sigmaWeight = s)
    expect_equal(round(moments(out), 4), hsz[[s]][[n]], info = paste(s, n))
  }
  aiyagari <- list( # rho = 0.6, sigma_eps^2 = 0.013
    conditional = list(`5` = c(0.1137, 0.1418), `9` = c(0.1140, 0.1425)),
    floden      = list(`5` = c(0.1139, 0.1424), `9` = c(0.1140, 0.1425)))
  for (s in names(aiyagari)) for (n in names(aiyagari[[s]])) {
    out <- tauchenHussey(0, sqrt(0.013), 0.6, as.integer(n), sigmaWeight = s)
    expect_equal(round(moments(out), 4), aiyagari[[s]][[n]], info = paste(s, n))
  }
})

test_that("tauchenHussey() special cases", {
  # rho = 0: every row is the Gauss-Hermite probability vector
  out <- tauchenHussey(alpha = 2, sigma = 0.5, rho = 0, size = 4)
  r <- markovchain:::.quadratureRule("gauss-hermite", 4, mean = 2, sd = 0.5)
  probs <- r$weights * dnorm(r$nodes, 2, 0.5)
  for (i in 1:4) {
    expect_equal(unname(out$chain@transitionMatrix[i, ]), probs, tolerance = 1e-13)
  }
  expect_equal(unname(out$states), r$nodes)
  # nodes centred on alpha
  expect_equal(mean(tauchenHussey(5, 1, 0.5, 5)$states), 5)
  expect_error(tauchenHussey(0, 1, 1, 5), "rho")
  expect_error(tauchenHussey(0, 1, 0.5, 5, sigmaWeight = "x"))
})

test_that("fitFunction() rejects invalid input", {
  k <- function(x, y) exp(-(x - y)^2)
  expect_error(fitFunction("f", 5), "function")
  expect_error(fitFunction(k, 1), "size")
  expect_error(fitFunction(k, 4, "simpson"), "odd")
  expect_error(fitFunction(k, 5, interval = c(1, 0)), "interval")
  expect_error(fitFunction(k, 5, "gauss-hermite", sd = 0), "sd")
  expect_error(fitFunction(k, 200, "gauss-hermite"), "150")
  expect_error(fitFunction(function(x, y) x - y, 5), "non-negative")
  expect_error(fitFunction(function(x, y) ifelse(x > 0.5, 0, 1), 5), "zero at every node")
  expect_error(fitFunction(k, 3, states = c("a", "a", "b")), "distinct")
})
