context("higher order models on lists of independent sequences")

.listSeqs <- function() {
  set.seed(21)
  Q <- matrix(c(.7, .2, .1, .1, .8, .1, .2, .2, .6), 3, byrow = TRUE,
              dimnames = list(letters[1:3], letters[1:3]))
  one <- function(n) {
    s <- c("a", "b")
    for (t in 3:n) {
      g <- sample(1:2, 1, prob = c(.7, .3))
      s[t] <- sample(letters[1:3], 1, prob = Q[s[t - g], ])
    }
    s
  }
  list(one(300), one(250), one(400))
}

test_that("a list with one sequence is the same as the sequence", {
  s <- .listSeqs()[[1]]
  for (o in 1:3) expect_equal(markovchain:::.seq2matHighList(list(s), o), seq2matHigh(s, o))
  expect_equal(fitHigherOrder(list(s), 2, method = "mle"), fitHigherOrder(s, 2, method = "mle"))
  expect_equal(fitMTD(list(s), 2), fitMTD(s, 2))
  expect_equal(higherOrderLogLik(list(s), fitHigherOrder(s, 2, method = "mle")),
               higherOrderLogLik(s, fitHigherOrder(s, 2, method = "mle")))
  expect_equal(assessOrder(list(s), verbose = FALSE)$statistic, assessOrder(s, verbose = FALSE)$statistic)
})

test_that("pooled lag matrices count the pairs within each sequence only", {
  seqs <- list(c("a", "b", "a", "b"), c("b", "b", "a"), c("a", "a"))
  Q1 <- markovchain:::.seq2matHighList(seqs, 1)
  # hand count of transitions from -> to: a->b 2, b->a 2, b->b 1, a->a 1
  expect_equal(unname(Q1["b", "a"]), 2 / 3)
  expect_equal(unname(Q1["a", "a"]), 1 / 3)
  expect_equal(unname(Q1["a", "b"]), 2 / 3)
  expect_equal(unname(Q1["b", "b"]), 1 / 3)
  # concatenating would add the spurious transitions b->b and a->b between sequences
  expect_false(isTRUE(all.equal(Q1, seq2matHigh(unlist(seqs), 1))))
  # lag 2: a->a and b->b in the first sequence, b->a in the second, none in the third
  Q2 <- markovchain:::.seq2matHighList(seqs, 2)
  expect_equal(unname(Q2["a", "a"]), 1)
  expect_equal(unname(Q2["a", "b"]), 1 / 2)
  expect_equal(unname(Q2["b", "b"]), 1 / 2)
  expect_equal(colSums(Q2), c(a = 1, b = 1))
})

test_that("fitHigherOrder on a list: mle weights maximize the pooled likelihood", {
  seqs <- .listSeqs()
  fit <- fitHigherOrder(seqs, 2, method = "mle")
  expect_equal(sum(fit$lambda), 1)
  expect_equal(fit$X, seq2freqProb(unlist(seqs)))
  ll <- function(l) higherOrderLogLik(seqs, list(lambda = c(l, 1 - l), Q = fit$Q))$logLik
  grid <- seq(0, 1, by = 0.005)
  best <- grid[which.max(vapply(grid, ll, numeric(1)))]
  expect_equal(unname(fit$lambda[1]), best, tolerance = 0.01)
  expect_gte(ll(fit$lambda[1]), max(vapply(grid, ll, numeric(1))) - 1e-8)
  skip_if_not_installed("Rsolnp")
  lsq <- fitHigherOrder(seqs, 2, method = "lsq")
  expect_equal(lsq$Q, fit$Q)
})

test_that("higherOrderLogLik sums over the sequences and honours start", {
  seqs <- .listSeqs()
  fit <- fitHigherOrder(seqs, 2, method = "mle")
  ll <- higherOrderLogLik(seqs, fit)
  parts <- lapply(seqs, higherOrderLogLik, fit = fit)
  expect_equal(ll$logLik, sum(vapply(parts, `[[`, numeric(1), "logLik")))
  expect_equal(ll$nobs, sum(lengths(seqs)) - 2 * length(seqs))
  ll4 <- higherOrderLogLik(seqs, fit, start = 4)
  expect_equal(ll4$nobs, sum(lengths(seqs) - 3))
  # a short sequence is skipped, not an error
  mixed <- c(seqs, list(c("a", "b")))
  expect_equal(higherOrderLogLik(mixed, fit)$logLik, ll$logLik)
  expect_error(higherOrderLogLik(list(c("a", "b")), fit), "no sequence")
})

test_that("fitMTD on a list stacks the observations of independent sequences", {
  seqs <- .listSeqs()
  fit <- fitMTD(seqs, 2)
  expect_equal(fit$nobs, sum(lengths(seqs)) - 2 * length(seqs))
  # the log-likelihood is the sum over the sequences evaluated at the pooled fit
  ll <- higherOrderLogLik(seqs, fit)
  expect_equal(ll$logLik, fit$logLikelihood, tolerance = 1e-8)
  expect_equal(ll$nobs, fit$nobs)
  # pooling cannot be mistaken for concatenation
  expect_false(isTRUE(all.equal(fit$logLikelihood, fitMTD(unlist(seqs), 2)$logLikelihood)))
  # the true model has weights (.7, .3)
  expect_equal(unname(fit$lambda[1]), 0.7, tolerance = 0.15)
  # sequences shorter than start are not used
  expect_equal(fitMTD(c(seqs, list(c("a", "b"))), 2)$logLikelihood, fit$logLikelihood)
  expect_error(fitMTD(list(c("a", "b"), c("b", "a")), 2), "at least three")
})

test_that("assessOrder on a list pools the triples within each sequence", {
  seqs <- .listSeqs()
  a <- assessOrder(seqs, verbose = FALSE)
  # reference: the same statistic from triples built by hand
  trip <- do.call(rbind, lapply(seqs, function(s) {
    n <- length(s); cbind(past = s[1:(n - 2)], pres = s[2:(n - 1)], fut = s[3:n])
  }))
  stat <- 0; dof <- 0
  for (p in unique(trip[, "pres"])) {
    m <- table(trip[trip[, "pres"] == p, "past"], trip[trip[, "pres"] == p, "fut"])
    if (nrow(m) > 1 && ncol(m) > 1) {
      stat <- stat + unname(suppressWarnings(chisq.test(m, correct = FALSE)$statistic))
      dof <- dof + (nrow(m) - 1) * (ncol(m) - 1)
    }
  }
  expect_equal(unname(a$statistic), stat, tolerance = 1e-10)
  expect_equal(unname(a$parameter), dof)
  expect_error(assessOrder(list(c("a", "b"), c("a")), verbose = FALSE), "four observations")
  expect_error(assessOrder(list(c("a", NA, "b", "a")), verbose = FALSE), "missing")
})
