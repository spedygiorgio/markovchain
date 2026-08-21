context("markovchainFit(method = \"map\") - Bayesian MAP estimation")

### Reference: with hyperparam = 1 for every cell (a flat/uniform Dirichlet(1,...,1)
### prior on each row), the MAP estimate reduces EXACTLY to the plain MLE
### n_ij / n_i+, because the mode of a Dirichlet(1,...,1) distribution's j-th
### component is (alpha_j - 1)/(sum(alpha) - k) = n_ij / n_i+. This identity
### holds regardless of the data, so it is checked directly below.

test_that("MAP with a flat (uniform) prior reduces exactly to the MLE", {
  s <- c("a", "b", "a", "a", "b", "a", "b", "b", "a", "a", "b", "a")
  hyp <- matrix(1, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))

  fit <- markovchainFit(s, method = "map", hyperparam = hyp)

  counts <- matrix(0, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  for (i in seq_len(length(s) - 1))
    counts[s[i], s[i + 1]] <- counts[s[i], s[i + 1]] + 1
  mle <- counts / rowSums(counts)

  expect_equal(unname(fit$estimate@transitionMatrix), unname(mle))
})

### Real, published data: Diaconis, P. and Rolles, S. W. W. (2006),
### "Bayesian analysis for reversible Markov chains", The Annals of
### Statistics, 34(3), 1270-1292, Table 3. The table reports the observed
### counts N_ij of each dinucleotide "ij" in a 3370-base subsequence of the
### human HLA-B gene. Row totals (620, 974, 1064, 711) match the k_a, k_c,
### k_g, k_t values the paper itself reports in equation (58), confirming
### the table was transcribed correctly.
###
### The paper's own "Example A" uses a Dirichlet(1,1,1,1) prior on the
### stationary nucleotide frequencies (a different, simpler model than the
### reversible-chain prior used in the rest of the paper). Applying that
### same flat-prior convention row-by-row to Table 3 -- exactly what
### markovchainFit(method = "map") with the package's default hyperparam
### does -- means the resulting MAP transition matrix must equal the
### row-wise MLE N_ij / N_i+, computable directly from the published table
### with no software beyond arithmetic. This gives an external, real-data
### check of the fitted values that does not depend on any number the
### package itself produced.
###
### Since markovchainFit() takes a sequence (or list of sequences) rather
### than a pre-tabulated count matrix, a sequence reproducing Table 3's
### transition counts exactly is reconstructed via an Eulerian trail
### (Hierholzer's algorithm) through the directed multigraph whose edge
### (i, j) has multiplicity N_ij. Its length (3370) matches the sequence
### length the paper reports.

.buildSequenceFromCounts <- function(counts, startState) {
  states <- rownames(counts)
  adj <- setNames(vector("list", length(states)), states)
  for (i in states) {
    dests <- character(0)
    for (j in states) if (counts[i, j] > 0) dests <- c(dests, rep(j, counts[i, j]))
    adj[[i]] <- dests
  }
  stack <- c(startState)
  path <- character(0)
  while (length(stack) > 0) {
    v <- stack[length(stack)]
    if (length(adj[[v]]) > 0) {
      u <- adj[[v]][1]
      adj[[v]] <- adj[[v]][-1]
      stack <- c(stack, u)
    } else {
      path <- c(v, path)
      stack <- stack[-length(stack)]
    }
  }
  path
}

dnaCounts <- matrix(
  c(91, 160, 261, 108,
    213, 351, 161, 249,
    251, 224, 388, 201,
    66, 239, 254, 152),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("a", "c", "g", "t"), c("a", "c", "g", "t"))
)

test_that("row totals of the published DNA table match the paper's own reported values", {
  # Diaconis & Rolles (2006), equation (58): k_a = 620, k_c = 974, k_g = 1064, k_t = 711
  expect_equal(unname(rowSums(dnaCounts)), c(620, 974, 1064, 711))
})

test_that("MAP estimate on real published DNA data matches the hand-computable MLE", {
  dnaSequence <- .buildSequenceFromCounts(dnaCounts, "t")
  expect_equal(length(dnaSequence), 3370)  # matches the sequence length reported in the paper

  hypFlat <- matrix(1, 4, 4, dimnames = list(c("a", "c", "g", "t"), c("a", "c", "g", "t")))
  fit <- markovchainFit(dnaSequence, method = "map", hyperparam = hypFlat)

  expectedMLE <- dnaCounts / rowSums(dnaCounts)
  expect_equal(unname(fit$estimate@transitionMatrix), unname(expectedMLE))
})

### Point estimate against the general (informative-prior) Dirichlet mode
### formula, hand-computed: mode_j = (alpha_j - 1) / (sum(alpha) - k), where
### alpha_j = n_ij + hyperparam_ij (Gelman et al., "Bayesian Data Analysis",
### 3rd ed., Section 3.4, gives the same Dirichlet-multinomial conjugate
### update and posterior mode formula in general form).

test_that("MAP point estimate matches the Dirichlet mode formula with an informative prior", {
  s <- c("a", "b", "a", "a", "b", "a", "b", "b", "a", "a")
  hyp <- matrix(c(2, 3, 4, 2), nrow = 2, byrow = TRUE,
               dimnames = list(c("a", "b"), c("a", "b")))

  fit <- markovchainFit(s, method = "map", hyperparam = hyp)

  counts <- matrix(0, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  for (i in seq_len(length(s) - 1))
    counts[s[i], s[i + 1]] <- counts[s[i], s[i + 1]] + 1
  alphaPost <- counts + hyp
  k <- 2
  expectedMode <- (alphaPost - 1) / (rowSums(alphaPost) - k)

  expect_equal(unname(fit$estimate@transitionMatrix), unname(expectedMode))
})

test_that("standard error matches the Beta(p, q) variance formula", {
  s <- c("a", "b", "a", "a", "b", "a", "b", "b", "a", "a")
  hyp <- matrix(c(2, 3, 4, 2), nrow = 2, byrow = TRUE,
               dimnames = list(c("a", "b"), c("a", "b")))
  fit <- markovchainFit(s, method = "map", hyperparam = hyp)

  counts <- matrix(0, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  for (i in seq_len(length(s) - 1))
    counts[s[i], s[i + 1]] <- counts[s[i], s[i + 1]] + 1
  alphaPost <- counts + hyp

  expectedSE <- matrix(0, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  for (i in 1:2) for (j in 1:2) {
    p <- alphaPost[i, j]; q <- sum(alphaPost[i, ]) - p
    expectedSE[i, j] <- sqrt(p * q / (p + q)^2 / (p + q + 1))
  }
  expect_equal(unname(fit$standardError), unname(expectedSE))
})

test_that("confidence interval is the equal-tailed Beta(p, q) interval", {
  s <- c("a", "b", "a", "a", "b", "a", "b", "b", "a", "a")
  hyp <- matrix(c(2, 3, 4, 2), nrow = 2, byrow = TRUE,
               dimnames = list(c("a", "b"), c("a", "b")))
  fit <- markovchainFit(s, method = "map", hyperparam = hyp, confidencelevel = 0.95)

  counts <- matrix(0, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  for (i in seq_len(length(s) - 1))
    counts[s[i], s[i + 1]] <- counts[s[i], s[i + 1]] + 1
  alphaPost <- counts + hyp

  p <- alphaPost["b", "a"]; q <- sum(alphaPost["b", ]) - p
  expect_equal(fit$confidenceInterval$lowerEndpointMatrix["b", "a"], qbeta(0.025, p, q),
              tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(fit$confidenceInterval$upperEndpointMatrix["b", "a"], qbeta(0.975, p, q),
              tolerance = 1e-6, check.attributes = FALSE)

  # exact nominal coverage under the fitted Beta(p, q) posterior, for any p, q
  coverage <- pbeta(fit$confidenceInterval$upperEndpointMatrix["b", "a"], p, q) -
              pbeta(fit$confidenceInterval$lowerEndpointMatrix["b", "a"], p, q)
  expect_equal(coverage, 0.95, tolerance = 1e-6)
})

test_that("standardError and confidenceInterval matrices are indexable by state name", {
  s <- c("a", "b", "a", "a", "b", "a", "b", "b", "a", "a")
  hyp <- matrix(c(2, 3, 4, 2), nrow = 2, byrow = TRUE,
               dimnames = list(c("a", "b"), c("a", "b")))
  fit <- markovchainFit(s, method = "map", hyperparam = hyp)

  expect_false(is.null(dimnames(fit$standardError)))
  expect_false(is.null(dimnames(fit$confidenceInterval$lowerEndpointMatrix)))
  expect_false(is.null(dimnames(fit$confidenceInterval$upperEndpointMatrix)))
  expect_true(is.finite(fit$standardError["a", "b"]))
})

test_that("hyperparam entries below 1 are rejected", {
  hypBad <- matrix(0.5, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  expect_error(
    markovchainFit(c("a", "b", "a", "b"), method = "map", hyperparam = hypBad)
  )
})
