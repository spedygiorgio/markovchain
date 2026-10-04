context("sanitize = \"absorbing\" for unobserved states (#213)")

# "a" and "b" alternate; "c" and "d" only exist through possibleStates and so
# have no observed outgoing transition at all.
.seq <- c("a", "b", "a", "b", "a", "b")
.possible <- c("a", "b", "c", "d")

test_that("the historical logical values of sanitize are unchanged", {
  zero <- createSequenceMatrix(.seq, TRUE, FALSE, .possible)
  expect_equal(unname(rowSums(zero)), c(1, 1, 0, 0))

  uniform <- createSequenceMatrix(.seq, TRUE, TRUE, .possible)
  expect_equal(as.numeric(uniform["c", ]), rep(0.25, 4))
  expect_equal(as.numeric(uniform["d", ]), rep(0.25, 4))

  # the counts matrix keeps the all-ones row of earlier versions
  counts <- createSequenceMatrix(.seq, FALSE, TRUE, .possible)
  expect_equal(as.numeric(counts["c", ]), rep(1, 4))
})

test_that("\"uniform\" is an exact synonym of TRUE", {
  expect_identical(createSequenceMatrix(.seq, TRUE, "uniform", .possible),
                   createSequenceMatrix(.seq, TRUE, TRUE, .possible))
  expect_identical(createSequenceMatrix(.seq, FALSE, "uniform", .possible),
                   createSequenceMatrix(.seq, FALSE, TRUE, .possible))
})

test_that("\"absorbing\" makes the unobserved states absorbing, as asked in #213", {
  probs <- createSequenceMatrix(.seq, TRUE, "absorbing", .possible)
  # a stochastic matrix, which sanitize = FALSE does not give
  expect_equal(unname(rowSums(probs)), rep(1, 4))
  # the unobserved states stay where they are
  expect_equal(as.numeric(probs["c", ]), c(0, 0, 1, 0))
  expect_equal(as.numeric(probs["d", ]), c(0, 0, 0, 1))
  # the observed rows are untouched
  expect_equal(as.numeric(probs["a", ]), c(0, 1, 0, 0))
  expect_equal(as.numeric(probs["b", ]), c(1, 0, 0, 0))

  # on counts, the diagonal count is the analogue of the all-ones row
  counts <- createSequenceMatrix(.seq, FALSE, "absorbing", .possible)
  expect_equal(as.numeric(counts["c", ]), c(0, 0, 1, 0))
  expect_equal(as.numeric(counts["a", ]),
               as.numeric(createSequenceMatrix(.seq, FALSE, FALSE, .possible)["a", ]))
})

test_that("absorbing and uniform agree when every state is observed", {
  # With no empty row there is nothing to sanitize, so all three modes must
  # return exactly the same matrix.
  full <- c("a", "b", "c", "a", "c", "b", "a")
  base <- createSequenceMatrix(full, TRUE, FALSE)
  expect_identical(createSequenceMatrix(full, TRUE, TRUE), base)
  expect_identical(createSequenceMatrix(full, TRUE, "absorbing"), base)
})

test_that("sanitize is validated", {
  expect_error(createSequenceMatrix(.seq, TRUE, "nonsense", .possible))
  expect_error(createSequenceMatrix(.seq, TRUE, NA, .possible), "sanitize")
  expect_error(createSequenceMatrix(.seq, TRUE, c(TRUE, FALSE), .possible), "sanitize")
  expect_error(createSequenceMatrix(.seq, TRUE, 1, .possible), "sanitize")
  expect_error(createSequenceMatrix(.seq, TRUE, NA_character_, .possible), "sanitize")
  # match.arg allows an unambiguous abbreviation
  expect_identical(createSequenceMatrix(.seq, TRUE, "absorb", .possible),
                   createSequenceMatrix(.seq, TRUE, "absorbing", .possible))
})

test_that("markovchainFit(sanitize = \"absorbing\") works for every method", {
  for (method in c("mle", "laplace", "bootstrap", "map")) {
    fit <- markovchainFit(.seq, method = method, possibleStates = .possible,
                          sanitize = "absorbing", nboot = 5)
    P <- fit$estimate@transitionMatrix
    expect_equal(unname(rowSums(P)), rep(1, 4), tolerance = 1e-12,
                 info = method)
    expect_equal(as.numeric(P["c", ]), c(0, 0, 1, 0), info = method)
    expect_equal(as.numeric(P["d", ]), c(0, 0, 0, 1), info = method)
  }
})

test_that("markovchainFit(sanitize = \"absorbing\") honours byrow = FALSE", {
  fit <- markovchainFit(.seq, method = "mle", byrow = FALSE,
                        possibleStates = .possible, sanitize = "absorbing")
  P <- fit$estimate@transitionMatrix
  expect_false(fit$estimate@byrow)
  # column-stochastic: the outgoing distribution of a state is its column
  expect_equal(unname(colSums(P)), rep(1, 4), tolerance = 1e-12)
  expect_equal(as.numeric(P[, "c"]), c(0, 0, 1, 0))
  expect_equal(as.numeric(P[, "d"]), c(0, 0, 0, 1))
})

test_that("markovchainFit keeps the earlier behaviour for logical sanitize", {
  uniformFit <- markovchainFit(.seq, method = "mle", possibleStates = .possible,
                               sanitize = TRUE)
  expect_equal(as.numeric(uniformFit$estimate@transitionMatrix["c", ]),
               rep(0.25, 4))
  zeroFit <- markovchainFit(.seq, method = "mle", possibleStates = .possible,
                            sanitize = FALSE)
  expect_equal(as.numeric(zeroFit$estimate@transitionMatrix["c", ]),
               rep(0, 4))
})

test_that("sanitize = \"absorbing\" agrees with naming the states explicitly", {
  viaSanitize <- markovchainFit(.seq, method = "mle",
                                possibleStates = .possible,
                                sanitize = "absorbing")
  viaArgument <- markovchainFit(.seq, method = "mle",
                                possibleStates = .possible,
                                absorbingStates = c("c", "d"))
  expect_equal(viaSanitize$estimate@transitionMatrix,
               viaArgument$estimate@transitionMatrix)
})

test_that("absorbingStates stays an MLE-only feature", {
  # sanitize = "absorbing" is not restricted that way, but naming the states
  # still is, as documented.
  expect_error(
    markovchainFit(.seq, method = "laplace", possibleStates = .possible,
                   absorbingStates = "c"),
    "only with method"
  )
  expect_silent(
    markovchainFit(.seq, method = "laplace", possibleStates = .possible,
                   sanitize = "absorbing")
  )
})

test_that("a state that is observed is never turned into an absorbing one", {
  # "c" appears with an outgoing transition, so sanitize must leave it alone
  # even in absorbing mode; only "d" is empty.
  seqWithC <- c("a", "b", "c", "a", "b", "c", "a")
  P <- createSequenceMatrix(seqWithC, TRUE, "absorbing",
                            c("a", "b", "c", "d"))
  expect_equal(as.numeric(P["c", ]), c(1, 0, 0, 0))
  expect_equal(as.numeric(P["d", ]), c(0, 0, 0, 1))
})

test_that("sanitize = \"absorbing\" works on a list of sequences", {
  sequences <- list(c("a", "b", "a"), c("b", "a", "b"))
  P <- createSequenceMatrix(sequences, TRUE, "absorbing", .possible)
  expect_equal(unname(rowSums(P)), rep(1, 4))
  expect_equal(as.numeric(P["c", ]), c(0, 0, 1, 0))
  expect_equal(as.numeric(P["d", ]), c(0, 0, 0, 1))
})
