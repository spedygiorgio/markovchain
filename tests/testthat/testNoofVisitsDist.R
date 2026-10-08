context("noofVisitsDist returns the mean occupation of the first N steps (#139)")

P <- matrix(c(.4, .6, 0,
              .3, .5, .2,
              .1, .1, .8), 3, byrow = TRUE,
            dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
mc <- new("markovchain", transitionMatrix = P)

test_that("it is (1/N) times the sum of the first N matrix powers", {
  for (N in c(1, 2, 7, 25)) {
    Pk <- diag(3); S <- matrix(0, 3, 3)
    for (k in seq_len(N)) { Pk <- Pk %*% P; S <- S + Pk }
    for (s in c("a", "c")) {
      expect_equal(unname(noofVisitsDist(mc, N, s)), unname(S[match(s, rownames(P)), ] / N),
                   tolerance = 1e-12)
    }
  }
})

test_that("it is a distribution with state names, and N = 1 is one step", {
  out <- noofVisitsDist(mc, 10, "b")
  expect_named(out, c("a", "b", "c"))
  expect_equal(sum(out), 1, tolerance = 1e-12)
  expect_equal(unname(noofVisitsDist(mc, 1, "b")), unname(P["b", ]))
})

test_that("N times the result is the expected number of visits", {
  # Monte Carlo check of E[V_j(N)], time 0 excluded
  set.seed(139)
  N <- 6; reps <- 20000
  sims <- replicate(reps, {
    x <- markovchainSequence(N, mc, t0 = "a")
    table(factor(x, levels = c("a", "b", "c")))
  })
  expect_equal(unname(rowMeans(sims)), unname(N * noofVisitsDist(mc, N, "a")),
               tolerance = 0.03)
})

test_that("it converges to the stationary distribution", {
  expect_equal(unname(noofVisitsDist(mc, 5000, "a")),
               unname(as.numeric(steadyStates(mc))), tolerance = 1e-3)
})
