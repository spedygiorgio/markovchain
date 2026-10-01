library(testthat)
library(markovchain)

st <- c("a", "b", "c")
P <- matrix(c(0.7, 0.2, 0.1,
              0.3, 0.4, 0.3,
              0.2, 0.45, 0.35),
            byrow = TRUE, nrow = 3, dimnames = list(st, st))
mc <- new("markovchain", states = st, transitionMatrix = P, name = "W")
mcCol <- new("markovchain", states = st, transitionMatrix = t(P),
             byrow = FALSE, name = "W")

test_that("redistribute() reproduces initial %*% P^t, including t = 0", {
  init <- c(0.2, 0.5, 0.3)
  traj <- redistribute(mc, 6, initial = init)
  expect_equal(dim(traj), c(7L, 3L))
  expect_equal(colnames(traj), st)
  expect_equal(rownames(traj), as.character(0:6))
  expect_equal(unname(traj[1, ]), init)
  Pt <- diag(3)
  for (t in 1:6) {
    Pt <- Pt %*% P
    expect_equal(unname(traj[t + 1, ]), as.numeric(init %*% Pt),
                 tolerance = 1e-12)
  }
  expect_equal(unname(rowSums(traj)), rep(1, 7), tolerance = 1e-12)
})

test_that("redistribute() defaults to the uniform distribution and converges", {
  traj <- redistribute(mc, 200)
  expect_equal(unname(traj[1, ]), rep(1 / 3, 3))
  expect_equal(traj[201, ], steadyStates(mc)[1, ], tolerance = 1e-10)
})

test_that("redistribute() accepts a state name and a named vector", {
  expect_equal(redistribute(mc, 3, "b"),
               redistribute(mc, 3, c(0, 1, 0)), ignore_attr = TRUE)
  expect_equal(redistribute(mc, 3, c(c = 0.5, b = 0, a = 0.5)),
               redistribute(mc, 3, c(0.5, 0, 0.5)))
})

test_that("redistribute(lastOnly = TRUE) returns the final named distribution", {
  last <- redistribute(mc, 5, "a", lastOnly = TRUE)
  expect_named(last, st)
  expect_equal(last, redistribute(mc, 5, "a")["5", ])
  expect_equal(unname(last), as.numeric((mc^5)@transitionMatrix["a", ]),
               tolerance = 1e-12)
})

test_that("redistribute() handles column-stochastic chains", {
  expect_equal(redistribute(mcCol, 5, "a"), redistribute(mc, 5, "a"))
})

test_that("redistribute() validates its arguments", {
  expect_error(redistribute(mc, -1), "steps")
  expect_error(redistribute(mc, 1.5), "steps")
  expect_error(redistribute(mc, c(1, 2)), "steps")
  expect_error(redistribute(mc, 2, c(0.5, 0.6)), "initial")
  expect_error(redistribute(mc, 2, c(0.5, 0.6, 0.6)), "initial")
  expect_error(redistribute(mc, 2, c(-0.5, 1, 0.5)), "initial")
  expect_error(redistribute(mc, 2, "z"), "initial")
})

test_that("topologicalEntropy() is log of the Perron root of the support", {
  expect_equal(topologicalEntropy(mc), log2(3), tolerance = 1e-12)
  expect_equal(topologicalEntropy(mc, base = exp(1)), log(3), tolerance = 1e-12)
  expect_equal(topologicalEntropy(mcCol), topologicalEntropy(mc))
  cyc <- new("markovchain", states = st,
             transitionMatrix = matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), 3,
                                       byrow = TRUE,
                                       dimnames = list(st, st)))
  expect_equal(topologicalEntropy(cyc), 0)
  # Golden-mean shift: forbidden "bb", Perron root of [[1,1],[1,0]] is phi.
  s2 <- c("a", "b")
  gm <- new("markovchain", states = s2,
            transitionMatrix = matrix(c(0.5, 0.5, 1, 0), 2, byrow = TRUE,
                                      dimnames = list(s2, s2)))
  expect_equal(topologicalEntropy(gm, base = exp(1)), log((1 + sqrt(5)) / 2),
               tolerance = 1e-12)
})

test_that("entropy rate never exceeds topological entropy; normalized is in [0,1]", {
  expect_lte(entropyRate(mc, 2), topologicalEntropy(mc, 2) + 1e-12)
  ne <- normalizedEntropyRate(mc)
  expect_gte(ne, 0)
  expect_lte(ne, 1)
  expect_equal(ne, entropyRate(mc, 2) / topologicalEntropy(mc, 2),
               tolerance = 1e-12)
  full <- new("markovchain", states = st,
              transitionMatrix = matrix(1 / 3, 3, 3, dimnames = list(st, st)))
  expect_equal(normalizedEntropyRate(full), 1, tolerance = 1e-12)
  cyc <- new("markovchain", states = st,
             transitionMatrix = matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), 3,
                                       byrow = TRUE,
                                       dimnames = list(st, st)))
  expect_equal(normalizedEntropyRate(cyc), 0)
})

test_that("relaxationTime() is 1/spectralGap and Inf for periodic chains", {
  expect_equal(relaxationTime(mc), 1 / spectralGap(mc), tolerance = 1e-12)
  expect_equal(relaxationTime(mc), 1 / 0.55, tolerance = 1e-10)
  s2 <- c("a", "b")
  flip <- new("markovchain", states = s2,
              transitionMatrix = matrix(c(0, 1, 1, 0), 2, byrow = TRUE,
                                        dimnames = list(s2, s2)))
  expect_equal(relaxationTime(flip), Inf)
  expect_error(relaxationTime(new("markovchain", states = s2,
    transitionMatrix = matrix(c(1, 0, 0, 1), 2, dimnames = list(s2, s2)))),
    "irreducible")
})

test_that("autoplot() supports the eigenvalues and flow types", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(ggplot2::autoplot(mc, type = "eigenvalues"), "ggplot")
  expect_s3_class(ggplot2::autoplot(mc, type = "flow", steps = 5,
                                    initial = "a"), "ggplot")
  expect_s3_class(ggplot2::autoplot(mcCol, type = "flow"), "ggplot")
  expect_error(ggplot2::autoplot(mc, type = "nope"))
  expect_error(ggplot2::autoplot(mc, type = "flow", steps = -1), "steps")
})

test_that("autoplot(type = 'comparison') compares chains by state name", {
  skip_if_not_installed("ggplot2")
  lazy <- lazyChain(mc, alpha = 0.5)
  lazy@name <- "Lazy"
  p <- ggplot2::autoplot(mc, type = "comparison", other = lazy)
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_equal(nlevels(d$chain), 2L)
  expect_equal(nrow(d), 2L * 9L)
  # plotted values are the transition probabilities of each chain
  at <- function(ch, f, t) d$probability[d$chain == ch & d$from == f & d$to == t]
  expect_equal(at("W", "a", "b"), P["a", "b"])
  expect_equal(at("Lazy", "a", "a"), 0.5 + 0.5 * P["a", "a"])

  # a permuted copy of the same chain on the same states gives identical values
  perm <- c("c", "a", "b")
  mcPerm <- new("markovchain", states = perm,
                transitionMatrix = P[perm, perm], name = "Perm")
  d2 <- ggplot2::autoplot(mc, type = "comparison", other = mcPerm)$data
  expect_equal(d2$probability[d2$chain == "W"], d2$probability[d2$chain == "Perm"])

  # column-stochastic storage and a named list
  p3 <- ggplot2::autoplot(mc, type = "comparison",
                          other = list(Col = mcCol, Lazy = lazy))
  expect_equal(levels(p3$data$chain), c("W", "Col", "Lazy"))
})

test_that("autoplot(type = 'comparison', what = 'stationary') shows steady states", {
  skip_if_not_installed("ggplot2")
  lazy <- lazyChain(mc, alpha = 0.5)
  p <- ggplot2::autoplot(mc, type = "comparison", other = lazy,
                         what = "stationary")
  expect_s3_class(p, "ggplot")
  d <- p$data
  # laziness does not change the stationary distribution
  expect_equal(d$probability[1:3], d$probability[4:6], tolerance = 1e-10)
  expect_equal(d$probability[1:3], as.numeric(steadyStates(mc)[1, ]),
               tolerance = 1e-10)
})

test_that("autoplot(type = 'comparison') validates its inputs", {
  skip_if_not_installed("ggplot2")
  expect_error(ggplot2::autoplot(mc, type = "comparison"), "other")
  expect_error(ggplot2::autoplot(mc, type = "comparison", other = 1), "other")
  s2 <- c("a", "b")
  small <- new("markovchain", states = s2,
               transitionMatrix = matrix(c(.5, .5, .5, .5), 2,
                                         dimnames = list(s2, s2)))
  expect_error(ggplot2::autoplot(mc, type = "comparison", other = small),
               "same set of states")
  absorb <- new("markovchain", states = st,
                transitionMatrix = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3,
                                          dimnames = list(st, st)),
                name = "Absorb")
  expect_error(ggplot2::autoplot(mc, type = "comparison", other = absorb,
                                 what = "stationary"), "irreducible")
})
