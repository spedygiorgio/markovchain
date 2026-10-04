context("Parallel RcppParallel path: thread safety, reproducibility, byrow consistency")

.mcListFixture <- function() {
  sn <- c("a", "b", "c")
  mcA <- new("markovchain", states = sn,
             transitionMatrix = matrix(c(.6, .3, .1,
                                         .2, .7, .1,
                                         .1, .1, .8),
                                       3, byrow = TRUE,
                                       dimnames = list(sn, sn)))
  mcB <- new("markovchain", states = sn,
             transitionMatrix = matrix(c(.3, .4, .3,
                                         .5, .4, .1,
                                         .2, .3, .5),
                                       3, byrow = TRUE,
                                       dimnames = list(sn, sn)))
  list(sn = sn, mcA = mcA, mcB = mcB,
       mlRow = new("markovchainList", markovchains = list(mcA, mcB)))
}

test_that("parallel rmarkovchain is reproducible under set.seed()", {
  fx <- .mcListFixture()
  set.seed(42)
  r1 <- rmarkovchain(500, fx$mlRow, what = "list",
                     parallel = TRUE, num.cores = 2, include.t0 = TRUE)
  set.seed(42)
  r2 <- rmarkovchain(500, fx$mlRow, what = "list",
                     parallel = TRUE, num.cores = 2, include.t0 = TRUE)
  expect_identical(r1, r2)

  # Same seed, different thread count should still give identical sequences:
  # the per-sequence seeding scheme makes the output independent of the TBB
  # sub-range split and so of the thread count. Protects against future
  # regressions that would re-introduce per-sub-range seeding.
  set.seed(42)
  r3 <- rmarkovchain(500, fx$mlRow, what = "list",
                     parallel = TRUE, num.cores = 4, include.t0 = TRUE)
  expect_identical(r1, r3)
})

test_that("parallel rmarkovchain respects the byrow slot (closes the parallel branch of #148)", {
  fx <- .mcListFixture()
  mcA_col <- new("markovchain", states = fx$sn, byrow = FALSE,
                 transitionMatrix = t(fx$mcA@transitionMatrix))
  mcB_col <- new("markovchain", states = fx$sn, byrow = FALSE,
                 transitionMatrix = t(fx$mcB@transitionMatrix))
  mlCol <- new("markovchainList", markovchains = list(mcA_col, mcB_col))

  set.seed(1); xRow <- rmarkovchain(1000, fx$mlRow, what = "list",
                                    parallel = TRUE, num.cores = 2,
                                    include.t0 = TRUE)
  set.seed(1); xCol <- rmarkovchain(1000, mlCol, what = "list",
                                    parallel = TRUE, num.cores = 2,
                                    include.t0 = TRUE)
  expect_identical(xRow, xCol)
})

test_that("parallel rmarkovchain tracks the theoretical transition probabilities", {
  fx <- .mcListFixture()
  set.seed(7)
  M <- rmarkovchain(20000, fx$mlRow, what = "matrix",
                    parallel = TRUE, num.cores = 2, include.t0 = TRUE)
  # Empirical frequencies of the first transition should match mcA.
  emp_a <- prop.table(table(factor(M[M[, 1] == "a", 2], levels = fx$sn)))
  expect_equal(as.numeric(emp_a),
               as.numeric(fx$mcA@transitionMatrix["a", ]),
               tolerance = 0.02)
  # Empirical frequencies of the second transition should match mcB.
  emp_b <- prop.table(table(factor(M[M[, 2] == "b", 3], levels = fx$sn)))
  expect_equal(as.numeric(emp_b),
               as.numeric(fx$mcB@transitionMatrix["b", ]),
               tolerance = 0.02)
})

test_that("parallel bootstrap is reproducible under set.seed()", {
  fx <- .mcListFixture()
  set.seed(21)
  seq_data <- rmarkovchain(600, fx$mcA, t0 = "a")

  set.seed(99)
  f1 <- markovchainFit(data = seq_data, method = "bootstrap", nboot = 100,
                       parallel = TRUE, num.cores = 2, confidencelevel = 0.95)
  set.seed(99)
  f2 <- markovchainFit(data = seq_data, method = "bootstrap", nboot = 100,
                       parallel = TRUE, num.cores = 2, confidencelevel = 0.95)
  expect_equal(f1$estimate@transitionMatrix,
               f2$estimate@transitionMatrix)
  # Changing num.cores must not change the output under the same seed.
  set.seed(99)
  f3 <- markovchainFit(data = seq_data, method = "bootstrap", nboot = 100,
                       parallel = TRUE, num.cores = 4, confidencelevel = 0.95)
  expect_equal(f1$estimate@transitionMatrix,
               f3$estimate@transitionMatrix)
})

test_that(".mcDesiredThreads honours options(), env vars and num.cores", {
  f <- markovchain:::.mcDesiredThreads

  # Snapshot and restore the options/env vars we touch, without adding a
  # dependency on withr (which is not in Suggests).
  opts_bak <- options(RcppParallel.numThreads = NULL, Ncpus = NULL)
  env_bak <- Sys.getenv(c("RCPP_PARALLEL_NUM_THREADS", "OMP_NUM_THREADS"),
                       names = TRUE, unset = NA)
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "", OMP_NUM_THREADS = "")
  on.exit({
    options(opts_bak)
    for (nm in names(env_bak)) {
      if (is.na(env_bak[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        args <- list(env_bak[[nm]]); names(args) <- nm
        do.call(Sys.setenv, args)
      }
    }
  })

  # Explicit num.cores wins over everything.
  expect_identical(f(num.cores = 3L), 3L)
  expect_identical(f(num.cores = "5"), 5L)

  options(RcppParallel.numThreads = 4L, Ncpus = 7L)
  expect_identical(f(), 4L)

  options(RcppParallel.numThreads = NULL, Ncpus = 6L)
  expect_identical(f(), 6L)

  options(RcppParallel.numThreads = NULL, Ncpus = NULL)
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "8")
  expect_identical(f(), 8L)
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "", OMP_NUM_THREADS = "5")
  expect_identical(f(), 5L)

  # Default: no option, no env var -> cap at 2 cores per CRAN policy.
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "", OMP_NUM_THREADS = "")
  expect_lte(f(), 2L)
  expect_gte(f(), 1L)
})
