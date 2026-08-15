context("Structural analysis - public textbook examples")

### "Land of Oz" weather chain (Kemeny, Snell & Thompson, "Introduction to
### Finite Mathematics", 1966; reproduced in Grinstead & Snell, "Introduction
### to Probability", 2nd ed., freely available at
### https://math.dartmouth.edu/~prob/prob/prob.pdf, chapter 11).
### This is also the exact "toy example" cited in the package's own paper
### (Spedicato, "Discrete Time Markov Chains with R", The R Journal, 2017,
### https://doi.org/10.32614/RJ-2017-036) as the "Land of Oz weather Matrix".
###
### States: Rainy, Nice, Snowy. The chain is regular (P^6 already has all
### positive entries, converging to the stationary row (.4, .2, .4) -- see
### Grinstead & Snell, Table 11.1), hence irreducible and aperiodic. This
### makes it a good public reference case for isRegular(), which currently
### has no test coverage anywhere in the package's test suite.

ozMatr <- matrix(c(1/2, 1/4, 1/4,
                    1/2,  0 , 1/2,
                    1/4, 1/4, 1/2),
                  nrow = 3, byrow = TRUE,
                  dimnames = list(c("R", "N", "S"), c("R", "N", "S")))
mcOz <- new("markovchain", states = c("R", "N", "S"), byrow = TRUE,
            transitionMatrix = ozMatr, name = "Land of Oz")

test_that("Land of Oz weather chain is regular, irreducible and aperiodic", {
  expect_true(isRegular(mcOz))
  expect_true(is.irreducible(mcOz))
  expect_equal(period(mcOz), 1)
  expect_equal(communicatingClasses(mcOz), list(c("R", "N", "S")))
  expect_equal(recurrentClasses(mcOz), list(c("R", "N", "S")))
  expect_equal(length(absorbingStates(mcOz)), 0)
  expect_equal(length(transientStates(mcOz)), 0)
})

test_that("Land of Oz stationary distribution matches the textbook value", {
  ## Grinstead & Snell, Section 11.3: w = (.4, .2, .4)
  expect_equal(unname(steadyStates(mcOz)[1, ]), c(0.4, 0.2, 0.4), tolerance = 1e-6)
})

### A simple periodic (2-cycle) irreducible chain: irreducible but NOT
### regular, since it never mixes -- a standard textbook counterexample
### to "irreducible implies regular" (e.g. Grinstead & Snell, Section 11.3,
### remark that regularity is strictly stronger than irreducibility).
### Useful as a negative control alongside the positive one above.

cycMatr <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE,
                   dimnames = list(c("a", "b"), c("a", "b")))
mcCycle <- new("markovchain", states = c("a", "b"), byrow = TRUE,
               transitionMatrix = cycMatr, name = "2-cycle")

test_that("An irreducible periodic chain is not regular", {
  expect_true(is.irreducible(mcCycle))
  expect_false(isRegular(mcCycle))
  expect_equal(period(mcCycle), 2)
})

### mathematicaMc (already used elsewhere in the test suite, sourced from
### Wolfram Mathematica's MarkovProcessProperties documentation examples,
### and reused in the package vignette). It is reducible, so it must not be
### regular; this strengthens the existing coverage, which so far checked
### recurrentClasses()/transientStates() but never isRegular() or
### communicatingClasses() directly, and only checked canonicForm() for
### self-consistency against the internal Rcpp call rather than against a
### known-correct reordering.

mathematicaMatr <- markovchain:::zeros(5)
mathematicaMatr[1, ] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2, ] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3, ] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4, ] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5, ] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
                     name = "Mathematica MC", states = statesNames)

test_that("mathematicaMc is reducible and therefore not regular", {
  expect_false(isRegular(mathematicaMc))
  expect_false(is.irreducible(mathematicaMc))
  expect_equal(
    communicatingClasses(mathematicaMc),
    list(c("a", "b"), c("c", "d"), c("e"))
  )
})

test_that("canonicForm places recurrent classes first, in block form, transient states last", {
  cf <- canonicForm(mathematicaMc)
  ## recurrent classes {c,d} and {e} first (in the order recurrentClasses()
  ## reports them), transient states a, b last, in their original order
  expect_equal(cf@states, c("c", "d", "e", "a", "b"))
  ## the recurrent block must be block-diagonal: no probability mass leaks
  ## from one recurrent class into a different one
  expect_equal(unname(cf@transitionMatrix["c", "e"]), 0)
  expect_equal(unname(cf@transitionMatrix["d", "e"]), 0)
  expect_equal(unname(cf@transitionMatrix["e", "c"]), 0)
  expect_equal(unname(cf@transitionMatrix["e", "d"]), 0)
})
