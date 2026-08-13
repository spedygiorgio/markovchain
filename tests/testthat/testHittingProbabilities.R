context("hittingProbabilities")

### Reference example taken directly from the roxygen docs of hittingProbabilities()
### (see probabilistic.R): a 5-state chain absorbing at "a" and "e".
### Exact fractions below were derived analytically via first-step analysis
### and cross-checked against the package's own output.

M <- markovchain:::zeros(5)
M[1, 1] <- M[5, 5] <- 1
M[2, 1] <- M[2, 3] <- 1 / 2
M[3, 2] <- M[3, 4] <- 1 / 2
M[4, 2] <- M[4, 5] <- 1 / 2
statesNames <- letters[1:5]
mcDoc <- new("markovchain", transitionMatrix = M, states = statesNames)

expectedDoc <- matrix(c(
  1,     0,     0,     0,     0,
  0.8,   0.375, 0.5,   1 / 3, 0.2,
  0.6,   0.75,  0.375, 2 / 3, 0.4,
  0.4,   0.5,   0.25,  1 / 6, 0.6,
  0,     0,     0,     0,     1
), nrow = 5, byrow = TRUE, dimnames = list(statesNames, statesNames))

test_that("hittingProbabilities matches the documented example exactly", {
  hp <- hittingProbabilities(mcDoc)
  expect_equal(hp, expectedDoc, tolerance = 1e-8)
})

test_that("hittingProbabilities satisfies its defining recurrence on a regular chain", {
  ## f(i, j) = p(i, j) + sum_{k != j} p(i, k) f(k, j), for all i, j
  ## (this is the same invariant markovchain:::.testthatAreHittingRcpp checks;
  ## calling it directly here avoids re-implementing the recurrence in R).
  hp <- hittingProbabilities(mcDoc)
  expect_true(markovchain:::.testthatAreHittingRcpp(M, hp, TRUE))
})

### Chain with two disjoint closed communicating classes (reuses the structure
### already exercised by "States are those that should be" in
### testStatesClassification.R, kept independent here to isolate failures).

mcMatrClosed <- matrix(c(
  0, 1,   0,   0,   0,   0,
  0.4, 0.6, 0, 0,   0,   0,
  0.3, 0,   0.4, 0.2, 0.1, 0,
  0, 0,   0,   0.3, 0.7, 0,
  0, 0,   0,   0.5, 0,   0.5,
  0, 0,   0,   0.3, 0,   0.7
), nrow = 6, byrow = TRUE)
mcClosed <- as(mcMatrClosed, "markovchain")

test_that("hittingProbabilities gives probability 1 within a closed class and 0 across disjoint closed classes", {
  hp <- hittingProbabilities(mcClosed)
  expect_true(markovchain:::.testthatRecurrentHittingRcpp(
    recurrentClasses(mcClosed), hp, states(mcClosed), TRUE))
})

### Regression test for the numerically ill-conditioned matrix reported as a
### bug: transition probabilities spanning a very wide dynamic range
### (down to ~1e-134) combined with several absorbing states used to make
### the per-target linear solve numerically singular, causing
### hittingProbabilities() to silently return values outside [0, 1] (and,
### more importantly, values that could be *wrong* even after being clipped
### back into [0, 1] -- clipping a negative value to 0 hides the symptom
### without fixing the underlying inaccuracy). The check below verifies
### actual correctness (the defining recurrence), not just that outputs
### happen to land in [0, 1].

states <- c("WT", "A", "B", "D", "E", "F", "A, B", "A, D", "A, E", "A, F", "B, D", "B, E", "B, F", "D, E", "D, F", "E, F", "A, B, D", "A, B, E", "A, B, F", "A, D, E", "A, D, F", "A, E, F", "B, D, E", "B, D, F", "B, E, F", "D, E, F", "A, B, D, E", "A, B, D, F", "A, B, E, F", "A, D, E, F", "B, D, E, F", "A, B, D, E, F")
P <- matrix(0, 32, 32, dimnames = list(states, states))
edges <- rbind(
  c("WT", "WT", 0.63449850680829512), c("WT", "A", 0.041793041492802649), c("WT", "B", 0.028441840284111373), c("WT", "D", 0.062717869820808361), c("WT", "E", 0.21727315706971301), c("WT", "F", 0.015275584524269642),
  c("A", "A", 0.57485088027933162), c("A", "A, B", 0.00041841619555324305), c("A", "A, D", 0.00060531147358965814), c("A", "A, E", 0.42412539205152555), c("A", "A, F", 5.5727123836777542e-55),
  c("B", "B", 0.66647990952108571), c("B", "A, B", 0.00022732644835666204), c("B", "B, D", 0.33284801940052894), c("B", "B, E", 0.00044474463002869947), c("B", "B, F", 6.5589784682323325e-57),
  c("D", "D", 0.51700377454046575), c("D", "A, D", 6.75560039898472e-06), c("D", "B, D", 0.0077652395471692779), c("D", "D, E", 0.47522423031196598), c("D", "D, F", 2.7983195202462127e-24),
  c("E", "E", 0.62040614659464222), c("E", "A, E", 0.0049273518113729601), c("E", "B, E", 7.730104830170793e-06), c("E", "D, E", 0.37465875158986894), c("E", "E, F", 1.9899285658071983e-08),
  c("F", "F", 0.99998362954317199), c("F", "A, F", 2.9757984014636251e-52), c("F", "B, F", 1.961158918919403e-52), c("F", "D, F", 5.0548037329469435e-20), c("F", "E, F", 1.6370456828126428e-05),
  c("A, B", "A, B", 3.0446260799250989e-14), c("A, B", "A, B, D", 0.99999999999949607), c("A, B", "A, B, E", 4.7355752243010144e-13), c("A, B", "A, B, F", 1.0457258192206783e-134),
  c("A, D", "A, D", 0.29890478052849356), c("A, D", "A, B, D", 0.68358435341884782), c("A, D", "A, D, E", 0.017510866052658704), c("A, D", "A, D, F", 1.3750128655650268e-52),
  c("A, E", "A, E", 0.71197625806376208), c("A, E", "A, B, E", 3.3248763467787424e-11), c("A, E", "A, D, E", 0.28802374190298907), c("A, E", "A, E, F", 1.0259561302774829e-45),
  c("A, F", "A, F", 0.21048650109365852), c("A, F", "A, B, F", 3.5306106557747918e-83), c("A, F", "A, D, F", 0.43260041173277847), c("A, F", "A, E, F", 0.3569130871735629),
  c("B, D", "B, D", 0.23612431837855821), c("B, D", "A, B, D", 0.28711333125825389), c("B, D", "B, D, E", 0.47676235036318776), c("B, D", "B, D, F", 9.8237991691823307e-51),
  c("B, E", "B, E", 0.17972411354402271), c("B, E", "A, B, E", 1.9267231831583575e-13), c("B, E", "B, D, E", 0.82027588645578464), c("B, E", "B, E, F", 5.8179252526038295e-54),
  c("B, F", "B, F", 0.21667137765401343), c("B, F", "A, B, F", 5.3936641286941072e-83), c("B, F", "B, D, F", 0.43184914933268037), c("B, F", "B, E, F", 0.35147947301330623),
  c("D, E", "D, E", 0.74472912418700465), c("D, E", "A, D, E", 0.0050770469649198836), c("D, E", "B, D, E", 0.25019382884807545), c("D, E", "D, E, F", 6.8057372942585063e-18),
  c("D, F", "D, F", 0.39768189471835552), c("D, F", "A, D, F", 1.8381568635755127e-33), c("D, F", "B, D, F", 1.1454793777784797e-33), c("D, F", "D, E, F", 0.60231810528164453),
  c("E, F", "E, F", 0.99999999999999889), c("E, F", "A, E, F", 1.8482961686940169e-48), c("E, F", "B, E, F", 1.0223505679033526e-48), c("E, F", "D, E, F", 1.0994420308182125e-15),
  c("A, B, D", "A, B, D", 0.066700029771387268), c("A, B, D", "A, B, D, E", 0.93329997022861266), c("A, B, D", "A, B, D, F", 2.454637913964187e-68),
  c("A, B, E", "A, B, E", 0.030760864616504714), c("A, B, E", "A, B, D, E", 0.96923913538349526), c("A, B, E", "A, B, E, F", 1.7543062778122525e-97),
  c("A, B, F", "A, B, F", 0.33333333333333331), c("A, B, F", "A, B, D, F", 0.33333333333333331), c("A, B, F", "A, B, E, F", 0.33333333333333331),
  c("A, D, E", "A, D, E", 0.031943390019762076), c("A, D, E", "A, B, D, E", 0.96805660998023801), c("A, D, E", "A, D, E, F", 5.0765936912656143e-29),
  c("A, D, F", "A, D, F", 0.39858537171400582), c("A, D, F", "A, B, D, F", 2.3037015754144223e-83), c("A, D, F", "A, D, E, F", 0.60141462828599423),
  c("A, E, F", "A, E, F", 0.26887179941219236), c("A, E, F", "A, B, E, F", 1.3126358255302064e-62), c("A, E, F", "A, D, E, F", 0.73112820058780759),
  c("B, D, E", "B, D, E", 0.11946386921162501), c("B, D, E", "A, B, D, E", 0.880536130788375), c("B, D, E", "B, D, E, F", 2.0009523093922217e-28),
  c("B, D, F", "B, D, F", 0.41684611031332863), c("B, D, F", "A, B, D, F", 1.9316209177149799e-77), c("B, D, F", "B, D, E, F", 0.58315388968667126),
  c("B, E, F", "B, E, F", 0.27902981223905138), c("B, E, F", "A, B, E, F", 3.5693670222903711e-83), c("B, E, F", "B, D, E, F", 0.72097018776094868),
  c("D, E, F", "D, E, F", 1), c("D, E, F", "A, D, E, F", 1.8380432217274186e-31), c("D, E, F", "B, D, E, F", 9.0480304335098159e-32),
  c("A, B, D, E", "A, B, D, E", 1), c("A, B, D, E", "A, B, D, E, F", 1.5449888073552199e-24),
  c("A, B, D, F", "A, B, D, F", 0.5), c("A, B, D, F", "A, B, D, E, F", 0.5),
  c("A, B, E, F", "A, B, E, F", 0.5), c("A, B, E, F", "A, B, D, E, F", 0.5),
  c("A, D, E, F", "A, D, E, F", 1), c("A, D, E, F", "A, B, D, E, F", 4.4492695668057762e-39),
  c("B, D, E, F", "B, D, E, F", 1), c("B, D, E, F", "A, B, D, E, F", 4.3568057771047043e-37),
  c("A, B, D, E, F", "A, B, D, E, F", 1)
)
for (k in seq_len(nrow(edges))) {
  P[edges[k, 1], edges[k, 2]] <- as.numeric(edges[k, 3])
}
P <- P / rowSums(P)
mcIllConditioned <- new("markovchain", transitionMatrix = P)

test_that("hittingProbabilities stays within [0, 1] on an ill-conditioned chain", {
  hp <- hittingProbabilities(mcIllConditioned)
  expect_true(all(hp >= -1e-9))
  expect_true(all(hp <= 1 + 1e-9))
})

test_that("hittingProbabilities satisfies its defining recurrence on an ill-conditioned chain", {
  ## This is the critical regression check: it catches values that are
  ## technically inside [0, 1] but mathematically wrong (e.g. a negative
  ## result silently clipped to 0), which a bounds-only check would miss.
  hp <- hittingProbabilities(mcIllConditioned)
  expect_true(markovchain:::.testthatAreHittingRcpp(P, hp, TRUE))
})

test_that("hittingProbabilities on an ill-conditioned chain agrees with a Monte Carlo estimate", {
  hp <- hittingProbabilities(mcIllConditioned)

  set.seed(1)
  n <- nrow(P)
  Poff <- P
  diag(Poff) <- 0
  rowSumsOff <- rowSums(Poff)
  absorbing <- rowSumsOff <= 1e-12
  J <- Poff
  J[!absorbing, ] <- Poff[!absorbing, ] / rowSumsOff[!absorbing]

  nSim <- 2e4
  visited <- setNames(numeric(n), states)
  startIdx <- match("WT", states)
  for (s in seq_len(nSim)) {
    cur <- startIdx
    visitedThisRun <- logical(n)
    while (!absorbing[cur]) {
      cur <- sample.int(n, 1, prob = J[cur, ])
      visitedThisRun[cur] <- TRUE
    }
    visited <- visited + visitedThisRun
  }
  monteCarloEstimate <- visited / nSim

  ## The starting state's own column is deliberately excluded: hp["WT", "WT"]
  ## is a *return* probability that includes the immediate self-loop
  ## (P(WT, WT) = 0.634...), but the jump chain above strips self-loops by
  ## construction (that's what makes it a jump chain) and so can only ever
  ## measure returns via an indirect path -- which for this chain is exactly
  ## 0, since none of the edges lead back to WT once it's left. Comparing
  ## that entry would be comparing two different definitions, not checking
  ## for a bug; every off-diagonal entry uses the same "ever visited" notion
  ## in both hp and the simulation and is fully comparable.
  off <- states != "WT"
  expect_equal(unname(hp["WT", off]), unname(monteCarloEstimate[off]), tolerance = 0.02)
})
