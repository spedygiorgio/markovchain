# Monte Carlo validation study for statistical tests in markovchain
#
# Run from the package root, for example:
#   Rscript inst/validation/validate_statistical_tests.R
#
# The study is deliberately outside R CMD check. It uses the base/recommended
# parallel package and a local PSOCK cluster, so it works on Windows as well as
# Unix-alike systems. By default it uses all detected logical CPUs minus one.
# Override with CORES, NREP, MC_B and SEED environment variables.

library(markovchain)
library(parallel)

nrep <- as.integer(Sys.getenv("NREP", "500"))
mc_B <- as.integer(Sys.getenv("MC_B", "199"))
seed <- as.integer(Sys.getenv("SEED", "20260815"))
requested_cores <- as.integer(Sys.getenv("CORES", "0"))

alpha <- c(0.01, 0.05, 0.10)
lengths <- c(50L, 100L, 250L, 500L, 1000L)
states_n <- c(2L, 3L, 5L)

random_transition_matrix <- function(r, concentration = 2) {
  P <- matrix(0, r, r)
  for (i in seq_len(r)) {
    P[i, ] <- rgamma(r, shape = concentration, rate = 1)
    P[i, ] <- P[i, ] / sum(P[i, ])
  }
  dimnames(P) <- list(as.character(seq_len(r)), as.character(seq_len(r)))
  P
}

make_mc <- function(P) {
  new("markovchain", states = rownames(P), transitionMatrix = P)
}

second_order_sequence <- function(n, P, delta = 0.35, t0 = NULL) {
  r <- nrow(P)
  if (is.null(t0)) t0 <- sample.int(r, 1L)
  x <- integer(n)
  x[1] <- t0
  x[2] <- sample.int(r, 1L, prob = P[x[1], ])
  for (t in 2:(n - 1L)) {
    target <- rep(0, r)
    target[x[t - 1L]] <- 1
    q <- (1 - delta) * P[x[t], ] + delta * target
    x[t + 1L] <- sample.int(r, 1L, prob = q)
  }
  as.character(x)
}

run_markov_property <- function(n, r, nrep, mc_B, task_seed) {
  set.seed(task_seed)
  P <- random_transition_matrix(r)
  mc <- make_mc(P)
  null_p_G <- null_p_P <- null_p_MC <- numeric(nrep)
  alt_p_G <- alt_p_P <- numeric(nrep)

  for (b in seq_len(nrep)) {
    x0 <- rmarkovchain(n, mc, t0 = rownames(P)[1L])
    null_p_G[b] <- verifyMarkovProperty(x0, method = "G", verbose = FALSE)$p.value
    null_p_P[b] <- verifyMarkovProperty(x0, method = "Pearson", verbose = FALSE)$p.value
    null_p_MC[b] <- verifyMarkovProperty(x0, method = "simulation",
                                          B = mc_B, seed = b,
                                          verbose = FALSE)$p.value

    x1 <- second_order_sequence(n, P)
    alt_p_G[b] <- verifyMarkovProperty(x1, method = "G", verbose = FALSE)$p.value
    alt_p_P[b] <- verifyMarkovProperty(x1, method = "Pearson", verbose = FALSE)$p.value
  }

  do.call(rbind, lapply(alpha, function(a) {
    data.frame(test = c("G", "Pearson", "MonteCarlo", "G", "Pearson"),
               scenario = c("null", "null", "null", "second-order", "second-order"),
               n = n, states = r, alpha = a,
               rejection_rate = c(mean(null_p_G < a), mean(null_p_P < a),
                                  mean(null_p_MC < a), mean(alt_p_G < a),
                                  mean(alt_p_P < a)))
  }))
}

run_empirical_theoretical <- function(n_per_state = 100, r = 3, nrep = 500,
                                      task_seed) {
  set.seed(task_seed)
  P <- random_transition_matrix(r)
  mc <- make_mc(P)
  p_G <- p_P <- p_MC <- numeric(nrep)
  for (b in seq_len(nrep)) {
    counts <- matrix(0, r, r, dimnames = dimnames(P))
    for (i in seq_len(r)) counts[i, ] <- rmultinom(1, n_per_state, P[i, ])
    p_G[b] <- verifyEmpiricalToTheoretical(counts, mc, "G", FALSE)$p.value
    p_P[b] <- verifyEmpiricalToTheoretical(counts, mc, "Pearson", FALSE)$p.value
    p_MC[b] <- verifyEmpiricalToTheoretical(counts, mc, "simulation",
                                             B = mc_B, seed = b,
                                             verbose = FALSE)$p.value
  }
  do.call(rbind, lapply(alpha, function(a)
    data.frame(test = c("G", "Pearson", "MonteCarlo"), scenario = "null",
               n = n_per_state, states = r, alpha = a,
               rejection_rate = c(mean(p_G < a), mean(p_P < a), mean(p_MC < a)))))
}

run_homogeneity <- function(n = 500, r = 3, groups = 3, nrep = 500,
                            task_seed) {
  set.seed(task_seed)
  P <- random_transition_matrix(r)
  mc <- make_mc(P)
  p_G <- p_P <- p_MC <- numeric(nrep)
  P_alt <- P
  P_alt[1, ] <- c(0.80, rep(0.20 / (r - 1), r - 1))
  if (r == 2) P_alt[1, ] <- c(0.80, 0.20)
  alt_G <- alt_P <- alt_MC <- numeric(nrep)

  for (b in seq_len(nrep)) {
    null_data <- lapply(seq_len(groups), function(g)
      rmarkovchain(n, mc, t0 = rownames(P)[1L]))
    p_G[b] <- verifyHomogeneity(null_data, "G", FALSE)$p.value
    p_P[b] <- verifyHomogeneity(null_data, "Pearson", FALSE)$p.value
    p_MC[b] <- verifyHomogeneity(null_data, "simulation", B = mc_B,
                                  seed = b, verbose = FALSE)$p.value

    alt_data <- lapply(seq_len(groups), function(g) {
      P_use <- if (g == 2L) P_alt else P
      rmarkovchain(n, make_mc(P_use), t0 = rownames(P)[1L])
    })
    alt_G[b] <- verifyHomogeneity(alt_data, "G", FALSE)$p.value
    alt_P[b] <- verifyHomogeneity(alt_data, "Pearson", FALSE)$p.value
    alt_MC[b] <- verifyHomogeneity(alt_data, "simulation", B = mc_B,
                                    seed = b, verbose = FALSE)$p.value
  }

  do.call(rbind, lapply(alpha, function(a)
    data.frame(test = c("G", "Pearson", "MonteCarlo", "G", "Pearson", "MonteCarlo"),
               scenario = c("null", "null", "null", "heterogeneous", "heterogeneous", "heterogeneous"),
               n = n, states = r, alpha = a,
               rejection_rate = c(mean(p_G < a), mean(p_P < a), mean(p_MC < a),
                                  mean(alt_G < a), mean(alt_P < a), mean(alt_MC < a)))))
}

# Build independent scenario tasks. Each task has its own deterministic seed,
# so parallel load balancing does not affect reproducibility.
tasks <- list()
task_id <- 0L
for (r in states_n) {
  for (n in lengths) {
    task_id <- task_id + 1L
    tasks[[task_id]] <- list(type = "markov", n = n, r = r,
                             seed = seed + task_id)
  }
}
task_id <- task_id + 1L
tasks[[task_id]] <- list(type = "empirical", n = 100L, r = 3L,
                        seed = seed + task_id)
task_id <- task_id + 1L
tasks[[task_id]] <- list(type = "homogeneity", n = 500L, r = 3L,
                        seed = seed + task_id)

run_task <- function(task) {
  message("Running ", task$type, " scenario: n=", task$n,
          ", r=", task$r)
  if (task$type == "markov")
    return(run_markov_property(task$n, task$r, nrep, mc_B, task$seed))
  if (task$type == "empirical")
    return(run_empirical_theoretical(task$n, task$r, nrep, task$seed))
  run_homogeneity(task$n, task$r, nrep = nrep, task_seed = task$seed)
}

available <- detectCores(logical = TRUE)
if (is.na(available)) available <- 1L
default_cores <- max(1L, available - 1L)
cores <- if (requested_cores > 0L) requested_cores else default_cores
cores <- max(1L, min(cores, length(tasks)))

message("Detected logical CPUs: ", available)
message("Validation workers: ", cores)
message("Outer replications per scenario: ", nrep)
message("Monte Carlo replicates per test: ", mc_B)

cl <- makeCluster(cores, type = "PSOCK", outfile = "")
on.exit(stopCluster(cl), add = TRUE)
clusterSetRNGStream(cl, iseed = seed)
clusterEvalQ(cl, library(markovchain))

# Load balancing is useful because the 15 Markov-property scenarios have
# different computational costs. Reproducibility is preserved by task seeds.
results <- parLapplyLB(cl, tasks, run_task)
results <- do.call(rbind, results)

out_dir <- file.path("inst", "validation", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(out_dir, "statistical_tests_validation.csv"), row.names = FALSE)
saveRDS(results, file.path(out_dir, "statistical_tests_validation.rds"))

print(results)
message("Validation completed. Results written to ", out_dir)
