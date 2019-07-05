#' markovchain: A package for Markov chains
#'
#' Package to represent DTMC and CTMC Markov chains, extract information from
#' them, related to their states, their stationarity and do inference
#'
#' @docType package
#' @name markovchain
NULL

#' @useDynLib markovchain, .registration = TRUE
#' @import igraph
#' @import Matrix
#' @import methods
#' @import parallel
#' @importFrom utils packageDescription
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats4 plot summary
#' @importFrom matlab zeros find eye size ones
#' @importFrom expm %^% logm
#' @importFrom stats sd rexp chisq.test pchisq predict aggregate
#' @importFrom grDevices colors
NULL