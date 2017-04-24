#' @name sales
#' 
#' @title Sales Demand Sequences
#' 
#' @description Sales demand sequences of five products (A, B, C, D, E).
#'              Each row corresponds to a sequence. First row corresponds to Sequence A, 
#'              Second row to Sequence B and so on.
#' 
#' @usage data("sales")
#' 
#' @details The example can be used to fit High order multivariate
#'          markov chain.
#' 
#' @examples 
#' data("sales")
#' # fitHighOrderMultivarMC(seqMat = sales, order = 2, Norm = 2)
#' 
"sales"