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
#' @name blanden
#' 
#' @title Mobility between income quartiles
#' 
#' @description This table show mobility between income quartiles for father and sons for the 1970 cohort born
#' 
#' @usage data(blanden)
#' 
#' @details The rows represent fathers' income quartile when the son is aged 16, whilst the columns represent sons' income quartiles when he is aged 30 (in 2000).
#' 
#' @source Personal reworking
#' 
#' @references Jo Blanden, Paul Gregg and Stephen Machin, Intergenerational Mobility in Europe and North America, Center for Economic Performances (2005)
#' 
#' @example 
#' data(blanden)
#' mobilityMc<-as(blanden, "markovchain")
#' 
"blanden"
