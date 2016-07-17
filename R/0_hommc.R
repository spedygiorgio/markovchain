#' An S4 class for representing High Order Multivariate Markovchain (HOMMC)
#' 
#' @slot order an integer equal to order of Multivariate Markovchain
#' @slot states a vector of states present in the HOMMC model
#' @slot P array of transition matrices
#' @slot Lambda a vector which stores the weightage of each transition matrices in P
#' @slot byrow if FALSE each column sum of transition matrix is 1 else row sum = 1
#' @slot name a name given to hommc
#' 
#' @author Giorgio Spedicato, Deepak Yadav
#' 
#' @examples 
#' statesName <- c("a", "b")
#' 
#' P <- array(0, dim = c(2, 2, 4), dimnames = list(statesName, statesName))
#' P[,,1] <- matrix(c(0, 1, 1/3, 2/3), byrow = FALSE, nrow = 2)
#' P[,,2] <- matrix(c(1/4, 3/4, 0, 1), byrow = FALSE, nrow = 2)
#' P[,,3] <- matrix(c(1, 0, 1/3, 2/3), byrow = FALSE, nrow = 2)
#' P[,,4] <- matrix(c(3/4, 1/4, 0, 1), byrow = FALSE, nrow = 2)
#' 
#' Lambda <- c(0.8, 0.2, 0.3, 0.7)
#' 
#' ob <- new("hommc", order = 1, states = statesName, P = P, 
#'           Lambda = Lambda, byrow = FALSE, name = "FOMMC")

hommc <- setClass("hommc",
                      slots = list(order = "numeric", states  =  "character",
                      P = "array", Lambda = "numeric", byrow = "logical",
                      name = "character")
)