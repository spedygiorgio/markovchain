#' @docType package
#' @name markovchain-package
#' @rdname markovchain
#' @title Easy Handling Discrete Time Markov Chains
#'
#' @description The package contains classes and method to create and manage
#'   (plot, print, export for example) discrete time Markov chains (DTMC). In
#'   addition it provide functions to perform statistical (fitting and drawing
#'   random variates) and probabilistic (analysis of DTMC proprieties) analysis
#'   
#' @author 
#' Giorgio Alfredo Spedicato 
#' Maintainer: Giorgio Alfredo Spedicato <spedicato_giorgio@yahoo.it>
#' @references Discrete-Time Markov Models, Bremaud, Springer 1999
#' @keywords package
#' 
#' @examples
#' # create some markov chains
#' statesNames=c("a","b")
#' mcA<-new("markovchain", transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),byrow=TRUE,
#'          nrow=2, dimnames=list(statesNames,statesNames)))
#'          
#' statesNames=c("a","b","c")
#' mcB<-new("markovchain", states=statesNames, transitionMatrix=
#'          matrix(c(0.2,0.5,0.3,0,1,0,0.1,0.8,0.1), nrow=3, 
#'          byrow=TRUE, dimnames=list(statesNames, statesNames)))
#'
#' statesNames=c("a","b","c","d")
#' matrice<-matrix(c(0.25,0.75,0,0,0.4,0.6,0,0,0,0,0.1,0.9,0,0,0.7,0.3), nrow=4, byrow=TRUE)
#' mcC<-new("markovchain", states=statesNames, transitionMatrix=matrice)
#' mcD<-new("markovchain", transitionMatrix=matrix(c(0,1,0,1), nrow=2,byrow=TRUE))
#'
#'
#' #operations with S4 methods
#' mcA^2
#' steadyStates(mcB)
#' absorbingStates(mcB)
#' markovchainSequence(n=20, markovchain=mcC, include=TRUE)
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
#' @importFrom expm %^% logm
#' @importFrom stats sd rexp chisq.test pchisq predict aggregate
#' @importFrom grDevices colors
NULL