





#' Matrix to create zeros
#'
#' @param n size of the matrix
#'
#' @return a square matrix of zeros
zeros <- function(n) {
  return(matrix(0,n,n))

}

#' Returns an Identity matrix
#'
#' @param n size of the matrix
#'
#' @return a identity matrix 
ones <- function(n) {
  return(diag(n))
}
