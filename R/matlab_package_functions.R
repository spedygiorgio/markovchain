#' Matrix to create zeros
#'
#' @param n size of the matrix
#'
#' @return a square matrix of zeros
# Small helpers named after their MATLAB counterparts (zeros(), ones()).
# Both build an n x n matrix filled with a constant.
zeros <- function(n) {
  # matrix(value, nrow, ncol): here an n x n matrix of zeros
  return(matrix(0,n,n))

}

#' Returns an Identity matrix
#'
#' @param n size of the matrix
#'
#' @return a identity matrix 
ones <- function(n) {
  # n x n matrix of ones (despite the roxygen title above, which says 'identity',
  # this is NOT the identity matrix: use diag(n) for that)
  return(matrix(1,n,n))
}
