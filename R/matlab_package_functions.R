## matlab brougt core functions
##-----------------------------------------------------------------------------
setClass("size_t",
         contains = "integer",
         prototype = as.integer(0))


##-----------------------------------------------------------------------------
size_t <- function(x) {
  new("size_t", as.integer(x))
}


##-----------------------------------------------------------------------------
is.size_t <- function(object) {
  return(data.class(object) == "size_t")
}


##-----------------------------------------------------------------------------
as.size_t <- function(object) {
  return(size_t(object))
}




##-----------------------------------------------------------------------------
find <- function(x) {
  expr <- if (is.logical(x)) {
    x
  } else {
    x != 0
  }
  return(which(expr))
}

### Create an identity matrix.
###


##-----------------------------------------------------------------------------
eye <- function(m, n) {
  if (is.size_t(m)) {
    m <- as.integer(m)
  }
  
  if (missing(n)) {
    len.m <- length(m)
    if (len.m == 1) {
      n <- m
    } else if (len.m > 1) {
      n <- m[-1]
      m <- m[1]
    }
  }
  
  if (!is.numeric(n)) {
    stop(sprintf("argument %s must be numeric", sQuote("n")))
  } else if (!(length(n) == 1)) {
    stop(sprintf("argument %s must be of length 1", sQuote("n")))
  } else if (!(n > 0)) {
    stop(sprintf("argument %s must be a positive quantity", sQuote("n")))
  }
  
  if (!is.numeric(m)) {
    stop(sprintf("argument %s must be numeric", sQuote("m")))
  } else if (!(length(m) == 1)) {
    stop(sprintf("argument %s must be of length 1", sQuote("m")))
  } else if (!(m > 0)) {
    stop(sprintf("argument %s must be a positive quantity", sQuote("m")))
  }
  
  return(diag(1, m, n))
}



#' Create a matrix of zeros of size
#'
#' @param ... typically the size of the matrix
#'
#' @return A matrix of zeros
#' @export
#'
#' @examples
#' zeros(c(2,2))
zeros <- function(...) {
  nargs <- length(dots <- list(...))
  dims <- as.integer(if (nargs == 1 && is.size_t(dots[[1]])) {
    dots[[1]]
  } else {
    unlist(dots)
  })
  
  if (length(dims) == 1) {
    dims[2] <- dims[1]
  }
  
  if (!(length(dims) > 1)) {
    stop("dimensions must be of length greater than 1")
  } else if (!(all(dims > 0))) {
    stop("dimensions must be a positive quantity")
  }
  
  return(array(0, dims))
}

##-----------------------------------------------------------------------------
ones <- function(...) {
  nargs <- length(dots <- list(...))
  dims <- as.integer(if (nargs == 1 && is.size_t(dots[[1]])) {
    dots[[1]]
  } else {
    unlist(dots)
  })
  
  if (length(dims) == 1) {
    dims[2] <- dims[1]
  }
  
  if (!(length(dims) > 1)) {
    stop("dimensions must be of length greater than 1")
  } else if (!(all(dims > 0))) {
    stop("dimensions must be a positive quantity")
  }
  
  return(array(1, dims))
}

###
### $Id: size.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Array dimensions.
###


##-----------------------------------------------------------------------------
setGeneric("size",
           function(X, dimen) {
             #cat("generic", match.call()[[1]], "\n")
             standardGeneric("size")
           })

setMethod("size",
          signature(X     = "vector",
                    dimen = "missing"),
          function(X, dimen) {
            #cat(match.call()[[1]], "(vector, missing)", "\n")
            #              return(as.size_t(length(X)))
            # :NOTE: Incompatible with previous implementation but consistent with MATLAB
            callGeneric(matrix(X, nrow = 1))
          })

setMethod("size",
          signature(X     = "matrix",
                    dimen = "missing"),
          function(X, dimen) {
            #cat(match.call()[[1]], "(matrix, missing)", "\n")
            return(as.size_t(dim(X)))
          })

setMethod("size",
          signature(X     = "array",
                    dimen = "missing"),
          function(X, dimen) {
            #cat(match.call()[[1]], "(array, missing)", "\n")
            return(as.size_t(dim(X)))
          })

setMethod("size",
          signature(X     = "vector",
                    dimen = "numeric"),
          function(X, dimen) {
            #cat(match.call()[[1]], "(vector, numeric)", "\n")
            callGeneric(matrix(X, nrow = 1), dimen)
          })

setMethod("size",
          signature(X     = "matrix",
                    dimen = "numeric"),
          function(X, dimen) {
            #cat(match.call()[[1]], "(matrix, numeric)", "\n")
            callGeneric(X, as.integer(dimen))
          })

setMethod("size",
          signature(X     = "matrix",
                    dimen = "integer"),
          function(X, dimen) {
            #cat(match.call()[[1]], "(matrix, integer)", "\n")
            return(getLengthOfDimension(X, dimen))
          })

setMethod("size",
          signature(X     = "array",
                    dimen = "numeric"),
          function(X, dimen) {
            #cat(match.call()[[1]], "(array, numeric)", "\n")
            callGeneric(X, as.integer(dimen))
          })

setMethod("size",
          signature(X     = "array",
                    dimen = "integer"),
          function(X, dimen) {
            #cat(match.call()[[1]],
            #    "(", data.class(X), ", ", data.class(dimen), ")", "\n")
            return(getLengthOfDimension(X, dimen))
          })

setMethod("size",
          signature(X     = "missing",
                    dimen = "ANY"),
          function(X, dimen) {
            #cat(match.call()[[1]], "(missing, ANY)", "\n")
            stop(sprintf("argument %s missing", sQuote("X")))
          })

##-----------------------------------------------------------------------------
getLengthOfDimension <- function(X, dimen) {
  if (!is.array(X)) {
    stop(sprintf("argument %s must be matrix or array", sQuote("X")))
  }
  
  if (!(length(dimen) == 1)) {
    stop(sprintf("argument %s must be of length 1", sQuote("dimen")))
  } else if (dimen < 1) {
    stop(sprintf("argument %s must be a positive quantity",
                 sQuote("dimen")))
  }
  
  len <- if (dimen <= length(dim(X))) {
    dim(X)[dimen]
  } else {
    1	# singleton dimension
  }
  
  return(as.integer(len))
}

