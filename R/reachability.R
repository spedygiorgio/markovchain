# Normalize the matrix orientation before calling the C++ reachability kernel.
# The kernel interprets the matrix as row-stochastic, while a markovchain
# object may represent transitions by columns when byrow = FALSE.
setMethod("is.accessible", c("markovchain", "missing", "missing"),
  function(object, from, to) {
    if (!object@byrow) {
      object <- object
      object@transitionMatrix <- t(object@transitionMatrix)
      object@byrow <- TRUE
    }

    .reachabilityMatrixRcpp(object)
  }
)
