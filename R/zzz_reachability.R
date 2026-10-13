# Normalize the matrix orientation before calling the C++ reachability kernel.
# The kernel interprets the matrix as row-stochastic, while a markovchain
# object may represent transitions by columns when byrow = FALSE.
# NOTE for beginners: this file is named zzz_* on purpose. R loads package files in
# alphabetical order, and setMethod() needs the generic 'is.accessible' (defined in
# another file) to already exist, so this file must be loaded after it.
setMethod("is.accessible", c("markovchain", "missing", "missing"),
  function(object, from, to) {
    # byrow = FALSE means the matrix stores probabilities by column (column-stochastic).
    # We transpose it so that the C++ code can always assume rows sum to one.
    if (!object@byrow) {
      object@transitionMatrix <- t(object@transitionMatrix)
      object@byrow <- TRUE
    }

    # Call the compiled routine (src/probabilistic.cpp) that returns the matrix of
    # who-can-reach-whom (TRUE/FALSE) between all pairs of states.
    .reachabilityMatrixRcpp(object)
  }
)
