# Fills in missing row/column names of a square matrix and resolves the
# corresponding states vector. Shared by the initialize methods of
# markovchain, ctmc and ictmc, which all followed this same convention
# (previously duplicated almost verbatim in each of them):
#   - if both row and column names are missing, states are taken from
#     `states` (or generated as "1", "2", ... when that is also NULL)
#   - if only one side is missing, it is copied from the other side
#   - if the two sides disagree, row names win
#
# Args:
#   matr: a square matrix (transition matrix or generator)
#   states: character vector of state names, or NULL if not supplied
#
# Returns:
#   a list with `matrix` (dimnames filled in) and `states`
.fillDimNames <- function(matr, states = NULL) {
  rowNames <- rownames(matr)
  colNames <- colnames(matr)

  if (is.null(rowNames) && is.null(colNames)) {
    stateNames <- if (is.null(states)) as.character(seq_len(nrow(matr))) else states
    rownames(matr) <- stateNames
    colnames(matr) <- stateNames
  } else if (is.null(rowNames)) {
    rownames(matr) <- colNames
  } else if (is.null(colNames)) {
    colnames(matr) <- rowNames
  } else if (!setequal(rowNames, colNames)) {
    colnames(matr) <- rowNames
  }

  if (is.null(states)) states <- rownames(matr)

  list(matrix = matr, states = states)
}

precomputeData <- function(mc) {
  list(
    object = mc,
    transitionMatrix = mc@transitionMatrix,
    states = mc@states,
    byrow = mc@byrow,
    irreducible = is.irreducible(mc),
    regular = is.regular(mc),
    canonicForm = canonicForm(mc),
    recurrentClasses = recurrentClasses(mc),
    transientClasses = transientClasses(mc),
    recurrentStates = recurrentStates(mc),
    transientStates = transientStates(mc),
    absorbingStates = absorbingStates(mc),
    hittingProbabilities = hittingProbabilities(mc),
    meanNumVisits = meanNumVisits(mc),
    meanRecurrenceTime = meanRecurrenceTime(mc),
    communicatingClasses = communicatingClasses(mc),
    steadyStates = steadyStates(mc),
    reachabilityMatrix = is.accessible(mc)
  )
}

precomputeSteadyStates <- function(mc) {
  list(
    object = mc,
    expected = steadyStates(mc)
  )
}