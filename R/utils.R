precomputeData <- function(mc) {
  list(
    object = mc,
    transitionMatrix = mc@transitionMatrix,
    states = mc@states,
    byrow = mc@byrow,
    irreducible = is.irreducible(mc),
    canonicForm = canonicForm(mc),
    recurrentClasses = recurrentClasses(mc),
    transientClasses = transientClasses(mc),
    recurrentStates = recurrentStates(mc),
    transientStates = transientStates(mc),
    absorbingStates = absorbingStates(mc),
    hittingProbabilities = hittingProbabilities(mc),
    communicatingClasses = communicatingClasses(mc),
    steadyStates = steadyStates(mc)
  )
}