# define higher order Markov Chain class

setClass("HigherOrderMarkovChain", #class name
         representation(
           states = "character", 
           order = "numeric",
           transitions = "list", 
           name = "character"
         )
#          , prototype(states = c("a","b"), byrow = TRUE, # prototypizing
#                    transitionMatrix=matrix(data = c(0,1,1,0),
#                                            nrow=2, byrow=TRUE, dimnames=list(c("a","b"), c("a","b"))),
#                    name="Unnamed Markov chain")
)

# setClass(
#   "HigerOrderMarkovChain",
#   representation(
#     states = "character",
#     order = "numeric",
#     transitions = "list",
#     lambda = "numeric",
#     logLikelihood = "numeric",
#     observations = "numeric",
#     start = "table",
#     end = "table",
#     transientStates = "character",
#     absorbingStates = "character",
#     absorbingProbabilities = "data.frame"
#   )
# )


