# function to simulate a ctmc
rctmc <- function(n, ctmc, initDist = numeric(), T = 0, include.T0 = TRUE, out.type = "list"){
  if(identical(initDist, numeric()))
    state <- sample(ctmc@states, 1) # sample state randomly
  else if(length(initDist) != dim(ctmc) | round(sum(initDist), 5) != 1)
    stop("Error! Provide a valid initial state probability distribution")
  else 
    state <- sample(ctmc@states, 1, prob = initDist) # if valid probability distribution,
  # sample accordingly
  
  trans <- generatorToTransitionMatrix(ctmc@generator)
  # obtain transition probability matrix from the generator matrix
  
  states <- c()
  time <- c()
  if (include.T0 == TRUE){
    states <- c(states, state)
    time <- c(time, 0)
  }
  
  t <- 0
  i <- 1
  while (i <= n){
    idx <- which(ctmc@states == state)
    t <- t + rexp(1, -ctmc@generator[idx, idx])
    state <- ctmc@states[sample(1:dim(ctmc), 1, prob = trans[idx, ])]
    
    if(T > 0 & t > T)
      break
    
    states <- c(states, state)
    time <- c(time, t)
    i <- i + 1
  }
  
  out <- list(states, time)
  if (out.type == "list")
    return(out)
  else if(out.type == "df"){
    df <- data.frame(matrix(unlist(out), nrow = length(states)))
    names(df) <- c("states", "time")
    return(df)
  }
  else
    stop("Not a valid output type")
}
