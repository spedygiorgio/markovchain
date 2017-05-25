require(MCMCpack)
require(markovchain)
#create the function to verify the inversion
checkInversion<-function(i,num){
  #simulate the priors
  priors.dirichlet<-runif(n=num)
  myStochasticMatrix<-rdirichlet(n=num,alpha=priors.dirichlet)
  #this code returns FALSE -> 0 if error in inversion 1 otherwise
  out<-tryCatch(steadyStates(as(myStochasticMatrix, "markovchain")),
                error=function(c) return(FALSE)
  )
  if(class(out)=="logical") return(0) else return(1)
}

dimensions2Test<-2:32
successRate<-numeric(length(dimensions2Test))


library(parallel)
no_cores <- detectCores() - 1
numSim=10000
cl <- makeCluster(no_cores)
clusterExport(cl, "checkInversion")
clusterEvalQ(cl, library(markovchain))
clusterEvalQ(cl, library(MCMCpack))
k=1
for (dimension in dimensions2Test){
  simulations<-parSapply(cl=cl,1:numSim,FUN=checkInversion,num=dimension)
  errorRate[k]<-mean(simulations)
  k=k+1
}
stopCluster(cl)

plot(x=dimensions2Test,y=errorRate,type="l");abline(h=0.5,col="red")