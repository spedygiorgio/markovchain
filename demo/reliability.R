#load required libraries
library(parallel)
require(MCMCpack)
require(markovchain)
dimensions2Test<-2:32
numSim=10000 
 
#helper function to create a random stochastic matrix 
 
createMatrix<-function(matr_size) { 
  out<-matrix(0, nrow=matr_size, ncol = matr_size) 
  for (i in 1:matr_size) { 
    priors.dirichlet<-runif(n=matr_size) 
    myStochasticRow<-rdirichlet(n=1,alpha=priors.dirichlet) 
    out[i,]<-myStochasticRow 
  } 
  return(out) 
}


createSparseMatrix<-function(matr_size, sparsity=0.75) {
  out <- matrix(0, nrow=matr_size, ncol = matr_size)
  nonzeroitems<-ceiling(matr_size*(1-sparsity))
  for (i in 1:matr_size) {
    priors.dirichlet<-runif(n=nonzeroitems)
    myStochasticRow<-rdirichlet(n=1,alpha=priors.dirichlet)
    columnsPositions<-sample(x = 1:matr_size,size = nonzeroitems,replace = FALSE)
    out[i,columnsPositions]<-myStochasticRow
  }
  return(out)
}


#first test: function to simulate the inversion of a matrix of size num
checkInversion<-function(i,num){
  #simulate the priors
  myStochasticMatrix<-createMatrix(matr_size = num)
  #this code returns FALSE -> 0 if error in inversion 1 otherwise
  out<-tryCatch(steadyStates(as(myStochasticMatrix, "markovchain")),
                error=function(c) return(FALSE)
  )
  if(class(out)=="logical") return(0) else return(1)
}

checkSparseMInversion<-function(i,num){
  #simulate the priors
  myStochasticMatrix<-createSparseMatrix(matr_size = num)
  #this code returns FALSE -> 0 if error in inversion 1 otherwise
  out<-tryCatch(steadyStates(as(myStochasticMatrix, "markovchain")),
                error=function(c) return(FALSE)
  )
  if(class(out)=="logical") return(0) else return(1)
}




#performing the simulation

successRate<-numeric(length(dimensions2Test))

#using parallel backend
no_cores <- detectCores() - 1 
cl <- makeCluster(no_cores)
clusterExport(cl, "checkInversion")
clusterExport(cl, "createMatrix")
clusterEvalQ(cl, library(markovchain))
clusterEvalQ(cl, library(MCMCpack))
k=1
for (dimension in dimensions2Test){
  simulations<-parSapply(cl=cl,1:numSim,FUN=checkInversion,num=dimension)
  successRate[k]<-mean(simulations)
  k=k+1
}
stopCluster(cl)
#summarising first test:
#plot(x=dimensions2Test,y=successRate,type="l",xlab="matrix size",ylab="success rate",main="Steady state computation success rate")
#abline(h=0.5,col="red")
#text(x=dimensions2Test,y=successRate,labels=round(successRate,digits=2),col="darkred",cex=0.7)
#dev.off()

dimensions2Test = 2^seq(from=3, to=8)
successRate<-numeric(length(dimensions2Test))

#using parallel backend
no_cores <- detectCores() - 1 
cl <- makeCluster(no_cores)
clusterExport(cl, "checkSparseMInversion")
clusterExport(cl, "createSparseMatrix")
clusterEvalQ(cl, library(markovchain))
clusterEvalQ(cl, library(MCMCpack))
k=1
for (dimension in dimensions2Test){
  simulations<-parSapply(cl=cl,1:numSim,FUN=checkSparseMInversion,num=dimension)
  successRate[k]<-mean(simulations)
  k=k+1
}
stopCluster(cl)
# 
# 
# plot(x=dimensions2Test,y=successRate,type="l",xlab="matrix size",ylab="success rate",main="Steady state computation success rate ??? sparse matrices")
# abline(h=0.5,col="red")
# text(x=dimensions2Test,y=successRate,labels=round(successRate,digits=2),col="darkred",cex=0.7)



#second test: simulating exponentiation

checkExponentiation<-function(i,num){
  #simulate the priors
  myStochasticMatrix<-createMatrix(matr_size = num)
  #this code returns FALSE -> 0 if error in inversion 1 otherwise
  out<-tryCatch((as(myStochasticMatrix, "markovchain"))^2,
                error=function(c) return(FALSE)
  )
  if(class(out)=="logical") return(0) else return(1)
}


#performing the simulation
successRate2<-numeric(length(dimensions2Test))


#using parallel backend
no_cores <- detectCores() - 1 
cl <- makeCluster(no_cores)
clusterExport(cl, "checkExponentiation")
clusterExport(cl, "createMatrix")
clusterEvalQ(cl, library(markovchain))
clusterEvalQ(cl, library(MCMCpack))
k=1
for (dimension in dimensions2Test){
  simulations<-parSapply(cl=cl,1:numSim,FUN=checkExponentiation,num=dimension)
  successRate2[k]<-mean(simulations)
  k=k+1
}
stopCluster(cl)
#summarising first test:


#par(mfrow=c(1,2))

# plot(x=dimensions2Test,y=successRate,type="l",xlab="matrix size",ylab="success rate",main="Steady state computation success rate")
# abline(h=0.5,col="red")
# text(x=dimensions2Test,y=successRate,labels=round(successRate,digits=2),col="darkred",cex=0.7)



# plot(x=dimensions2Test,y=successRate2,type="l",xlab="matrix sixe",ylab="success rate",main="Exponentiation computation success rate")
# abline(h=0.5,col="red")
# text(x=dimensions2Test,y=successRate2,labels=round(successRate2,digits=2),col="darkred",cex=0.7)

