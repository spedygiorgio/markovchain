#PSEUDO BAYESIAN ESTIMATOR


mcFitPseudoBayesian<-function(data, apriori) {
  rawTransitionMatrix<-createSequenceMatrix(stringchar = data)
  #obtain the colsums
  vi<-rowSums(rawTransitionMatrix)
  yij2<- rawTransitionMatrix^2
}

apriori=as(matrix(c(0.5,0.5,0.5,0.5),byrow=TRUE,ncol=2),"markovchain")
names(apriori)<-c("a","b")