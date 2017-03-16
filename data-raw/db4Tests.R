#esempi di check fitting markov chains x tests
ciao<-c("a","a","b","b","a",NA,"b","a","b","a","a")
simpleMcCiaoFit<-markovchainFit(ciao)
data(rain)
checksAlofiRawTransitions<-createSequenceMatrix(rain$rain)

#check by matrix

data(holson)
myHolson<-as.matrix(holson[,-1]); rownames(myHolson)<-holson$id
checkmarkovchainFitList<-markovchainListFit(data=myHolson)


devtools::use_data(simpleMcCiaoFit,checksAlofiRawTransitions,checkmarkovchainFitList,internal = TRUE,overwrite = TRUE)

