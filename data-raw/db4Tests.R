#esempi di check fitting markov chains x tests
ciao<-c("a","a","b","b","a",NA,"b","a","b","a","a")
simpleMcCiaoFit<-markovchainFit(ciao)
data(rain)
checksAlofiRawTransitions<-createSequenceMatrix(rain$rain)
devtools::use_data(simpleMcCiaoFit,checksAlofiRawTransitions,internal = TRUE)