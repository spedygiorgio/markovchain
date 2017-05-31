library(microbenchmark)
library(markovchain)

#using the rain data sequence
data(rain)
rainSequence<-rain$rain
#choosing different sample size
sizes<-c(10,50,100,250,500,1096)
#estimating microseconds
microseconds<-numeric(length(sizes))
for(i in 1:length(sizes)) {
  mydim<-sizes[i]
  mysequence<-rainSequence[1:mydim]
  
  out<-microbenchmark(
    myFit<-markovchainFit(data=mysequence)
  )
  microseconds[i]<-mean(out$time)
}
#plot(sizes, microseconds,type="o",col="steelblue",xlab="character sequence size",ylab="microseconds",main="Computational time vs size")
