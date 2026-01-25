#up and down markov chains

mcUpDown<-new("markovchain", states=c("up","down"), 
transitionMatrix=matrix(c(0,1,1,0),nrow=2, byrow=TRUE),name="UpDown")

#gamblers ruin

gamblerRuinMarkovChain<-function(moneyMax, prob=0.5) {
  matr<-markovchain::zeros(moneyMax+1)
  states<-as.character(seq(from=0, to=moneyMax, by=1))
  rownames(matr)=states; colnames(matr)=states
  matr[1,1]=1;matr[moneyMax+1,moneyMax+1]=1
  for(i in 2:moneyMax)
  {
    matr[i,i-1]=1-prob;matr[i,i+1]=prob
  }
  out<-new("markovchain",  
           transitionMatrix=matr, 
           name=paste("Gambler ruin",moneyMax,"dim",sep=" ")
           )
  return(out)
}