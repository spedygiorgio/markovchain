
verifyMarkovProperty<-function(mc) {
  n<-length(mc)
  u<-unique(mc)
  stateNames<-u
  nelements<-length(stateNames)
  mat<-zeros(nrow=nelements, ncol=2)
  dimnames(mat)<-list(stateNames, c("SSO", "TSO-SSO"))
  SSO<-numeric()
  for(i in 1:nelements) {
    sname<-stateNames[i]
    SSO[sname]<-0
  }
  TSO<-SSO
  out<-list()
  for(present in stateNames) {
    for(future in stateNames) {
      # print(paste0(present,'->',future))
      for(i in 1:nelements) TSO[i]<-SSO[i]<-0
      for(i in 1:(n-1))
      {
        past<-mc[i]
        if(mc[i+1] == present) {
          TSO[past] <- TSO[past] + 1
          if((i < n - 1) && (mc[i+2] == future)) {
            for(s in stateNames) {
              if(s == past) {
                SSO[s] <- SSO[s] + 1
              }
            }
          }
        }
      }
      for(i in 1:(length(SSO))) {
        mat[i,1]<-SSO[i]
        mat[i,2]<-TSO[i] - SSO[i]
      }
      # chi-squared test
      res<-chisq.test(mat)
      out[[paste0(present,future)]]<-res
    }
  }
  return(out)
}

#the out should be a function with following slots:
#stateTransitionSequenceTable (as in page 25)
#statistic: the chi - square statistic 
#p-value: the p-value of the statistic with appropriate degrees of freedom (see p 25-27 on the calculation)

#also the following should run: data(blanden); myMc<-as(blanden,"markovchain");sequenza<-rmarkovchain(n = 100,myMc)
#verifyMarkovProperty(sequenza)
#http://stats.stackexchange.com/questions/37386/check-memoryless-property-of-a-markov-chain 

assessOrder<-function(mc) {
  n<-length(mc)
  u<-unique(mc)
  states<-u
  nelements<-length(states)
  mat<-zeros(nelements)
  dimnames(mat)<-list(states, states)
  # print(mat)
  SSO<-numeric()
  TSO<-SSO
  out<-list()
  for(present in states) {
    mat<-zeros(nelements)
    for(i in 1:(n - 2)) {
      if(present == mc[i + 1]) {
        # print(paste0(mc[i],'->',mc[i+2]))
        # mat[mc[i], mc[i+2]] <- mat[mc[i], mc[i+2]] + 1
      }
      # print(paste0(present,'->',future))
    }
    #       # chi-squared test
    #       res<-chisq.test(mat)
    #       out[[paste0(present,future)]]<-res
  }
  # print(out)
  return(out)
}

assessStationarity<-function(object) {
  
}

divergenceTest<-function(object) {
  
}