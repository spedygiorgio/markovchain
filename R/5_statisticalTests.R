
verifyMarkovProperty<-function(object) {
  n<-length(object)
  u<-unique(object)
  stateNames<-u
  nelements<-length(stateNames)
  mat<-zeros(nelements)
  dimnames(mat)<-list(stateNames, c("SSO", "TSO-SSO"))
  SSO<-numeric()
  for(i in 1:nelements) {
    sname<-stateNames[i]
    SSO[sname]<-0
  }
  # print(SSO)
  TSO<-SSO
  # print(TSO)
  prob<-SSO
  present<-"a"
  future<-"b"
  for(i in 1:(n-2))
  {
    # print(object[i])
    past<-object[i]
    if(object[i+1] == present) {
      TSO[past] <- TSO[past] + 1
      if(object[i+2] == future) {
        for(s in stateNames) {
          # print(s)
          if(s == past) {
            # print(paste0(s,"->",present,"->",future))
            SSO[s] <- SSO[s] + 1
          }
        }
      }
    }
  }
  for(i in 1:(length(SSO))) {
    prob[i]<-SSO[i]/TSO[i]
    mat[i,0]<-SSO[i]
    mat[i,1]<-TSO[i] - SSO[i]
  }
  # print(mat)
  
  # chi-squared test
  res<-chisq.test(mat)
  # print(res)
  
  return(res)
}

assessOrder<-function(object) {
  return(1)
}

assessStationarity<-function(object) {
  
}

divergenceTest<-function(object) {
  
}