
verifyMarkovProperty<-function(mc) {
  n<-length(mc)
  u<-unique(mc)
  stateNames<-u
  nelements<-length(stateNames)
  mat<-zeros(nelements)
  dimnames(mat)<-list(stateNames, c("SSO", "TSO-SSO"))
  SSO<-numeric()
  for(i in 1:nelements) {
    sname<-stateNames[i]
    SSO[sname]<-0
  }
  TSO<-SSO
  prob<-SSO
  out<-list()
  for(present in stateNames) {
    for(future in stateNames) {
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
        prob[i]<-SSO[i]/TSO[i]
        mat[i,0]<-SSO[i]
        mat[i,1]<-TSO[i] - SSO[i]
      }
      
      # chi-squared test
      res<-chisq.test(mat)
      out[[paste0(present,future)]]<-res
      # out[[length(out) + 1]]<-res
    }
  }
  
  return(out)
}

assessOrder<-function(object) {
  return(1)
}

assessStationarity<-function(object) {
  
}

divergenceTest<-function(object) {
  
}