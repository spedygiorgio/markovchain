
verifyMarkovProperty<-function(object) {
  significanceLevel<-0.05
  # print(object)
  n<-length(object)
  u<-unique(object)
  # print(u)
  stateNames<-u
  nelements<-length(stateNames)
  mat<-zeros(nelements)
  dimnames(mat)<-list(stateNames, stateNames)
  # print(mat)
  SSO<-c("a"=0,"b"=0)
  TSO<-c("a"=0,"b"=0)
  prob<-c("a"=0,"b"=0)
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
#   print(SSO)
#   print(TSO)
  for(i in 1:(length(SSO))) {
    prob[i]<-SSO[i]/TSO[i]
  }
  # print(prob)
  
  # chi-square test
  Q<-0
  for(i in 1:2) {
    for(j in 1:2) {
      # Q<-Q+(N[i,j] - n*(rowsum[i]/n)*(colsum[j]))^2 / (n*(rowsum[i]/n)*(colsum[j]/n))
    }
  }
  df<-length(stateNames)
  # print(df)
  # res<-chisq.test(object)
  # print(res)
  # chisquare<-P_X(Q)
  # if(chisquare > significanceLevel) return(TRUE)
  
  # P<-object$estimate@transitionMatrix
  # P<-object@transitionMatrix
  # stateNames<-states(object)
  #   print(P)
  #   n<-nrow(P)
  #   j<-1
  #   m<-2
  #   for (i in 1:n) {
  #     if(P[i,j] > 0) {
  #       if(P[j,m] > 0) {
  #         r<-P[i,j] * P[j,m] / P[i,j]
  #         print(r)
  #       }
  #     }
  #   }
  #   colnames(outMatr)<-stateNames
  #   rownames(outMatr)<-1:n
  return(FALSE)
}

assessOrder<-function(object) {
  return(1)
}

assessStationarity<-function(object) {
  
}

divergenceTest<-function(object) {
  
}