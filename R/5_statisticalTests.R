
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
  states<-unique(mc)
  nelements<-length(states)
  # mat<-zeros(nelements)
  # dimnames(mat)<-list(states, states)
  out<-list()
  for(present in states) {
    mat<-zeros(nelements)
    dimnames(mat)<-list(states, states)
    for(i in 1:(n - 2)) {
      if(present == mc[i + 1]) {
        past<-mc[i]
        future<-mc[i+2]
        # print(paste0(past,'->',future))
        mat[past, future] <- mat[past, future] + 1 
      }
    }
    # chi-squared test
    res<-chisq.test(mat)
    out[[present]]<-res
  }
  return(out)
}

assessStationarity<-function(mc) {
  n<-length(mc)
  states<-unique(mc)
  nstates<-length(states)
  out<-list()
  times<-numeric(nstates)
  names(times)<-states
  matlist<-list()
  for(state in states) {
    mat<-zeros(n-1,nstates)
    dimnames(mat)<-list(1:(n-1), states)
    matlist[[state]]<-mat
  }
  for(t in 1:(n - 1)) {
    past<-mc[t]
    present<-mc[t + 1]
    for(state in states) {
      mat<-matlist[[state]]
      for(s in states) {
        if(t > 1) 
          mat[t, s] <- mat[t-1, s]
        if((state == past) && (s == present)) {
          if(t == 1) mat[t, s] <- 1
          else mat[t, s] <- mat[t-1, s] + 1
        }
      }
      matlist[[state]]<-mat
    }
  }
  for(s in states) {
    mat<-matlist[[s]]
    rowsums<-rowSums(mat)
    indices<-which(rowsums == 0)
    mat<-mat/rowsums
    for(i in indices) mat[i,]<-1/nstates
    # chi-squared test
    res<-chisq.test(mat)
    out[[s]]<-res
  }
  return(out)
}

divergenceTest<-function(m1, m2, mc) {
  n<-length(mc)
  M<-nrow(m1)
  v<-numeric()
  out<-2*n/.phi2(1)
  sum<-0
  for(i in 1:M) {
    sum2<-0
    sum3<-0
    for(j in 1:M) {
      sum2<-sum2+m2[i,j]*.phi(m1[i,j]/m2[i,j])
      if(j > 1 && mc[j-1] == i)
        sum3<-sum3 + 1
    }
    v[i]<-sum3
    sum<-v[i]/n*sum2
  }
  out<-out*sum
  return (out)
}

.phi<-function(x) {
  out<-x*log(x)-x+1
  return(out)
}

.phi2<-function(x) {
  out<-1/x
  return(out)
}
