#PSEUDO BAYESIAN ESTIMATOR

getAlphaVector <- function(raw, apriori){
  v_i <- rowSums(raw) 
  K_i <- numeric(nrow(raw))
  sumSquaredY <- rowSums(raw^2)
  #get numerator
  K_i_num <- v_i^2-sumSquaredY
  #get denominator
  VQ <- matrix(0,nrow= nrow(apriori),ncol=ncol(apriori))
  for (i in 1:nrow(VQ)) {
    VQ[i,]<-v_i[i]*apriori[i,]
  }
  
  K_i_den<-rowSums((raw - VQ)^2)
  
  K_i <- K_i_num/K_i_den
  
  #get the alpha vector
  alpha <- K_i / (v_i+K_i)
  
  #empirical transition matrix
  Emp<-raw/rowSums(raw)
  
  #get the estimate
  out<-matrix(0, nrow= nrow(raw),ncol=ncol(raw))
  for (i in 1:nrow(out)) {
    out[i,]<-alpha[i]*apriori[i,]+(1-alpha[i])*Emp[i,]
  }
  return(out)
}