#moved in Rcpp

# .commclassesKernel <- function(P){
# 	m <- ncol(P)
# 	stateNames <- rownames(P)
# 	T <- zeros(m) 
# 	i <- 1
# 	while (i<=m) { 
# 		a=i 
# 		b<-zeros(1,m)
# 		b[1,i]<-1
# 		old<-1
# 		new<-0
# 		while (old != new) {
# 			old <- sum(find(b>0))
# 			n <- size(a)[2]
# 			matr <- matrix(as.numeric(P[a,]),ncol=m,nrow=n) #fix
# 			c <- colSums(matr)
# 			d <- find(c)
# 			n <- size(d)[2]
# 			b[1,d] <- ones(1,n)
# 			new <- sum(find(b>0))
# 			a<-d
# 		}
# 		T[i,] <- b
# 		i <- i+1 
# 	}
# 	F <- t(T)  
# 	C <- (T>0)&(F>0)
# 	v <- (apply(t(C)==t(T),2,sum)==m)
# 	colnames(C) <- stateNames
# 	rownames(C) <- stateNames
# 	names(v) <- stateNames
# 	out <- list(C=C,v=v)
# 	return(out)
# }

#@ Tae: to be fully moved in Rcpp
#returns the underlying communicating classes
# .communicatingClasses<-function(adjMatr)
# {
#   len <- dim(adjMatr)[1]
#   classesList <- list()
#   for(i in 1:len)
#   {
#     row2Check <- adjMatr[i,]
#     proposedCommClass <- names(which(row2Check==TRUE))
#     if(i>1) 
#     {
#       for(j in 1:(length(classesList)))
#       {
#         check <- FALSE
#         check <- setequal(classesList[[j]],proposedCommClass)
#         if(check==TRUE) {proposedCommClass<-NULL; break}
#       }
#     }
#     if(!is.null(proposedCommClass)) classesList[[length(classesList)+1]]<-proposedCommClass
#   }
#   return(classesList)
# }

#@ Tae: to be fully moved in Rcpp
#communicating states
# .commStatesFinder<-function(matr)
# {
#   #Reachability matrix
#   dimMatr<-dim(matr)[1]
#   Z<-sign(matr)
#   temp<-(eye(dimMatr)+Z)%^%(dimMatr-1)
#   R<-sign(temp)
#   return(R)
# }

#moved in Rcpp

is.accessible<-function(object, from, to)
{
  out<-FALSE
  statesNames<-states(object)
  fromPos<-which(statesNames==from)
  toPos<-which(statesNames==to)
#   R<-.commStatesFinder(object@transitionMatrix)
  R<-.commStatesFinderRcpp(object@transitionMatrix)
  if(R[fromPos,toPos]==TRUE) out<-TRUE
  return(out)
}

#moved in Rcpp

#a markov chain is irreducible if is composed by only one communicating class
is.irreducible<-function(object)
{
  out<-FALSE
#   tocheck<-.communicatingClasses(.commclassesKernel(object@transitionMatrix)$C)
  tocheck<-.communicatingClassesRcpp(object)
  if(length(tocheck)==1) out<-TRUE
  return(out)
}

#moved in Rcpp

# .summaryKernel<-function(object)
# {
#   matr<-object@transitionMatrix
#   temp<-.commclassesKernelRcpp(matr)
#   communicatingClassList<-.communicatingClassesRcpp(temp$C)
#   transientStates<-names(which(temp$v==FALSE))
#   closedClasses<-list()
#   transientClasses<-list()
#   for(i in 1:length(communicatingClassList))
#   {
#     class2Test<-communicatingClassList[[i]]
#     if(length(intersect(class2Test,transientStates))>0) transientClasses[[length(transientClasses)+1]]<-class2Test else closedClasses[[length(closedClasses)+1]]<-class2Test
#   }
#   summaryMc<-list(closedClasses=closedClasses, 
#                   transientClasses=transientClasses)
#   return(summaryMc)
# }

#moved in Rcpp
#here the kernel function to compute the first passage
# .firstpassageKernel<-function(P,i,n){
#   G<-P
#   H <- matrix(NA, ncol=dim(P)[2], nrow=n) #here Thoralf suggestion
#   H[1,]<-P[i,] #initializing the first row
#   E<-1-diag(size(P)[2])
#   for (m in 2:n) {
#     G<-P%*%(G*E)
#     #H<-rbind(H,G[i,]) #removed thanks to Thoralf 
#     H[m,] <- G[i,] #here Thoralf suggestion
#   }
#   return(H)
# }

firstPassage<-function(object,state,n)
{
  P<-object@transitionMatrix
  stateNames<-states(object)
  i<-which(stateNames==state)
#   outMatr<-.firstpassageKernel(P=P,i=i,n=n)
  outMatr<-.firstpassageKernelRcpp(P=P,i=i,n=n)
  colnames(outMatr)<-stateNames
  rownames(outMatr)<-1:n
  return(outMatr)
}
#periodicity


# greatest common denominator: to be moved in Rcpp
# .gcd = function(f,s) {
#   
#   f <- abs(f)
#   s <- abs(s)
#   
# 	n <- min(f,s)
# 	N <- max(f,s)
#   
# 	if (n==0) {
# 		g=N
# 	}
# 	else {
# 		u=1
# 		while (u!=0) {
# 			u=N%%n
# 			if (u==0) {
# 				g=n
# 			}
# 			N=n
# 			n=u
# 		}
# 	}
# 	return(g)
# }

#moved in Rcpp

#function to  get the period of a DTMC
# period<-function(object) {
# 	check<-is.irreducible(object)
# 	if(check==FALSE){
# 		warning("The matrix is not irreducible")
# 		return(0)
# 	} else {
# 	P<-object@transitionMatrix
# 	n=size(P,2)
# 	v=zeros(1,n)
# 	v[1,1]=1
# 	w=numeric()
# 	d=0
# 	T=c(1)
# 	m=size(T,2)
# 	while (m>0 & d!=1) {
# 		i <- T[1]
# 		T <- T[-1]
# 		w <- c(w,i)
# 		j <- 1
# 		while (j<=n) {
# 			if (P[i,j]>0) {
# 				r=c(w,T)
# 				k=sum(r==j)
# 				if (k>0) {
# 					b=v[1,i]+1-v[1,j]
# 					d=.gcdRcpp(d,b)
# 				}
# 				else {
# 					T=c(T,j)
# 					v[1,j]=v[1,i]+1
# 				}
# 			}
# 			j=j+1
# 		}
# 		m=size(T,2)
# 	}
# 	v=v%%d
# 	return(d)
# 	}
# }

communicatingClasses<-function(object) {
  out<-.communicatingClassesRcpp(object)
  return(out)
}

recurrentClasses<-function(object) {
  out<-.recurrentClassesRcpp(object)
  return(out)
}

#' @title Check if a DTMC is regular
#' 
#' @description Function to check wether a DTCM is regular
#' 
#' @details A regular Markov chain has $A^n$ strictly positive for some n. 
#' So we check: if there is only one eigenvector; if the steadystate vector is striclty positive.
#' 
#' @param object a markovchain object
#' 
#' @return A boolean value
#' 
#' @examples 
#' P=matrix(c(0.5,.25,.25,.5,0,.5,.25,.25,.5),nrow = 3)
#' colnames(P)<-rownames(P)<-c("R","N","S")
#' ciao<-as(P,"markovchain")
#' is.regular(ciao)
#' 
#' @seealso \code{\link{is.irreducible}}


is.regular<-function(object) {
  eigenValues<-steadyStates(object = object)
  minDim<-min(dim(eigenValues))
  out <- minDim==1 & all(eigenValues>0)
  return(out)
}