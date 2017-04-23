#FITTING MARKOV CHAINS

#esempi di check fitting markov chains x tests
ciao<-c("a","a","b","b","a",NA,"b","a","b","a","a")
simpleMcCiaoFit<-markovchainFit(ciao)
data(rain)
checksAlofiRawTransitions<-createSequenceMatrix(rain$rain)

#check by matrix

data(holson)
myHolson<-as.matrix(holson[,-1]); rownames(myHolson)<-holson$id
checkmarkovchainFitList<-markovchainListFit(data=myHolson)
#devtools::use_data(simpleMcCiaoFit,checksAlofiRawTransitions,checkmarkovchainFitList,internal = TRUE,overwrite = TRUE)

#STEADY STATE ANALYSIS


mc1<-as(matrix(c(0.5 , 0.5 , 0 , 0,0.3 , 0.3 , 0.4 ,0,0 , 0.5 , 0.5 , 0,0 , 0 , 0.5 , 0.5),byrow=TRUE,nrow=4),"markovchain")
steadyStates1<-steadyStates(mc1)

mc2<-as(matrix(c(0,0.5,0.50,0.0,0,0.00,0,0.0,0.99,0.0,0,0.01,0,0.8,0.00,0.2,0,0.00,0,0.0,0.00,
                 0.0,1,0.00,0,0.0,0.00,1.0,0,0.00,0,0.0,0.00,0.0,0,1.00), nrow=6, byrow=TRUE),"markovchain")
steadyStates2<-steadyStates(mc2)

mc3<-matrix(0, nrow=5, ncol=5)
mc3[1:2,1:2]<-matrix(c(0.4,0.6,.5,.5),nrow=2, byrow=TRUE)
mc3[3:5,3:5]<-matrix(c(0.3,0.7,0,.5,.4,.1,0,.8,.2),nrow=3, byrow=TRUE)
mc3<-as(mc3, "markovchain")
steadyStates3<-steadyStates(mc3)

data(rain)
temp<-markovchainFit(data=rain$rain)
steadyStates4<-steadyStates(temp$estimate)

mc5<-matrix(0, nrow=4,ncol=4)
mc5[1,] <- c(0.5,0.5,0,0)
mc5[2,] <- c(0.5,0.5,0,0)
mc5[3,] <- c(1/3,1/6,1/6,1/3)
mc5[4,4] <- 1
mc5<-as(mc5,"markovchain")
steadyStates5<-steadyStates(mc5)


mathematicaMatr <- matrix(0, nrow=5,ncol=5)
mathematicaMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
                     name = "Mathematica MC", states = statesNames)


steadyStates6<-steadyStates(mathematicaMc)

#SAVING

devtools::use_data(simpleMcCiaoFit,checksAlofiRawTransitions,checkmarkovchainFitList,
                   steadyStates1,steadyStates2,steadyStates3,steadyStates4,steadyStates5,steadyStates6,
                   internal = TRUE,overwrite = TRUE)