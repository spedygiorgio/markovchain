context("Checking that steadyStates works as expected")

#load & prepare data
data(rain)

mc1<-as(matrix(c(0.5 , 0.5 , 0 , 0,0.3 , 0.3 , 0.4 ,0,0 , 0.5 , 0.5 , 0,0 , 0 , 0.5 , 0.5),byrow=TRUE,nrow=4),"markovchain")
mc2<-as(matrix(c(0,0.5,0.50,0.0,0,0.00,0,0.0,0.99,0.0,0,0.01,0,0.8,0.00,0.2,0,0.00,0,0.0,0.00,
                 0.0,1,0.00,0,0.0,0.00,1.0,0,0.00,0,0.0,0.00,0.0,0,1.00), nrow=6, byrow=TRUE),"markovchain")
mc3<-matrix(0, nrow=5, ncol=5)
mc3[1:2,1:2]<-matrix(c(0.4,0.6,.5,.5),nrow=2, byrow=TRUE)
mc3[3:5,3:5]<-matrix(c(0.3,0.7,0,.5,.4,.1,0,.8,.2),nrow=3, byrow=TRUE)
mc3<-as(mc3,"markovchain")
mcRain<-steadyStates(markovchainFit(data=rain$rain)$estimate)


mc5<-matrix(0, nrow=4,ncol=4)
mc5[1,] <- c(0.5,0.5,0,0)
mc5[2,] <- c(0.5,0.5,0,0)
mc5[3,] <- c(1/3,1/6,1/6,1/3)
mc5[4,4] <- 1
mc5<-as(mc5,"markovchain")


mathematicaMatr <- matrix(0, nrow=5,ncol=5)
mathematicaMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
                     name = "Mathematica MC", states = statesNames)



test_that("Check steadyStates", {
  expect_equal(steadyStates(mc1), steadyStates1)
  expect_equal(steadyStates(mc2), steadyStates2)
  expect_equal(steadyStates(mc3), steadyStates3)
  expect_equal(mcRain, steadyStates4)
  expect_equal(steadyStates(mc5), steadyStates5)
  expect_equal(steadyStates(mathematicaMc), steadyStates6)
})
