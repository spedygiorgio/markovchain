context("Checking that meanFirstPassageTime works as expected")

#load & prepare data

scr <- c("s","c","r")

Pmat <- matrix( c(6,3,1,  
                  2,3,5, 
                  4,1,5)/10,3,byrow=T)
P <- new("markovchain", states=scr, transitionMatrix=Pmat)

# Analytic solutions
P_r    <- c(s=50,c=30)/11
P_full <- matrix( c( 0,    15/4, 50/11,
                     10/3, 0,    30/11,
                     8/3,  5,    0 ), byrow=T, ncol=3)
rownames(P_full) <- scr
colnames(P_full) <- scr


Poz <- new("markovchain", states=scr, 
           transitionMatrix=matrix(c(2,1,1, 
                                     2,0,2, 
                                     1,1,2)/4, byrow=T, ncol=3)) 

Poz_full <- matrix( c( 0,  4, 10/3,
                     8/3,  0, 8/3,
                     10/3, 4, 0   ), byrow=T, ncol=3)
rownames(Poz_full) <- scr
colnames(Poz_full) <- scr

test_that("Check meanFirstPassageTime", {
  expect_equal(meanFirstPassageTime(P,"r"), P_r)
  expect_equal(meanFirstPassageTime(P),     P_full)
  expect_equal(meanFirstPassageTime(Poz),   Poz_full)
})
