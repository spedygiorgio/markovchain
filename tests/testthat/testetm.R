library(markovchain)
library(etm)

data(sir.cont)

# Modification for patients entering and leaving a state
# at the same date
# Change on ventilation status is considered
# to happen before end of hospital stay
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
for (i in 2:nrow(sir.cont)) {
  if (sir.cont$id[i]==sir.cont$id[i-1]) {
    if (sir.cont$time[i]==sir.cont$time[i-1]) {
      sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
    }
  }
}

### Computation of the transition probabilities
# Possible transitions.
tra <- matrix(ncol=3,nrow=3,FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE
# print(tra)

# etm
tr.prob <- etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)

tr.prob
# df<-tr.prob$trans
# df
# summary(tr.prob)

etm2mc<-as(tr.prob, "markovchain")
# etm2mc

test_that("Conversion of objects", 
          {
            expect_equal(class(etm2mc)=="markovchain",TRUE)
          })
