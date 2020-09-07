#library(markovchain)
context("Checking conversion of objects to etm")

check_etm_availability <- function(){
  is_available <-  require("etm")
  if (!is_available) {
    skip("etm package unavailable")
  }
}

get_etm_transition_matrix <- function() {
  require(etm)
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
  tr.prob <- etm::etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)
  return(tr.prob)
}

# tr.prob
# df<-tr.prob$trans
# df
# summary(tr.prob)
# etm2mc

test_that("Conversion of objects", 
          {
            check_etm_availability() #check package availablity
            obj_to_test <- get_etm_transition_matrix() #obtain the etm obj
            etm2mc<-as(obj_to_test, "markovchain") #try to convert
            expect_equal(class(etm2mc)=="markovchain",TRUE)
          })
