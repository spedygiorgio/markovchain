#library(markovchain)

context("Bootstrap parallel")

#could be very long
seqs<-sample(x = letters[1:4],size = 200,replace = TRUE)

# mcFit<-markovchainFit(data = seqs[1:100],nboot = 50,method="bootstrap",parallel=FALSE)
mcFit<-markovchainFit(data = seqs[1:100],nboot = 100,method="bootstrap",parallel=TRUE)
# print(mcFit)

# library(rbenchmark)
# res<-benchmark(
#   markovchainFit(data = seqs[1:100],nboot = 100,method="bootstrap",parallel=TRUE),
#   markovchainFit(data = seqs[1:100],nboot = 100,method="bootstrap",parallel=FALSE)
#   , replications=10
# )
# print(res)

# library(microbenchmark)
# res<-microbenchmark(
#   markovchainFit(data = seqs[1:100],nboot = 10,method="bootstrap",parallel=FALSE),
#   markovchainFit(data = seqs[1:100],nboot = 10,method="bootstrap",parallel=TRUE)
# )
# print(res)

# data(rain)
# sequs<-rain$rain
# mcBoot<-markovchainFit(data = sequs,nboot = 10,method="bootstrap") #ok
# mcBoot2<-markovchainFit(data = sequs,nboot = 200,method="bootstrap") #ok but slower
# mcBoot<-markovchainFit(data = sequs[1:100],nboot = 100,method="bootstrap",parallel=TRUE) # ok
# mcBoot<-markovchainFit(data = sequs,nboot = 100,method="bootstrap",parallel=TRUE)
