#Bard PPT examples
#Page 8
#require(matlab)
statesNames=as.character(0:3)
pg8<-markovchain::zeros(4)
pg8[1,c(2,3)]<-0.5
pg8[2,1]<-1
pg8[3,4]<-1
pg8[4,c(1,4)]<-0.5
pg8Mc<-new("markovchain",transitionMatrix=pg8,states=statesNames, name="Page 8")
summary(pg8Mc)
#Page 9
statesNames=as.character(0:4)
pg9<-markovchain::zeros(5)
pg9[c(1,2),c(1,2)]<-0.5
pg9[3,3]<-1
pg9[4,c(3,4)]<-1/2
pg9[5,1]<-1
pg9Mc<-new("markovchain",
           transitionMatrix=pg9,
           states=statesNames, 
           name="Page 9")
summary(pg9Mc)
#Page 10
statesNames=as.character(0:3)
pg10<-markovchain::zeros(4)
pg10[1,c(2,3)]<-1/2
pg10[c(2,3),4]<-1
pg10[4,1]<-1
pg10Mc<-new("markovchain",
           transitionMatrix=pg10,
           states=statesNames, 
           name="Page 10")
summary(pg10Mc)
#Page 11
statesNames=as.character(1:5)
pg11<-markovchain::zeros(5)
pg11[1,c(1,2)]<-c(0.4,0.6)
pg11[2,c(1,2)]<-c(0.5,0.5)
pg11[3,c(3,4)]<-c(0.3,0.7)
pg11[4,c(3,4,5)]<-c(5,4,1)/10
pg11[5,c(4,5)]<-c(0.8,0.2)

pg11Mc<-new("markovchain",
            transitionMatrix=pg11,
            states=statesNames, 
            name="Page 10")
summary(pg11Mc)
