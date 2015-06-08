#create the first Mathematica 9 example matrix: STRUCTURAL PROPRIETIES
mathematicaAMatr <- matrix(c(0, 1/3, 0, 2/3, 0, 
                            1/2, 0, 0, 0, 1/2, 
                            0, 0, 1/2, 1/2, 0, 
                            0, 0, 1/2, 1/2, 0, 
                            0, 0, 0, 0, 1), byrow = TRUE,
                          nrow=5)
mathematicaAMc<-as(mathematicaAMatr, "markovchain")
summary(mathematicaAMc)
canonicForm(mathematicaAMc)
is.irreducible(mathematicaAMc)
transientStates(mathematicaAMc)
absorbingStates(mathematicaAMc)

#functions above should behave like 
#http://www.wolfram.com/mathematica/new-in-9/markov-chains-and-queues/structural-properties-of-finite-markov-processes.html

#first passage time
#@TAE: check this: 
mathematicaBMatr= matrix(c(  0, 1/2, 1/2,  1/2, 0, 1/2,  1/2, 1/2, 0),byrow=TRUE, nrow=3) ;
mathematicabMc<-as(mathematicaBMatr, "markovchain")
firstPassage(mathematicabMc, "s3",3) #if you repeat this more thime results change. No good
#should behave like
#http://www.wolfram.com/mathematica/new-in-9/markov-chains-and-queues/distribution-of-times-to-reach-a-target-state.html