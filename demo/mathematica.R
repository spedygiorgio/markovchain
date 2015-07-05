#STRUCTURAL PROPRIETIES
#see mathematica 9
mathematicaAMatr <- matrix(c(0, 1/3, 0, 2/3, 0, 
                            1/2, 0, 0, 0, 1/2, 
                            0, 0, 1/2, 1/2, 0, 
                            0, 0, 1/2, 1/2, 0, 
                            0, 0, 0, 0, 1), byrow = TRUE,
                          nrow=5)
#should behave like
#http://www.wolfram.com/mathematica/new-in-9/markov-chains-and-queues/structural-properties-of-finite-markov-processes.html
mathematicaAMc<-as(mathematicaAMatr, "markovchain")
summary(mathematicaAMc)
canonicForm(mathematicaAMc)
is.irreducible(mathematicaAMc)
transientStates(mathematicaAMc)
absorbingStates(mathematicaAMc)

communicatingClasses(mathematicaAMc)
recurrentClasses(mathematicaAMc)

#first passage time
#should behave like
#http://www.wolfram.com/mathematica/new-in-9/markov-chains-and-queues/distribution-of-times-to-reach-a-target-state.html
mathematicaBMatr= matrix(c(  0, 1/2, 1/2,  1/2, 0, 1/2,  1/2, 1/2, 0),byrow=TRUE, nrow=3) ;
mathematicabMc<-as(mathematicaBMatr, "markovchain")
firstPassage(mathematicabMc, "s3",9) #if you repeat this more thime results change. No good

communicatingClasses(mathematicabMc)
recurrentClasses(mathematicabMc)