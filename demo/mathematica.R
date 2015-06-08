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

#first passage time

mathematicaBMatr= matrix(c(  0, 1/2, 1/2,  1/2, 0, 1/2,  1/2, 1/2, 0),byrow=TRUE, nrow=3) ;
mathematicabMc<-as(mathematicaBMatr, "markovchain")
firstPassage(mathematicabMc, "s3",3)