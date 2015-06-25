setClass("ctmc",
         representation(states = "character",
                        byrow = "logical",
                        generator = "matrix",
                        name = "character"),
         prototype(states = c("a", "b"), byrow = TRUE,
                   generator = matrix(data = c(-1, 1, 1, -1), byrow = TRUE, nrow = 2,
                                      dimnames = list(c("a", "b"), c("a", "b"))),
                   name = "Unnamed CTMC")
)

setMethod("initialize",
          signature(.Object = "ctmc"),
          function (.Object, states, byrow, generator,name,...) 
          {
            # put the standard markovchain 
            if(missing(generator)) generator=matrix(data=c(-1,1,1,-1), #create a naive matrix
                                                    nrow=2,
                                                    byrow=TRUE, 
                                                    dimnames=list(c("a","b"), c("a","b"))
            )
            
            # check names of transition matrix
            if(all(is.null(rownames(generator)), is.null(colnames(generator)))==TRUE) { #if all names are missing it initializes them to "1", "2",...
              if(missing(states)) {
                nr=nrow(generator)
                stateNames<-as.character(seq(1:nr))
              } else {stateNames=states}
              
              rownames(generator)=stateNames
              colnames(generator)=stateNames
            } else if(is.null(rownames(generator))) { #fix when rownames null
              rownames(generator)=colnames(generator)
            } else if(is.null(colnames(generator))) { #fix when colnames null
              colnames(generator)=rownames(generator)
            } else if(!setequal(rownames(generator),colnames(generator)))  colnames(generator)=rownames(generator) #fix when different
            if(missing(states)) states=rownames(generator) #assign
            if(missing(byrow)) byrow=TRUE #set byrow as true by default
            if(missing(name)) name="Unnamed Markov chain"
            callNextMethod(.Object, states = states, byrow = byrow, generator=generator,name=name,...)
          }
)

setMethod("states","ctmc", 
          function(object) {
            out <- object@states
            return(out)
          }
)

setMethod("dim","ctmc", 
          function(x) {
            out <- nrow(x@generator)
            return(out)
          }
)

setValidity("ctmc",
            function(object) {
              check<-NULL
              # performs a set of check whose results are saved in check
              if (.isGenRcpp(object@generator)==FALSE) check <- "Error! Not a generator matrix" 
              if (object@byrow==TRUE) {
                if(any(round(rowSums(object@generator),5)!=0)) check <- "Error! Row sums not equal to zero"
              } else {
                if(any(round(colSums(object@generator),5)!=0)) check <- "Error! Col sums not equal to zero"
              } 
              if (nrow(object@generator)!=ncol(object@generator)) check <- "Error! Not squared matrix" #check if squalre matrix
              if (!setequal(colnames(object@generator),object@states)) check <- "Error! Colnames <> states" #checks if 
              if (!setequal(rownames(object@generator),object@states)) check <- "Error! Rownames <> states"
              if ( is.null(check) ) return(TRUE) else return(check)
            }
)
