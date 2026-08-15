# Easy Handling Discrete Time Markov Chains

The package contains classes and method to create and manage (plot,
print, export for example) discrete time Markov chains (DTMC). In
addition it provide functions to perform statistical (fitting and
drawing random variates) and probabilistic (analysis of DTMC
proprieties) analysis

## References

Discrete-Time Markov Models, Bremaud, Springer 1999

## See also

Useful links:

- <https://github.com/spedygiorgio/markovchain/>

- Report bugs at <https://github.com/spedygiorgio/markovchain/issues>

## Author

Giorgio Alfredo Spedicato Maintainer: Giorgio Alfredo Spedicato
\<spedicato_giorgio@yahoo.it\>

## Examples

``` r
# create some markov chains
statesNames=c("a","b")
mcA<-new("markovchain", transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),byrow=TRUE,
         nrow=2, dimnames=list(statesNames,statesNames)))
         
statesNames=c("a","b","c")
mcB<-new("markovchain", states=statesNames, transitionMatrix=
         matrix(c(0.2,0.5,0.3,0,1,0,0.1,0.8,0.1), nrow=3, 
         byrow=TRUE, dimnames=list(statesNames, statesNames)))

statesNames=c("a","b","c","d")
matrice<-matrix(c(0.25,0.75,0,0,0.4,0.6,0,0,0,0,0.1,0.9,0,0,0.7,0.3), nrow=4, byrow=TRUE)
mcC<-new("markovchain", states=statesNames, transitionMatrix=matrice)
mcD<-new("markovchain", transitionMatrix=matrix(c(0,1,0,1), nrow=2,byrow=TRUE))


#operations with S4 methods
mcA^2
#> Unnamed Markov chain^2 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  a, b 
#>  The transition matrix  (by rows)  is defined as follows: 
#>      a    b
#> a 0.52 0.48
#> b 0.16 0.84
#> 
steadyStates(mcB)
#>      a b c
#> [1,] 0 1 0
absorbingStates(mcB)
#> [1] "b"
markovchainSequence(n=20, markovchain=mcC, include=TRUE)
#>  [1] "a" "a" "b" "b" "b" "b" "b" "b" "a" "a" "a" "b" "b" "b" "b" "b" "b" "b" "b"
#> [20] "a" "b"
```
