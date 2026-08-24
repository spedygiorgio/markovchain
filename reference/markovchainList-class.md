# Non homogeneus discrete time Markov Chains class

A class to handle non homogeneous discrete Markov chains

## Slots

- `markovchains`:

  Object of class `"list"`: a list of markovchains

- `name`:

  Object of class `"character"`: optional name of the class

## Note

The class consists in a list of `markovchain` objects. It is aimed at
working with non homogeneous Markov chains.

## Objects from the Class

A `markovchainlist` is a list of `markovchain` objects. They can be used
to model non homogeneous discrete time Markov Chains, when transition
probabilities (and possible states) change by time.

## Methods

- \[\[:

  `signature(x = "markovchainList")`: extract the i-th `markovchain`

- dim:

  `signature(x = "markovchainList")`: number of `markovchain` underlying
  the matrix

- predict:

  `signature(object = "markovchainList")`: predict from a
  `markovchainList`

- print:

  `signature(x = "markovchainList")`: prints the list of markovchains

- show:

  `signature(object = "markovchainList")`: same as `print`

## References

A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall
2010

## See also

[`markovchain`](markovchain-class.md)

## Author

Giorgio Spedicato

## Examples

``` r
showClass("markovchainList")
#> Class "markovchainList" [package "markovchain"]
#> 
#> Slots:
#>                                 
#> Name:  markovchains         name
#> Class:         list    character
#define a markovchainList
statesNames=c("a","b")

mcA<-new("markovchain",name="MCA", 
         transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),
                          byrow=TRUE, nrow=2, 
                          dimnames=list(statesNames,statesNames))
        )
                                                           
mcB<-new("markovchain", states=c("a","b","c"), name="MCB",
         transitionMatrix=matrix(c(0.2,0.5,0.3,0,1,0,0.1,0.8,0.1),
         nrow=3, byrow=TRUE))
 
mcC<-new("markovchain", states=c("a","b","c","d"), name="MCC",
         transitionMatrix=matrix(c(0.25,0.75,0,0,0.4,0.6,
                                   0,0,0,0,0.1,0.9,0,0,0.7,0.3), 
                                 nrow=4, byrow=TRUE)
)
mcList<-new("markovchainList",markovchains=list(mcA, mcB, mcC), 
           name="Non - homogeneous Markov Chain")
```
