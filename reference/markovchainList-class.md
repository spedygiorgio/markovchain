# Non-homogeneous discrete time Markov chain class

The S4 class of non-homogeneous discrete time Markov chains: a sequence
of [`markovchain`](markovchain-class.md) objects, the \\k\\-th giving
the transition probabilities from time \\k - 1\\ to time \\k\\.

## Slots

- `markovchains`:

  A list of [`markovchain`](markovchain-class.md) objects.

- `name`:

  Optional name of the object.

## Note

The class consists in a list of `markovchain` objects. It is aimed at
working with non homogeneous Markov chains.

## Objects from the Class

Objects are created with
`new("markovchainList", markovchains = list(...))` or estimated by
[`markovchainListFit()`](markovchainListFit.md). The transition
probabilities, and the states, may change from one time to the next.

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
