# Markov Chain class

The S4 class that describes `markovchain` objects.

## Slots

- `states`:

  Name of the states. Must be the same of `colnames` and `rownames` of
  the transition matrix

- `byrow`:

  TRUE or FALSE indicating whether the supplied matrix is either
  stochastic by rows or by columns

- `transitionMatrix`:

  Square transition matrix

- `name`:

  Optional character name of the Markov chain

## Note

1.  `markovchain` object are backed by S4 Classes.

2.  Validation method is used to assess whether either columns or rows
    totals to one. Rounding is used up to `.Machine$double.eps * 100`.
    If state names are not properly defined for a probability `matrix`,
    coercing to `markovchain` object leads to overriding states name
    with artificial "s1", "s2", ... sequence. In addition, operator
    overloading has been applied for \\+,\*,^,==,!=\\ operators.

## Creation of objects

Objects can be created by calls of the form
`new("markovchain", states, byrow, transitionMatrix, ...)`.

## Methods

- \*:

  `signature(e1 = "markovchain", e2 = "markovchain")`: multiply two
  `markovchain` objects

- \*:

  `signature(e1 = "markovchain", e2 = "matrix")`: markovchain by matrix
  multiplication

- \*:

  `signature(e1 = "markovchain", e2 = "numeric")`: markovchain by
  numeric vector multiplication

- \*:

  `signature(e1 = "matrix", e2 = "markovchain")`: matrix by markov chain

- \*:

  `signature(e1 = "numeric", e2 = "markovchain")`: numeric vector by
  `markovchain` multiplication

- \[:

  `signature(x = "markovchain", i = "ANY", j = "ANY", drop = "ANY")`:
  ...

- ^:

  `signature(e1 = "markovchain", e2 = "numeric")`: power of a
  `markovchain` object

- ==:

  `signature(e1 = "markovchain", e2 = "markovchain")`: equality of two
  `markovchain` object

- !=:

  `signature(e1 = "markovchain", e2 = "markovchain")`: non-equality of
  two `markovchain` object

- absorbingStates:

  `signature(object = "markovchain")`: method to get absorbing states

- canonicForm:

  `signature(object = "markovchain")`: return a `markovchain` object
  into canonic form

- coerce:

  `signature(from = "markovchain", to = "data.frame")`: coerce method
  from markovchain to `data.frame`

- conditionalDistribution:

  `signature(object = "markovchain")`: returns the conditional
  probability of subsequent states given a state

- coerce:

  `signature(from = "data.frame", to = "markovchain")`: coerce method
  from `data.frame` to `markovchain`

- coerce:

  `signature(from = "table", to = "markovchain")`: coerce method from
  `table` to `markovchain`

- coerce:

  `signature(from = "msm", to = "markovchain")`: coerce method from
  `msm` to `markovchain`

- coerce:

  `signature(from = "msm.est", to = "markovchain")`: coerce method from
  `msm.est` (but only from a Probability Matrix) to `markovchain`

- coerce:

  `signature(from = "etm", to = "markovchain")`: coerce method from
  `etm` to `markovchain`

- coerce:

  `signature(from = "sparseMatrix", to = "markovchain")`: coerce method
  from `sparseMatrix` to `markovchain`

- coerce:

  `signature(from = "markovchain", to = "igraph")`: coercing to `igraph`
  objects

- coerce:

  `signature(from = "markovchain", to = "matrix")`: coercing to `matrix`
  objects

- coerce:

  `signature(from = "markovchain", to = "sparseMatrix")`: coercing to
  `sparseMatrix` objects

- coerce:

  `signature(from = "matrix", to = "markovchain")`: coercing to
  `markovchain` objects from `matrix` one

- dim:

  `signature(x = "markovchain")`: method to get the size

- names:

  `signature(x = "markovchain")`: method to get the names of states

- names\<-:

  `signature(x = "markovchain", value = "character")`: method to set the
  names of states

- initialize:

  `signature(.Object = "markovchain")`: initialize method

- plot:

  `signature(x = "markovchain", y = "missing")`: plot method for
  `markovchain` objects

- predict:

  `signature(object = "markovchain")`: predict method

- print:

  `signature(x = "markovchain")`: print method.

- show:

  `signature(object = "markovchain")`: show method.

- sort:

  `signature(x = "markovchain", decreasing=FALSE)`: sorting the
  transition matrix.

- states:

  `signature(object = "markovchain")`: returns the names of states (as
  `names`.

- steadyStates:

  `signature(object = "markovchain")`: method to get the steady vector.

- summary:

  `signature(object = "markovchain")`: method to summarize structure of
  the markov chain

- transientStates:

  `signature(object = "markovchain")`: method to get the transient
  states.

- t:

  `signature(x = "markovchain")`: transpose matrix

- transitionProbability:

  `signature(object = "markovchain")`: transition probability

## References

A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall
2010

## See also

[`markovchainSequence`](markovchainSequence.md),[`markovchainFit`](markovchainFit.md)

## Author

Giorgio Spedicato

## Examples

``` r
#show markovchain definition
showClass("markovchain")
#> Class "markovchain" [package "markovchain"]
#> 
#> Slots:
#>                                                                           
#> Name:            states            byrow transitionMatrix             name
#> Class:        character          logical           matrix        character
#create a simple Markov chain
transMatr<-matrix(c(0.4,0.6,.3,.7),nrow=2,byrow=TRUE)
simpleMc<-new("markovchain", states=c("a","b"),
              transitionMatrix=transMatr, 
              name="simpleMc")
#power
simpleMc^4
#> simpleMc^4 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  a, b 
#>  The transition matrix  (by rows)  is defined as follows: 
#>        a      b
#> a 0.3334 0.6666
#> b 0.3333 0.6667
#> 
#some methods
steadyStates(simpleMc)
#>              a         b
#> [1,] 0.3333333 0.6666667
absorbingStates(simpleMc)
#> character(0)
simpleMc[2,1]
#> [1] 0.3
t(simpleMc)
#> Unnamed Markov chain 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  a, b 
#>  The transition matrix  (by cols)  is defined as follows: 
#>     a   b
#> a 0.4 0.3
#> b 0.6 0.7
#> 
is.irreducible(simpleMc)
#> [1] TRUE
#conditional distributions
conditionalDistribution(simpleMc, "b")
#>   a   b 
#> 0.3 0.7 
#example for predict method
sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit<-markovchainFit(data=sequence)
predict(mcFit$estimate, newdata="b",n.ahead=3)
#> [1] "a" "b" "a"
#direct conversion
myMc<-as(transMatr, "markovchain")

#example of summary
summary(simpleMc)
#> simpleMc  Markov chain that is composed by: 
#> Closed classes: 
#> a b 
#> Recurrent classes: 
#> {a,b}
#> Transient classes: 
#> NONE 
#> The Markov chain is irreducible 
#> The absorbing states are: NONE
if (FALSE) plot(simpleMc) # \dontrun{}
```
