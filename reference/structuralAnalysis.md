# Various function to perform structural analysis of DTMC

These functions return absorbing and transient states of the
`markovchain` objects.

## Usage

``` r
period(object)

communicatingClasses(object)

recurrentClasses(object)

transientClasses(object)

transientStates(object)

recurrentStates(object)

absorbingStates(object)

canonicForm(object)
```

## Arguments

- object:

  A `markovchain` object.

## Value

- `period`:

  returns a integer number corresponding to the periodicity of the
  Markov chain (if it is irreducible)

- `absorbingStates`:

  returns a character vector with the names of the absorbing states in
  the Markov chain

- `communicatingClasses`:

  returns a list in which each slot contains the names of the states
  that are in that communicating class

- `recurrentClasses`:

  analogously to `communicatingClasses`, but with recurrent classes

- `transientClasses`:

  analogously to `communicatingClasses`, but with transient classes

- `transientStates`:

  returns a character vector with all the transient states for the
  Markov chain

- `recurrentStates`:

  returns a character vector with all the recurrent states for the
  Markov chain

- `canonicForm`:

  returns the Markov chain reordered by a permutation of states so that
  we have blocks submatrices for each of the recurrent classes and a
  collection of rows in the end for the transient states

## References

Feres, Matlab listing for markov chain.

## See also

[`markovchain`](markovchain-class.md)

## Author

Giorgio Alfredo Spedicato, Ignacio Cordón

## Examples

``` r
statesNames <- c("a", "b", "c")
mc <- new("markovchain", states = statesNames, transitionMatrix =
          matrix(c(0.2, 0.5, 0.3,
                   0,   1,   0,
                   0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
                 dimnames = list(statesNames, statesNames))
         )

communicatingClasses(mc)
#> [[1]]
#> [1] "a" "c"
#> 
#> [[2]]
#> [1] "b"
#> 
recurrentClasses(mc)
#> [[1]]
#> [1] "b"
#> 
recurrentClasses(mc)
#> [[1]]
#> [1] "b"
#> 
absorbingStates(mc)
#> [1] "b"
transientStates(mc)
#> [1] "a" "c"
recurrentStates(mc)
#> [1] "b"
canonicForm(mc)
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  b, a, c 
#>  The transition matrix  (by rows)  is defined as follows: 
#>     b   a   c
#> b 1.0 0.0 0.0
#> a 0.5 0.2 0.3
#> c 0.8 0.1 0.1
#> 

# periodicity analysis
A <- matrix(c(0, 1, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0, 1, 0), 
            nrow = 4, ncol = 4, byrow = TRUE)
mcA <- new("markovchain", states = c("a", "b", "c", "d"), 
          transitionMatrix = A,
          name = "A")

is.irreducible(mcA) #true
#> [1] TRUE
period(mcA) #2
#> [1] 2

# periodicity analysis
B <- matrix(c(0, 0, 1/2, 1/4, 1/4, 0, 0,
                   0, 0, 1/3, 0, 2/3, 0, 0,
                   0, 0, 0, 0, 0, 1/3, 2/3,
                   0, 0, 0, 0, 0, 1/2, 1/2,
                   0, 0, 0, 0, 0, 3/4, 1/4,
                   1/2, 1/2, 0, 0, 0, 0, 0,
                   1/4, 3/4, 0, 0, 0, 0, 0), byrow = TRUE, ncol = 7)
mcB <- new("markovchain", transitionMatrix = B)
period(mcB)
#> [1] 3
```
