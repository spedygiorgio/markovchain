# Stationary states of a `markovchain` object

This method returns the stationary vector in matricial form of a
markovchain object.

## Usage

``` r
steadyStates(object)
```

## Arguments

- object:

  A discrete `markovchain` object

## Value

A matrix corresponding to the stationary states

## Note

The steady states are identified starting from which eigenvectors
correspond to identity eigenvalues and then normalizing them to sum up
to unity. When negative values are found in the matrix, the eigenvalues
extraction is performed on the recurrent classes submatrix.

## References

A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall
2010

## See also

[`markovchain`](markovchain-class.md)

## Author

Giorgio Spedicato

## Examples

``` r
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, transitionMatrix =
                matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
                byrow = TRUE, dimnames=list(statesNames,statesNames)),
               name = "A markovchain Object" 
)       
steadyStates(markovB)
#>      a b c
#> [1,] 0 1 0
```
