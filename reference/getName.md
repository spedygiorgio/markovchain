# Method to retrieve name of markovchain object

This method returns the name of a markovchain object

## Usage

``` r
name(object)

# S4 method for class 'markovchain'
name(object)
```

## Arguments

- object:

  A markovchain object

## Author

Giorgio Spedicato, Deepak Yadav

## Examples

``` r
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, transitionMatrix =
                matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
                byrow = TRUE, dimnames=list(statesNames,statesNames)),
                name = "A markovchain Object" 
)
name(markovB)
#> [1] "A markovchain Object"
```
