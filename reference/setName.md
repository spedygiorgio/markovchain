# Method to set name of markovchain object

This method modifies the existing name of markovchain object

## Usage

``` r
name(object) <- value

# S4 method for class 'markovchain'
name(object) <- value
```

## Arguments

- object:

  A markovchain object

- value:

  New name of markovchain object

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
name(markovB) <- "dangerous mc"
```
