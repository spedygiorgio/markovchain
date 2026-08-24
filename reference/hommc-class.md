# An S4 class for representing High Order Multivariate Markovchain (HOMMC)

An S4 class for representing High Order Multivariate Markovchain (HOMMC)

## Usage

``` r
hommc
```

## Slots

- `order`:

  an integer equal to order of Multivariate Markovchain

- `states`:

  a vector of states present in the HOMMC model

- `P`:

  array of transition matrices

- `Lambda`:

  a vector which stores the weightage of each transition matrices in P

- `byrow`:

  if FALSE each column sum of transition matrix is 1 else row sum = 1

- `name`:

  a name given to hommc

## Author

Giorgio Spedicato, Deepak Yadav

## Examples

``` r
statesName <- c("a", "b")

P <- array(0, dim = c(2, 2, 4), dimnames = list(statesName, statesName))
P[,,1] <- matrix(c(0, 1, 1/3, 2/3), byrow = FALSE, nrow = 2)
P[,,2] <- matrix(c(1/4, 3/4, 0, 1), byrow = FALSE, nrow = 2)
P[,,3] <- matrix(c(1, 0, 1/3, 2/3), byrow = FALSE, nrow = 2)
P[,,4] <- matrix(c(3/4, 1/4, 0, 1), byrow = FALSE, nrow = 2)

Lambda <- c(0.8, 0.2, 0.3, 0.7)

ob <- new("hommc", order = 1, states = statesName, P = P, 
          Lambda = Lambda, byrow = FALSE, name = "FOMMC")
          
```
