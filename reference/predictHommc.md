# Simulate a higher order multivariate markovchain

This function provides a prediction of states for a higher order
multivariate markovchain object

## Usage

``` r
predictHommc(hommc,t,init)
```

## Arguments

- hommc:

  a hommc-class object

- t:

  no of iterations to predict

- init:

  matrix of previous states size of which depends on hommc

## Value

The function returns a matrix of size s X t displaying t predicted
states in each row coressponding to every categorical sequence.

## Details

The user is required to provide a matrix of giving n previous
coressponding every categorical sequence. Dimensions of the init are s X
n, where s is number of categorical sequences and n is order of the
homc.

## Author

Vandit Jain
