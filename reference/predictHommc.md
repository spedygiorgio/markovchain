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
homc. The last column of `init` holds the most recent state of each
sequence.

At each step the next state of sequence \\j\\ is drawn from \\\sum_k
\sum_h \lambda\_{jkh} P_h^{(jk)} x^{(k)}\_{t-h+1}\\, where
\\x^{(k)}\_{t-h+1}\\ is the state of sequence \\k\\ \\h - 1\\ steps
before the current one (Ching et al., 2008). The matrices are read by
column (`P[to, from]`), as returned by
[`fitHighOrderMultivarMC`](fitHighOrderMultivarMC.md), or by row when
`byrow = TRUE`, and all sequences are drawn from the same past before it
is updated.

## Author

Vandit Jain
