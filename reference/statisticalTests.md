# Various functions to perform statistical inference of DTMC

These functions verify the Markov property, assess the order and
stationarity of the Markov chain.

This function tests whether an empirical transition matrix is
statistically compatible with a theoretical one. It is a chi-square
based test. In case a cell in the empirical transition matrix is \>0
that is 0 in the theoretical transition matrix the null hypothesis is
rejected.

Verifies that the s elements in the input list belongs to the same DTMC

Tests the null hypothesis that the conditional distribution of the next
state depends only on the current state, against a second-order Markov
alternative. Equivalently, it tests conditional independence of the
previous and next states given the current state.

`assessStationarity` tests whether transition probabilities are constant
across the specified blocks. Structural zeros can be supplied explicitly
with `structural.zeros`; sampling zeros remain part of the contingency
table.

## Usage

``` r
verifyMarkovProperty(
  sequence,
  method = c("G", "Pearson", "simulation"),
  B = 9999,
  seed = NULL,
  verbose = TRUE
)

assessOrder(sequence, verbose = TRUE)

assessStationarity(sequence, nblocks, structural.zeros = NULL, verbose = TRUE)

verifyEmpiricalToTheoretical(
  data,
  object,
  method = c("G", "Pearson", "simulation"),
  B = 9999,
  seed = NULL,
  verbose = TRUE
)

verifyHomogeneity(
  inputList,
  method = c("G", "Pearson", "simulation"),
  B = 9999,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- sequence:

  An empirical sequence of states.

- method:

  Test statistic: \`"G"\`, \`"Pearson"\`, or \`"simulation"\`.

- B:

  Number of Monte Carlo replicates for \`method = "simulation"\`.

- seed:

  Optional random seed.

- verbose:

  Should results be printed?

- nblocks:

  Number of blocks.

- structural.zeros:

  Optional logical transition matrix. \`TRUE\` marks a structurally
  impossible transition. Sampling zeros are not removed.

- data:

  An empirical sequence or a matrix of transition counts.

- object:

  A \`markovchain\` object specifying theoretical probabilities.

- inputList:

  A list whose elements are empirical sequences or matrices.

## Value

Verification result

a list with following slots: statistic (the chi - square statistic), dof
(degrees of freedom), and corresponding p-value. In case a cell in the
empirical transition matrix is \>0 that is 0 in the theoretical
transition matrix the null hypothesis is rejected. In that case a
p-value of 0 and statistic and dof of NA are returned.

a list of transition matrices?

A list with statistic, degrees of freedom, p-value, observed and
expected triple-transition counts, and the fitted first-order transition
matrix. For simulation, simulated statistics and \`B\` are also
returned.

A list containing statistic, degrees of freedom, p-value, method,
observed counts, expected counts, and theoretical transition matrix.

A list containing statistic, degrees of freedom, p-value, method, pooled
transition counts, and individual transition counts.

## References

Anderson and Goodman.

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. \*The Annals of Mathematical Statistics\*, 28(1), 89–110.

Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for Contingency
Tables and Markov Chains. \*Technometrics\*, 4(4), 573–608.

Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for Contingency
Tables and Markov Chains. \*Technometrics\*, 4(4), 573–608.

## See also

`markovchain`

## Author

Tae Seung Kang, Giorgio Alfredo Spedicato

## Examples

``` r
sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b",
              "a", "b", "a", "a", "b", "b", "b", "a")
mcFit <- markovchainFit(data = sequence, byrow = FALSE)
verifyMarkovProperty(sequence)
#> 
#>  Likelihood-ratio test
#> 
#> data:  sequence
#> G-squared = 0.59916, df = 2, p-value = 0.7411
#> 
assessOrder(sequence)
#> 
#>  Pearson's Chi-squared test for Markov order
#> 
#> data:  sequence
#> X-squared = 0.63, df = 2, p-value = 0.7298
#> 
assessStationarity(sequence, 2)
#> 
#>  Pearson's Chi-squared test for time-homogeneity
#> 
#> data:  sequence
#> X-squared = 1.345, df = 2, p-value = 0.5104
#> 

# Structural zeros: a -> a is impossible.
structural.zeros <- matrix(FALSE, 2, 2,
                           dimnames = list(c("a", "b"), c("a", "b")))
structural.zeros["a", "a"] <- TRUE
assessStationarity(c("a", "b", "a", "b", "a", "b", "b", "a"),
                   nblocks = 2, structural.zeros = structural.zeros)
#> 
#>  Pearson's Chi-squared test for time-homogeneity
#> 
#> data:  c("a", "b", "a", "b", "a", "b", "b", "a")
#> X-squared = 1.3333, df = 1, p-value = 0.2482
#> 



#Example taken from Kullback Kupperman Tests for Contingency Tables and Markov Chains

sequence<-c(0,1,2,2,1,0,0,0,0,0,0,1,2,2,2,1,0,0,1,0,0,0,0,0,0,1,1,
2,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,2,1,0,
0,2,1,0,0,0,0,0,0,1,1,1,2,2,0,0,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,0,2,
0,1,1,0,0,0,1,2,2,0,0,0,0,0,0,2,2,2,1,1,1,1,0,1,1,1,1,0,0,2,1,1,
0,0,0,0,0,2,2,1,1,1,1,1,2,1,2,0,0,0,1,2,2,2,0,0,0,1,1)

mc=matrix(c(5/8,1/4,1/8,1/4,1/2,1/4,1/4,3/8,3/8),byrow=TRUE, nrow=3)
rownames(mc)<-colnames(mc)<-0:2; theoreticalMc<-as(mc, "markovchain")

verifyEmpiricalToTheoretical(data=sequence,object=theoreticalMc)
#> 
#>  Likelihood-ratio test
#> 
#> data:  sequence
#> G-squared = 6.5184, df = 6, p-value = 0.3677
#> 


data(kullback)
verifyHomogeneity(inputList=kullback,verbose=TRUE)
#> 
#>  Likelihood-ratio test
#> 
#> data:  kullback
#> G-squared = 160.96, df = 29, p-value < 2.2e-16
#> 
```
