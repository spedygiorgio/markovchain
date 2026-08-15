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

assessStationarity(sequence, nblocks, verbose = TRUE)

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

verifyMarkovProperty(
  sequence,
  method = c("G", "Pearson", "simulation"),
  B = 9999,
  seed = NULL,
  verbose = TRUE
)

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
#> Testing compatibility with a first-order Markov model
#>  Method: G 
#> Statistic: 0.5991613 
#> Degrees of freedom: 2 
#> p-value: 0.741129 
assessOrder(sequence)
#> Warning: The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.
#> Warning: Chi-squared approximation may be incorrect
#> Warning: Chi-squared approximation may be incorrect
#> The assessOrder test statistic is:  1.55307e-31 
#> The Chi-Square d.f. are:  2 
#> The p-value is:  1 
assessStationarity(sequence, 1)
#> Warning: The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.
#> Warning: Chi-squared approximation may be incorrect
#> Warning: Chi-squared approximation may be incorrect
#> The assessStationarity test statistic is:  0.1960191 
#> The Chi-Square d.f. are:  0 
#> The p-value is:  0 



#Example taken from Kullback Kupperman Tests for Contingency Tables and Markov Chains

sequence<-c(0,1,2,2,1,0,0,0,0,0,0,1,2,2,2,1,0,0,1,0,0,0,0,0,0,1,1,
2,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,2,1,0,
0,2,1,0,0,0,0,0,0,1,1,1,2,2,0,0,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,0,2,
0,1,1,0,0,0,1,2,2,0,0,0,0,0,0,2,2,2,1,1,1,1,0,1,1,1,1,0,0,2,1,1,
0,0,0,0,0,2,2,1,1,1,1,1,2,1,2,0,0,0,1,2,2,2,0,0,0,1,1)

mc=matrix(c(5/8,1/4,1/8,1/4,1/2,1/4,1/4,3/8,3/8),byrow=TRUE, nrow=3)
rownames(mc)<-colnames(mc)<-0:2; theoreticalMc<-as(mc, "markovchain")

verifyEmpiricalToTheoretical(data=sequence,object=theoreticalMc)
#> Testing compatibility of empirical and theoretical transition probabilities
#>  Method: G 
#> Statistic: 6.518384 
#> Degrees of freedom: 6 
#> p-value: 0.3676877 


data(kullback)
verifyHomogeneity(inputList=kullback,verbose=TRUE)
#> Testing homogeneity of transition probabilities across sequences
#>  Method: G 
#> Statistic: 160.9606 
#> Degrees of freedom: 29 
#> p-value: 3.087972e-20 
```
