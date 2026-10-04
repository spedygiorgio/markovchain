# Topological entropy of a Markov chain

Computes the topological entropy of the graph of a discrete-time Markov
chain: the exponential growth rate of the number of distinct admissible
paths, ignoring their probabilities.

## Usage

``` r
topologicalEntropy(object, base = 2)

# S4 method for class 'markovchain'
topologicalEntropy(object, base = 2)
```

## Arguments

- object:

  A `markovchain` object.

- base:

  A finite numeric scalar strictly greater than one. The default, `2`,
  returns bits, consistently with [`entropyRate`](entropyRate.md).

## Value

A non-negative numeric scalar in units determined by `base`.

## Details

If \\A\\ is the 0/1 adjacency matrix with \\A\_{ij}=1\\ exactly when
\\p\_{ij}\>0\\, the topological entropy is \$\$h\_{top} = \log_b
\rho(A),\$\$ where \\\rho(A)\\ is the spectral radius (Perron root) of
\\A\\.

Only the pattern of positive entries matters: the value depends on which
transitions are possible, not on how likely they are. It is the upper
bound of the entropy rate over all the Markov chains sharing that graph
(the variational principle, see Parry, 1964), so
`entropyRate(object) <= topologicalEntropy(object)` for an irreducible
chain. The bound is attained by the maximal-entropy (Parry) chain on the
same graph, and also, for instance, by a chain whose every row is
uniform over a common number of successors. A chain that is a single
cycle (deterministic dynamics) has `topologicalEntropy = 0`.

No irreducibility is needed: for a reducible chain the result is the
largest value over its communicating classes. A probability that is
positive but numerically tiny counts as a transition, exactly as in
[`is.irreducible`](is.irreducible.md).

The cost is one eigenvalue computation, \\O(n^3)\\ time and \\O(n^2)\\
memory for a dense chain. It mirrors PyDTMC's `topological_entropy`,
which uses the natural logarithm.

## References

Parry, W. (1964). Intrinsic Markov chains. *Transactions of the American
Mathematical Society*, 112, 55-66.

Cover, T. M. and Thomas, J. A. (2006). *Elements of Information Theory*,
2nd edition. Wiley.

## See also

[`entropyRate`](entropyRate.md),
[`normalizedEntropyRate`](normalizedEntropyRate.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
topologicalEntropy(mc)
#> [1] 1
```
