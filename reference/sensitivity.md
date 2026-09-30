# Sensitivity of the stationary distribution to a state's transition row

Computes, for a finite, irreducible discrete-time Markov chain, the
coefficients needed to obtain the first-order change in the stationary
distribution caused by an infinitesimal perturbation of one row of the
transition matrix.

## Usage

``` r
sensitivity(object, state)

# S4 method for class 'markovchain'
sensitivity(object, state)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

- state:

  A single state name (character) or state index (single positive
  integer), identifying the row \\k\\ of the transition matrix to be
  perturbed.

## Value

An \\n\times n\\ numeric matrix \\S\\, with both dimensions named after
`states(object)`. Row \\l\\ of \\S\\ corresponds to the perturbation
direction "increase \\p\_{kl}\\" (paired with a compensating decrease
elsewhere in row \\k\\); column \\j\\ corresponds to the affected
stationary probability \\\pi_j\\. See Details for how to read individual
entries.

## Details

A single entry \\p\_{kl}\\ of a stochastic matrix cannot be perturbed on
its own without leaving row \\k\\: some other entry (or entries) of that
row must move to compensate, so that the row still sums to one. Any
admissible perturbation of row \\k\\ is therefore a direction vector
\\d\in\mathbb{R}^n\\ with \\\sum_l d_l = 0\\, giving the perturbed
matrix \\P(\varepsilon) = P + \varepsilon\\ e_k d^{\mathsf T}\\ for
small \\\varepsilon\\. This function returns the \\n\times n\\ matrix
\\S\\ such that, for *every* such \\d\\ and every state \\j\\,
\$\$\left.\frac{d\pi_j}{d\varepsilon}\right\|\_{\varepsilon=0} = \sum_l
d_l\\ S\_{lj} = \left(d^{\mathsf T} S\right)\_j.\$\$

**Closed form.** Let \\Z=(I-P+\mathbf 1\pi^{\mathsf T})^{-1}\\ be the
fundamental matrix already used by
[`kemenyConstant`](kemenyConstant.md). Then \$\$S\_{lj} =
\pi_k\left(Z\_{lj} - \pi_j\right).\$\$ This particular centering
(subtracting \\\pi_j\\, the same constant for every row \\l\\) is what
makes \\S\\ usable directly with *any* zero-sum direction \\d\\, because
\\\sum_l d_l \pi_j = \pi_j\sum_l d_l = 0\\ drops out of the sum above –
adding any other per-column constant to \\S\\ would give the same
directional derivatives, but this one has the convenient side effect
that `sensitivity(object, state)[state, ]` is the sensitivity of
"leaving row `state` unchanged", which is informative on its own (it
need not be zero: the \*direction\* \\d=e\_{\mathrm{state}}\\ is
generally not itself a valid zero-sum perturbation by itself, only
differences of rows are).

**The common two-state case.** The usual textbook question – "increase
\\p\_{k,\mathrm{to}}\\ by \\\varepsilon\\, decrease
\\p\_{k,\mathrm{from}}\\ by \\\varepsilon\\, how does \\\pi\\ move?" –
is answered by taking the difference of two rows of \\S\\:
\$\$\left.\frac{d\pi}{d\varepsilon}\right\|\_{\varepsilon=0} =
S\_{\mathrm{to}, \cdot} - S\_{\mathrm{from}, \cdot}.\$\$ See the second
example below, which checks this against a direct finite-difference
recomputation of the stationary distribution.

Only irreducibility is required, not aperiodicity: \\Z\\ and \\\pi\\ are
well defined for any irreducible chain regardless of periodicity.

The implementation calls [`steadyStates`](steadyStates.md) once and then
solves one dense linear system for \\Z\\; both are \\O(n^3)\\ time and
\\O(n^2)\\ memory for a dense \\n\\-state transition matrix, the same
cost as [`kemenyConstant`](kemenyConstant.md). It supports both row- and
column-stochastic storage; \\S\\ is always returned with rows/columns
indexed by state name in the chain's own state order.

## References

Schweitzer, P. J. (1968). Perturbation theory and finite Markov chains.
*Journal of Applied Probability*, 5(2), 401-413.

Meyer, C. D. (1980). The condition of a finite Markov chain and
perturbation bounds for the limiting probabilities. *SIAM Journal on
Algebraic and Discrete Methods*, 1(3), 273-283.

Cho, G. E. and Meyer, C. D. (2001). Comparison of perturbation bounds
for the stationary distribution of a Markov chain. *Linear Algebra and
its Applications*, 335(1-3), 137-150.

## See also

[`kemenyConstant`](kemenyConstant.md),
[`steadyStates`](steadyStates.md), [`is.irreducible`](is.irreducible.md)

## Examples

``` r
statesNames <- c("a", "b", "c")
mc <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0.5, 0.3, 0.2,
                              0.2, 0.6, 0.2,
                              0.1, 0.1, 0.8), byrow = TRUE, nrow = 3,
                            dimnames = list(statesNames, statesNames)))
S <- sensitivity(mc, "a")
S
#>             a           b          c
#> a  0.27879009 -0.01093294 -0.2678571
#> b -0.02733236  0.29518950 -0.2678571
#> c -0.10386297 -0.16399417  0.2678571

# Check against a finite-difference recomputation of the stationary
# distribution: increase p("a"->"c") and decrease p("a"->"b") by eps.
eps <- 1e-6
P2 <- mc@transitionMatrix
P2["a", "c"] <- P2["a", "c"] + eps
P2["a", "b"] <- P2["a", "b"] - eps
mc2 <- new("markovchain", states = statesNames, transitionMatrix = P2)
(steadyStates(mc2) - steadyStates(mc)) / eps   # finite difference
#>                a          b         c
#> [1,] -0.07653058 -0.4591835 0.5357141
S["c", ] - S["b", ]                            # closed-form prediction
#>           a           b           c 
#> -0.07653061 -0.45918367  0.53571429 
```
