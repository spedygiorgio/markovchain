# Build a population-genetics Markov chain (Moran or Wright-Fisher)

Constructs the Markov chain tracking the number of copies of a mutant
allele in a population of `n` individuals under mutation and viability
selection, for either the Moran or the Wright-Fisher model of population
genetics.

## Usage

``` r
populationGeneticsModel(
  model = c("moran", "wright-fisher"),
  n,
  s = 0,
  u = 1e-09,
  v = 1e-09,
  states = NULL
)
```

## Arguments

- model:

  Either `"moran"` or `"wright-fisher"`, selecting which reproduction
  scheme generates the transition matrix.

- n:

  A single integer of at least `2`: the (constant) population size. The
  chain has `n + 1` states, \\0,1,\ldots,n\\, the possible mutant-allele
  counts.

- s:

  A single finite number greater than `-1`: the selection coefficient.
  The mutant allele's fitness relative to the wild-type is \\1+s\\
  (`s = 0` is neutral drift, `s > 0` favours the mutant, `-1 < s < 0`
  disfavours it).

- u:

  A single number in \\\[0,1\]\\: the backward mutation rate, mutant to
  wild-type. Defaults to `1e-9`, matching the effectively mutation-free
  chain conventionally used for the neutral/absorbing case.

- v:

  A single number in \\\[0,1\]\\: the forward mutation rate, wild-type
  to mutant. Defaults to `1e-9`.

- states:

  An optional character vector of `n + 1` state names. Defaults to
  `as.character(0:n)`.

## Value

A new, row-stochastic `markovchain` object on `n + 1` states. States
`"0"` and `as.character(n)` (loss and fixation of the mutant allele) are
absorbing, following the standard textbook treatment of both models;
every interior state has a transition row determined by `model`, `s`,
`u` and `v` as detailed below.

## Details

**Moran model.** At each step one individual is chosen to reproduce
(with probability proportional to its type's relative fitness, mutant
fitness \\r=1+s\\ against wild-type fitness \\1\\) and one individual,
chosen uniformly at random, dies and is replaced by the offspring, which
mutates with probability `u` or `v` depending on the parent's type. From
interior state \\i\\ (\\0\<i\<n\\), writing \\r_i=(1+s)i\\ and
\\m_i=n-i\\, \$\$P\_{i,i-1}=\frac{i}{n}\cdot\frac{r_i v +
m_i(1-u)}{r_i+m_i},\qquad P\_{i,i+1}=\frac{m_i}{n}\cdot\frac{r_i(1-u) +
m_i v}{r_i+m_i},\$\$ with \\P\_{ii}\\ the remainder.

**Wright-Fisher model.** Generations are discrete and non-overlapping:
from interior state \\i\\, the mutant-allele frequency \\k=i/n\\ is
first updated for mutation, \\p=k(1-u)+(1-k)v\\, then for viability
selection, \\p'=\min\bigl(p(1+s)/(1+ps),\\1\bigr)\\, and the next
generation's mutant count is \\\mathrm{Binomial}(n,p')\\-distributed:
\\P\_{ij}=\binom{n}{j}p'^{\\j}(1-p')^{n-j}\\.

Both constructions follow PyDTMC's `population_genetics_model`, with one
deliberate difference: this function uses mutant relative fitness
\\1+s\\ in *both* models (PyDTMC's Moran implementation instead uses
\\1-s\\, so that a positive `s` there favours the wild-type rather than
the mutant, the opposite of its own Wright-Fisher convention and of the
usual textbook one). `s = 0` (neutral drift) and the mutation-rate terms
are unaffected by this choice and match PyDTMC exactly; away from
`s = 0` the two packages' Moran chains differ by construction, by
design, to keep `s`'s meaning consistent between `"moran"` and
`"wright-fisher"` within this package. As `u, v` shrink to `0`, both
reduce to the classical drift-only chains with absorbing loss/fixation
states, going back to Wright (1931) and Moran (1958); see Ewens (2004)
for a modern textbook treatment of both.

## References

Wright, S. (1931). Evolution in Mendelian populations. *Genetics*,
16(2), 97-159.

Moran, P. A. P. (1958). Random processes in genetics. *Mathematical
Proceedings of the Cambridge Philosophical Society*, 54(1), 60-71.

Ewens, W. J. (2004). *Mathematical Population Genetics I: Theoretical
Introduction* (2nd ed.). Springer.

## See also

[`birthDeath`](birthDeath.md), [`gamblersRuin`](gamblersRuin.md)

## Examples

``` r
# Neutral drift (s = 0): absorption ("fixation") probability of the
# mutant allele from state i equals i/n, the classical result.
neutral <- populationGeneticsModel(model = "moran", n = 8, s = 0)
ap <- absorptionProbabilities(neutral)
ap["4", "8"] # close to 4/8 = 0.5
#> [1] 0.5

# Positive selection increases the fixation probability from any
# interior starting count.
favoured <- populationGeneticsModel(model = "wright-fisher", n = 8, s = 0.5)
absorptionProbabilities(favoured)["4", "8"]
#> [1] 0.9639052
```
