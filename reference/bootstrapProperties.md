# Bootstrap confidence intervals for properties of an estimated chain

Estimates a Markov chain from one or more observed sequences and gives
bootstrap standard errors and confidence intervals for properties of the
chain, such as the second largest eigenvalue modulus, the spectral gap,
the relaxation time or the mixing time, which the corresponding
functions return as point values only.

## Usage

``` r
bootstrapProperties(
  data,
  statistics = c("slem", "spectralGap", "relaxationTime", "mixingTime"),
  nboot = 1000L,
  method = c("multinomial", "parametric"),
  level = 0.95,
  type = c("percentile", "basic"),
  possibleStates = character(),
  epsilon = 0.25,
  seed = NULL
)
```

## Arguments

- data:

  An observed sequence (character, factor or numeric vector), a list of
  independent sequences, or a square matrix of transition counts with
  the states as row and column names (only with
  `method = "multinomial"`).

- statistics:

  Either a character vector of built-in statistics, among `"slem"`,
  `"spectralGap"`, `"relaxationTime"`, `"mixingTime"` (the default
  four), `"entropyRate"`, `"kemenyConstant"` and `"steadyStates"`, or a
  named list whose elements are such names or functions of a
  `markovchain` object returning a numeric vector of fixed length (an
  unnamed built-in name is named after itself).

- nboot:

  A single integer of at least `2`: the number of bootstrap replicates.

- method:

  How the replicates are drawn: `"multinomial"` (the default) resamples,
  for each state, its observed number of departures from the estimated
  row; `"parametric"` simulates new sequences of the observed lengths
  from the estimated chain, each starting from the observed first state.
  See Details.

- level:

  The confidence level, a single number in \\(0, 1)\\.

- type:

  The kind of interval: `"percentile"` (the default) or `"basic"`.

- possibleStates:

  Optional character vector of states to include even if they are not
  observed. Every state must be left at least once in the data,
  otherwise its transition probabilities cannot be estimated.

- epsilon:

  The total variation threshold of [`mixingTime`](mixingTime.md).

- seed:

  An optional whole number. When given, the replicates are drawn after
  `set.seed(seed)` and the caller's random number stream is left
  untouched.

## Value

An object of class `"markovchainBootstrap"`, a list with:

- `estimate`:

  The statistics of the estimated chain.

- `ci`:

  A matrix with one row per statistic and the lower and upper limits of
  the interval.

- `se`, `bias`:

  Bootstrap standard errors and bias estimates (mean of the replicates
  minus the estimate).

- `failures`:

  For each statistic, the number of replicates for which it could not be
  computed (they are left out of the intervals).

- `replicates`:

  The `nboot` by statistics matrix of bootstrap values.

- `chain`:

  The estimated `markovchain`.

- `counts`:

  The matrix of observed transition counts.

- `method`, `type`, `level`, `nboot`:

  The settings used.

## Details

The chain is estimated by maximum likelihood, \\\hat P\_{ij} =
n\_{ij}/n_i\\. With `method = "multinomial"` a replicate draws, for
every state \\i\\, the counts \\(n^\*\_{i1}, \ldots, n^\*\_{ir})\\ from
a multinomial distribution with \\n_i\\ trials and probabilities \\\hat
P\_{i\cdot}\\, and re-estimates \\P^\*\_{ij} = n^\*\_{ij}/n_i\\. This is
the conditional bootstrap of Basawa et al. (1990), valid for finite
ergodic chains, which keeps the number of visits to each state fixed.
With `method = "parametric"` a replicate simulates a new path of each
observed length from \\\hat P\\, as in Kulperger and Prakasa Rao (1989);
the numbers of visits then vary too, and a replicate in which some state
is never left cannot be estimated and counts as a failure. The
parametric version is slower but closer to the actual sampling scheme of
a single path.

When all the statistics are among `"slem"`, `"spectralGap"`,
`"relaxationTime"` and `"mixingTime"`, the replicates are computed in
C++, with one eigen-decomposition per replicate; otherwise each
replicate is turned into a `markovchain` object and the statistics are
evaluated in R. The results are the same either way, and the random
draws come from R's generator, so `set.seed` and `seed` control them.

The percentile interval is formed by the \\(1-\code{level})/2\\ and
\\(1+\code{level})/2\\ quantiles of the replicates; the basic interval
reflects them around the estimate, \\2\hat\theta - q\_{1-\alpha/2}\\ and
\\2\hat\theta - q\_{\alpha/2}\\.

Some care is needed in reading the results. The SLEM and the spectral
gap are not smooth functions of the matrix where the second largest
eigenvalue modulus is reached by more than one eigenvalue, the mixing
time is an integer, and replicates of a chain with rare transitions may
be reducible or periodic, in which case some statistics fail; the
bootstrap distribution of all of them can be skewed. The intervals are
therefore approximate, and their coverage should be checked by
simulation for the case at hand, as done in the package's vignette. The
asymptotic theory of these estimators is developed by Hsu et al. (2019)
for the spectral gap of reversible chains.

## References

Basawa, I. V., Green, T. A., McCormick, W. P. and Taylor, R. L. (1990).
Asymptotic bootstrap validity for finite Markov chains. *Communications
in Statistics - Theory and Methods*, 19(4), 1493-1510.

Kulperger, R. J. and Prakasa Rao, B. L. S. (1989). Bootstrapping a
finite state Markov chain. *Sankhya, Series A*, 51(2), 178-191.

Hsu, D., Kontorovich, A., Levin, D. A., Peres, Y., Szepesvari, C. and
Wolfer, G. (2019). Mixing time estimation in reversible Markov chains
from a single sample path. *The Annals of Applied Probability*, 29(4),
2439-2480. [doi:10.1214/18-AAP1457](https://doi.org/10.1214/18-AAP1457)

## See also

[`slem`](slem.md), [`spectralGap`](spectralGap.md),
[`relaxationTime`](relaxationTime.md), [`mixingTime`](mixingTime.md),
[`markovchainFit`](markovchainFit.md) (whose `method = "bootstrap"`
gives bootstrap standard errors of the transition probabilities
themselves)

## Examples

``` r
data(rain)
set.seed(1)
boot <- bootstrapProperties(rain$rain, nboot = 200)
boot
#> Bootstrap of a Markov chain on 3 states (multinomial, 200 replicates)
#> 
#>                estimate      se   2.5%  97.5%
#> slem             0.4106 0.03039 0.3501 0.4710
#> spectralGap      0.5894 0.03039 0.5290 0.6499
#> relaxationTime   1.6970 0.08766 1.5390 1.8910
#> mixingTime       2.0000 0.12190 2.0000 2.0000
#> 
#> Percentile 95% intervals.
boot$ci
#>                     2.5%     97.5%
#> slem           0.3501461 0.4710463
#> spectralGap    0.5289537 0.6498539
#> relaxationTime 1.5388076 1.8905255
#> mixingTime     2.0000000 2.0000000

# any function of the chain, here the stationary probability of no rain
bootstrapProperties(rain$rain, nboot = 200, seed = 2,
  statistics = list(dry = function(mc) steadyStates(mc)[1, "0"]))
#> Bootstrap of a Markov chain on 3 states (multinomial, 200 replicates)
#> 
#>     estimate      se   2.5%  97.5%
#> dry   0.5009 0.02326 0.4528 0.5445
#> 
#> Percentile 95% intervals.
```
