# Higher order Markov chains

## Higher Order Markov Chains

Continuous time Markov chains are discussed in the CTMC vignette which
is a part of the package.

An experimental `fitHigherOrder` function has been written in order to
fit a higher order Markov chain (Ching et al.
([2008](#ref-ching2008higher))). `fitHigherOrder` takes three inputs

1.  sequence: a categorical data sequence.
2.  order: order of Markov chain to fit with default value 2.
3.  method: how the weights are estimated, `"lsq"` (default) or `"mle"`,
    see below.

The output will be a `list` which consists of

1.  lambda: model parameter(s).
2.  Q: a list of transition matrices. $Q_{i}$ is the $ith$ step
    transition matrix stored column-wise.
3.  X: frequency probability vector of the given sequence.

The model is a mixture of the empirical lag-$i$ transition matrices
$Q_{i}$ with weights $\lambda_{i} \geq 0$ summing to one, and the
matrices $Q_{i}$ are the same for both methods. With the default
`method = "lsq"` the weights solve a quadratic programming problem, the
minimization of the squared distance between the stationary distribution
and its image under the mixture, which is solved using `solnp` function
of the Rsolnp package ([Ghalanos and Theussl 2014](#ref-pkg:Rsolnp)).
With `method = "mle"` the weights maximize the log-likelihood of the
observations; for fixed $Q_{i}$ this is a concave problem, so its
maximum is global, and it is solved by the EM algorithm for mixture
weights, without the Rsolnp package.

``` r
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  data(rain)
  rain_small <- rain$rain[1:150]
  fitHigherOrder(rain_small, 2)
}
#> $lambda
#> [1] 0.77333 0.22667
#> 
#> $Q
#> $Q[[1]]
#>             0       1-5      6+
#> 0   0.5538462 0.4230769 0.21875
#> 1-5 0.2923077 0.4038462 0.37500
#> 6+  0.1538462 0.1730769 0.40625
#> 
#> $Q[[2]]
#>             0       1-5      6+
#> 0   0.5384615 0.3725490 0.34375
#> 1-5 0.2615385 0.3921569 0.43750
#> 6+  0.2000000 0.2352941 0.21875
#> 
#> 
#> $X
#>         0       1-5        6+ 
#> 0.4333333 0.3466667 0.2200000
```

## Comparing models of different orders

`higherOrderLogLik` evaluates the log-likelihood of a sequence under the
model returned by `fitHigherOrder`, and derives the deviance ($-2$ times
the log-likelihood), the AIC and the BIC, so that models of different
orders can be compared. The probability of moving to state $x_{t}$ is
the mixture
$\sum_{i = 1}^{k}\lambda_{i}Q_{i}\lbrack x_{t},x_{t - i}\rbrack$, the
form of the mixture transition distribution model of ([Raftery
1985](#ref-raftery1985model)), and the log-likelihood is the sum of the
logarithms of these probabilities over the observations that a model of
order $k$ can predict, that is from observation $k + 1$ onwards. To
compare orders on exactly the same observations the argument `start`
must be set to one plus the largest order compared.

Two caveats apply. First, with the default `method = "lsq"` the weights
$\lambda$ are chosen by least squares on the stationary distribution,
not by maximum likelihood, so the value returned is the log-likelihood
*of the fitted model* and not the maximum attainable one. It can
therefore be lower for a higher order than for a lower one, which cannot
happen for maximum-likelihood fits of nested models. With
`method = "mle"` the weights maximize this log-likelihood for the
observations that a model of that order can predict, and in the examples
below the log-likelihood then never decreases with the order. Second,
the number of parameters used for the criteria is
$(r - 1)\,(1 + k(r - 1))$, with $r$ the number of states. Each lag
matrix has $r(r - 1)$ free probabilities and there are $k - 1$ free
weights, but in a mixture of lag matrices the weights are not
identifiable (a distribution common to all departure states can be moved
from one lag to another without changing any transition probability), so
the model can represent a set of transition laws of dimension
$(r - 1)(1 + k(r - 1))$, that is $r(k - 1)$ less than the naive count
$k\, r(r - 1) + (k - 1)$; the two coincide for $k = 1$.

The example compares orders one to three on the Alofi Island daily
rainfall and on the preproglucacon DNA sequence, both analysed by ([P.
J. Avery and D. A. Henderson 1999](#ref-averyHenderson)), with both
estimation methods. Maximum likelihood weights raise the log-likelihood
of the higher orders, by up to about 11 units for the rainfall and 28
for the DNA sequence, whereas the least squares weights can leave it
below that of the first-order model. Even so, the gain is too small to
pay for the additional parameters, and both criteria select the
first-order model for both sequences with both methods. The values are
computed by the package; they are not claimed to reproduce those of the
original paper.

``` r
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  compareOrders <- function(sequence, orders = 1:3) {
    byMethod <- lapply(c("lsq", "mle"), function(method) {
      fits <- lapply(orders, function(k) fitHigherOrder(sequence, k, method = method))
      sapply(fits, function(f)
        unlist(higherOrderLogLik(sequence, f, start = max(orders) + 1)[c("logLik", "BIC")]))
    })
    out <- rbind(byMethod[[1]], byMethod[[2]])
    rownames(out) <- c("logLik (lsq)", "BIC (lsq)", "logLik (mle)", "BIC (mle)")
    colnames(out) <- paste("order", orders)
    round(out, 1)
  }
  data(rain)
  print(compareOrders(rain$rain))
  data(preproglucacon)
  print(compareOrders(preproglucacon$preproglucacon))
}
#>              order 1 order 2 order 3
#> logLik (lsq) -1038.1 -1047.8 -1047.8
#> BIC (lsq)     2118.1  2165.5  2193.5
#> logLik (mle) -1038.1 -1036.5 -1036.5
#> BIC (mle)     2118.1  2142.9  2170.9
#>              order 1 order 2 order 3
#> logLik (lsq) -2026.0 -2024.6 -2051.6
#> BIC (lsq)     4140.3  4203.7  4323.9
#> logLik (mle) -2026.0 -2024.0 -2024.0
#> BIC (mle)     4140.3  4202.5  4268.7
```

### Selecting the order of a full Markov chain

The models above restrict the dependence on the past to a mixture of
lags. `selectOrder` instead fits the fully parameterized Markov chains
of order $0,1,\ldots,K$ (order 0 is independence) by maximum likelihood,
all on the same observations (from `start`, by default $K + 1$), and
selects the order that minimizes the BIC or the AIC ([Tong
1975](#ref-tong1975determination); [Katz 1981](#ref-katz1981some)). An
order-$k$ chain on $r$ states has $r^{k}(r - 1)$ free parameters; with
`parameters = "observed"` only the transition probabilities not
estimated as zero are counted, as in Berchtold and Raftery
([2002](#ref-berchtold2002mixture)). The table also gives the
likelihood-ratio statistic of each order against the previous one, with
its asymptotic chi-squared p-value ([Anderson and Goodman
1957](#ref-anderson1957statistical)). BIC is a consistent estimator of
the order ([Csiszár and Shields 2000](#ref-csiszar2000consistency)),
whereas AIC tends to choose higher orders, and both become unreliable
when $r^{K}$ is not small compared with the number of observations.

``` r
data(rain)
rainOrder <- selectOrder(rain$rain, maxOrder = 3)
rainOrder$order
#> [1] 1
print(rainOrder$table, digits = 4, row.names = FALSE)
#>  order logLik npar  AIC  BIC     LR df   p.value
#>      0  -1134    2 2272 2282     NA NA        NA
#>      1  -1038    6 2088 2118 191.53  4 2.487e-40
#>      2  -1025   18 2086 2176  25.94 12 1.096e-02
#>      3  -1006   54 2119 2389  39.06 36 3.338e-01
selectOrder(rain$rain, maxOrder = 3, criterion = "AIC")$order
#> [1] 2
```

For the Alofi rainfall the BIC selects the first-order chain, whereas
the AIC slightly prefers the second-order one (2086.2 against 2088.1),
which needs 18 parameters instead of 6. For the preproglucacon sequence
(`selectOrder(preproglucacon$preproglucacon, 3)`) the BIC of
independence and of the first-order chain differ by about one unit, and
the AIC selects order one. These values are computed by the package.

### Reproducing a published comparison

([Berchtold and Raftery 2002](#ref-berchtold2002mixture)) compare, by
log-likelihood and BIC, independence, Markov chains of order one to
three and mixture transition distribution (MTD) models for the hourly
wind direction at Koeberg (South Africa; 744 observations recoded into
four directions, originally from ([MacDonald and Zucchini
1997](#ref-macdonald1997hidden))) and for a daily series of epileptic
seizures (204 observations). The authors kindly provided the two series,
which are distributed with the package in `inst/extdata`. Their
convention is to condition every model on the first 14 observations, so
that all models are evaluated on the same $n - 14$ observations, which
in `higherOrderLogLik` corresponds to `start = 15`; the BIC uses
$n - 14$ as sample size and counts only the parameters that are not
forced to zero.

`fitHigherOrder` does not estimate the MTD model of that paper: it fits
a different transition matrix for each lag, with weights chosen by least
squares or, with `method = "mle"`, by maximum likelihood given those
matrices, whereas the MTD model uses a single matrix $Q$ for all lags
and estimates it together with the weights by maximum likelihood (that
model is fitted by `fitMTD`, described in the next section). The
comparison below evaluates, with `higherOrderLogLik`, the first-order
chain estimated on the observations entering the likelihood and the
MTD(2) model with the weights and the matrix $Q$ printed in the paper
(whose rows are the departure states, hence the transposition).

``` r
koeberg <- as.character(read.csv(system.file("extdata", "koeberg_wind.csv",
                                             package = "markovchain"))$state)
start <- 15
n_eff <- length(koeberg) - (start - 1)

# first-order Markov chain, estimated on the transitions that enter the likelihood
Q1 <- seq2matHigh(koeberg[(start - 1):length(koeberg)], 1)
mc1 <- higherOrderLogLik(koeberg, list(lambda = 1, Q = list(Q1)), start = start)

# MTD(2) with lambda and Q as printed in Section 1.3 of the paper
Qpaper <- matrix(c(0.8301, 0.0689, 0.0077, 0.0933,
                   0.0369, 0.9012, 0.0619, 0.0000,
                   0.0155, 0.1553, 0.8070, 0.0222,
                   0.0779, 0.0000, 0.0528, 0.8693), 4, 4, byrow = TRUE)
Q <- t(Qpaper)
dimnames(Q) <- list(as.character(1:4), as.character(1:4))
mtd2 <- higherOrderLogLik(koeberg, list(lambda = c(0.7569, 0.2431), Q = list(Q, Q)),
                          start = start)

# parameters not forced to zero: 11 for the chain (one empty transition),
# 4 * 3 - 2 + (2 - 1) = 11 for the MTD(2) (two structural zeros in Q)
comparison <- data.frame(
  model = c("Markov chain, order 1", "MTD, order 2"),
  logLik = c(mc1$logLik, mtd2$logLik),
  BIC = c(-2 * mc1$logLik + 11 * log(n_eff), -2 * mtd2$logLik + 11 * log(n_eff)),
  logLik_published = c(-413.3, -393.4),
  BIC_published = c(899.1, 859.3))
print(comparison, digits = 4, row.names = FALSE)
#>                  model logLik   BIC logLik_published BIC_published
#>  Markov chain, order 1 -413.3 899.1           -413.3         899.1
#>           MTD, order 2 -393.4 859.3           -393.4         859.3
```

The values coincide with Table 2 of the paper up to its rounding to one
decimal, and the BIC prefers the MTD(2) model to the first-order chain,
as in the paper. The same agreement is obtained for the Markov chains of
order two and three and for the seizure series (Table 3); these checks
are part of the unit tests of the package. The rows for independence and
for the Markov chains of order one to three are obtained in one call
with
`selectOrder(koeberg, maxOrder = 3, start = 15, parameters = "observed")`:

``` r
print(selectOrder(koeberg, maxOrder = 3, start = 15, parameters = "observed")$table[, 1:5],
      digits = 5, row.names = FALSE)
#>  order  logLik npar     AIC     BIC
#>      0 -954.83    3 1915.66 1929.44
#>      1 -413.30   11  848.61  899.13
#>      2 -374.92   27  803.84  927.85
#>      3 -346.19   39  770.39  949.52
```

## The mixture transition distribution model

A Markov chain of order $k$ on $r$ states has $r^{k}(r - 1)$ free
transition probabilities, a number that grows so quickly with $k$ that
high orders can rarely be estimated. The mixture transition distribution
(MTD) model of ([Raftery 1985](#ref-raftery1985model)) replaces the full
transition array with a mixture of contributions of the individual lags,
all governed by the same transition matrix:
$$P(X_{t} = j \mid X_{t - 1} = i_{1},\ldots,X_{t - k} = i_{k}) = \sum\limits_{g = 1}^{k}\lambda_{g}\, q_{i_{g}j},$$
where $Q = (q_{ij})$ is an $r \times r$ transition matrix (rows are
departure states) and the lag weights $\lambda_{g}$ sum to one. The
model has only $r(r - 1) + k - 1$ parameters, one more for each
additional lag, and for $k = 1$ it is the first-order Markov chain.
([Berchtold and Raftery 2002](#ref-berchtold2002mixture)) review the
model, its extensions and its applications.

`fitMTD(sequence, order, start, nstart, tol, maxit)` estimates $Q$ and
$\lambda$ by maximum likelihood. As in most applications, the weights
are constrained to be non-negative (Raftery’s original formulation also
admits negative weights, provided that every transition probability
stays in $\lbrack 0,1\rbrack$; this case is not supported). Under this
constraint the model is a mixture in which an unobserved lag generates
each observation, and the likelihood is maximized by the EM algorithm of
([Lèbre and Bourguignon 2008](#ref-lebre2008em)), implemented in C++:
the E-step computes the posterior probability of each lag for each
observation, the M-step updates the weights as the average of these
probabilities and $Q$ as the transition counts weighted by them. Each
iteration increases the likelihood, but the MTD likelihood can have
several local maxima ([Berchtold 2001](#ref-berchtold2001estimation)),
so the models of order $1,\ldots,k$ are fitted in turn on the same
observations and each order is started both from equal weights and from
the fit of the previous order, extended with a zero and with a small
positive weight for the new lag. The first extension has the likelihood
of the previous order, so the likelihood returned never decreases with
the order, as it must for nested models; `nstart - 1` further random
starting points can be added, and the best fit is returned.

The likelihood is conditional on the observations before `start`, by
default `order + 1`; as for `higherOrderLogLik`, models of different
orders are comparable only when they share `start`. The function returns
the weights, the matrix $Q$ as a `markovchain` object (`estimate`), the
maximized log-likelihood with AIC and BIC (based on $r(r - 1) + k - 1$
parameters), and, in the element `Q`, the matrix in the column layout
used by `fitHigherOrder`, so that the fit can also be passed to
`higherOrderLogLik`.

The chunk below estimates the MTD models of order two and three on both
series of ([Berchtold and Raftery 2002](#ref-berchtold2002mixture)),
with their conventions (`start = 15`, and the elements of $Q$ estimated
as zero excluded from the number of parameters of the BIC), and compares
the results with those published in Tables 2 and 3 of the paper.

``` r
readSeries <- function(file, column)
  read.csv(system.file("extdata", file, package = "markovchain"))[[column]]
series <- list(Koeberg = readSeries("koeberg_wind.csv", "state"),
               seizures = readSeries("epileptic_seizures.csv", "seizure"))
published <- data.frame(series = rep(c("Koeberg", "seizures"), each = 2),
                        order = c(2, 3, 2, 3),
                        logLik_published = c(-393.4, -393.2, -119.5, -117.7),
                        BIC_published = c(859.3, 865.6, 254.7, 256.4))
fits <- Map(function(s, k) fitMTD(series[[s]], order = k, start = 15),
            published$series, published$order)
names(fits) <- paste(published$series, published$order)
estimated <- t(sapply(fits, function(fit) {
  zeros <- sum(fit$estimate@transitionMatrix < 1e-8)
  c(logLik = fit$logLikelihood,
    BIC = -2 * fit$logLikelihood + (fit$npar - zeros) * log(fit$nobs))
}))
print(cbind(published, round(estimated, 1)), row.names = FALSE)
#>    series order logLik_published BIC_published logLik   BIC
#>   Koeberg     2           -393.4         859.3 -393.4 859.3
#>   Koeberg     3           -393.2         865.6 -393.2 865.6
#>  seizures     2           -119.5         254.7 -119.5 254.7
#>  seizures     3           -117.7         256.4 -117.7 256.4
```

All the published values are reproduced. For the wind series the BIC
selects the MTD(2) model, with 11 parameters, over the first-order chain
and over the Markov chains of order two and three, which need up to 39
parameters (see the previous section and Table 2 of the paper). The
estimated weights and transition matrix of the MTD(2) model agree with
those printed in Section 1.3 of the paper to the third decimal:

``` r
fits[["Koeberg 2"]]$lambda
#>      lag1      lag2 
#> 0.7568342 0.2431658
round(fits[["Koeberg 2"]]$estimate@transitionMatrix, 4)
#>        1      2      3      4
#> 1 0.8307 0.0687 0.0077 0.0930
#> 2 0.0370 0.9011 0.0619 0.0000
#> 3 0.0154 0.1552 0.8071 0.0223
#> 4 0.0779 0.0000 0.0527 0.8694
```

The weight of the first lag is about three times that of the second: the
wind direction depends mostly on the previous hour, but the hour before
still adds information, at the cost of a single extra parameter. More
general MTD variants (different matrices for each lag, covariates,
hidden states) are described in ([Berchtold and Raftery
2002](#ref-berchtold2002mixture)) and are outside the scope of `fitMTD`.

### Prediction and simulation

`higherOrderPredict(fit, history)` returns the distribution of the next
state given the most recent states (oldest first), and
`higherOrderSimulate(n, fit, t0)` draws a sequence from the fitted
model, starting from the states in `t0`. Both accept the output of
`fitMTD` and of `fitHigherOrder`, and use the same transition
probabilities as `higherOrderLogLik`, so that summing the logarithms of
the predicted probabilities of the observed states gives back the
log-likelihood.

``` r
fit <- fits[["Koeberg 2"]]
# next direction after two hours from directions 1 and then 2, and after 2 and then 2
round(higherOrderPredict(fit, rbind("1 then 2" = c(1, 2), "2 then 2" = c(2, 2))), 3)
#>              1     2     3     4
#> 1 then 2 0.230 0.699 0.049 0.023
#> 2 then 2 0.037 0.901 0.062 0.000
set.seed(123)
higherOrderSimulate(24, fit, t0 = c(2, 2))
#>  [1] "2" "2" "2" "2" "3" "3" "3" "2" "2" "2" "3" "3" "3" "3" "3" "2" "2" "2" "2"
#> [20] "3" "2" "2" "2" "1"
```

## Higher Order Multivariate Markov Chains

### Introduction

HOMMC model is used for modeling behaviour of multiple categorical
sequences generated by similar sources. The main reference is ([Ching et
al. 2008](#ref-ching2008higher)). Assume that there are s categorical
sequences and each has possible states in M. In nth order MMC the state
probability distribution of the jth sequence at time $t = r + 1$ depend
on the state probability distribution of all the sequences (including
itself) at times $t = r,r - 1,...,r - n + 1$.

\[ x\_{r+1}^{(j)} =
*{k=1}^({s}\_{h=1}){n}*{jk}^({(h)}P\_{h}){(jk)}x\_{r-h+1}^{(k)}, j = 1,
2, …, s, r = n-1, n, … \]

with initial distribution
$x_{0}^{(k)},x_{1}^{(k)},...,x_{n - 1}^{(k)}(k = 1,2,...,s)$. Here

\[ *{jk}^{(h)} , 1j, ks, 1hn and* {k=1}^({s}\_{h=1}){n} \_{jk}^{(h)} =
1, j = 1, 2, 3, … , s. \]

Now we will see the simpler representation of the model which will help
us understand the result of `fitHighOrderMultivarMC` method.

Let
$X_{r}^{(j)} = ((x_{r}^{(j)})^{T},(x_{r - 1}^{(j)})^{T},...,(x_{r - n + 1}^{(j)})^{T})^{T}for\mspace{9mu} j = 1,2,3,...,s.$
Then

\[ $$\begin{pmatrix}
X_{r + 1}^{(1)} \\
X_{r + 1}^{(2)} \\
. \\
. \\
. \\
X_{r + 1}^{(s)}
\end{pmatrix}$$ = $$\begin{pmatrix}
B^{11} & B^{12} & . & . & B^{1s} & \\
B^{21} & B^{22} & . & . & B^{2s} & \\
. & . & . & . & . & \\
. & . & . & . & . & \\
. & . & . & . & . & \\
B^{s1} & B^{s2} & . & . & B^{ss} & \\
 & & & & & 
\end{pmatrix}\begin{pmatrix}
X_{r}^{(1)} \\
X_{r}^{(2)} \\
. \\
. \\
. \\
X_{r}^{(s)}
\end{pmatrix}$$

\]

\[B^{ii} = $$\begin{pmatrix}
{\lambda_{ii}^{(1)}P_{1}^{(ii)}} & {\lambda_{ii}^{(2)}P_{2}^{(ii)}} & . & . & {\lambda_{ii}^{(n)}P_{n}^{(ii)}} & \\
I & 0 & . & . & 0 & \\
0 & I & . & . & 0 & \\
. & . & . & . & . & \\
. & . & . & . & . & \\
0 & . & . & I & 0 & 
\end{pmatrix}$$

\_{mn\*mn} \]

\[ B^{ij} = $$\begin{pmatrix}
{\lambda_{ij}^{(1)}P_{1}^{(ij)}} & {\lambda_{ij}^{(2)}P_{2}^{(ij)}} & . & . & {\lambda_{ij}^{(n)}P_{n}^{(ij)}} & \\
0 & 0 & . & . & 0 & \\
0 & 0 & . & . & 0 & \\
. & . & . & . & . & \\
. & . & . & . & . & \\
0 & . & . & 0 & 0 & 
\end{pmatrix}$$

\_{mn\*mn} ij. \]

### Representation of parameters in the code

$P_{h}^{(ij)}$ is represented as $Ph(i,j)$ and $\lambda_{ij}^{(h)}$ as
Lambdah(i,j). For example: $P_{2}^{(13)}$ as $P2(1,3)$ and
$\lambda_{45}^{(3)}$ as Lambda3(4,5).

### Definition of HOMMC class

``` r
showClass("hommc")
#> Class "hommc" [package "markovchain"]
#> 
#> Slots:
#>                                                                   
#> Name:      order    states         P    Lambda     byrow      name
#> Class:   numeric character     array   numeric   logical character
```

Any element of `hommc` class is comprised by following slots:

1.  states: a character vector, listing the states for which transition
    probabilities are defined.
2.  byrow: a logical element, indicating whether transition
    probabilities are shown by row or by column.
3.  order: order of Multivariate Markov chain.
4.  P: an array of all transition matrices.
5.  Lambda: a vector to store the weightage of each transition matrix.
6.  name: optional character element to name the HOMMC

### How to create an object of class HOMMC

``` r
states <- c('a', 'b')
P <- array(dim = c(2, 2, 4), dimnames = list(states, states))
P[ , , 1] <- matrix(c(1/3, 2/3, 1, 0), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 2] <- matrix(c(0, 1, 1, 0), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 3] <- matrix(c(2/3, 1/3, 0, 1), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 4] <- matrix(c(1/2, 1/2, 1/2, 1/2), byrow = FALSE, nrow = 2, ncol = 2)

Lambda <- c(.8, .2, .3, .7)

hob <- new("hommc", order = 1, Lambda = Lambda, P = P, states = states, 
           byrow = FALSE, name = "FOMMC")
hob
#> Order of multivariate markov chain = 1 
#> states = a b 
#> 
#> List of Lambda's and the corresponding transition matrix (by cols) :
#> Lambda1(1,1) : 0.8
#> P1(1,1) : 
#>           a b
#> a 0.3333333 1
#> b 0.6666667 0
#> 
#> Lambda1(1,2) : 0.2
#> P1(1,2) : 
#>   a b
#> a 0 1
#> b 1 0
#> 
#> Lambda1(2,1) : 0.3
#> P1(2,1) : 
#>           a b
#> a 0.6666667 0
#> b 0.3333333 1
#> 
#> Lambda1(2,2) : 0.7
#> P1(2,2) : 
#>     a   b
#> a 0.5 0.5
#> b 0.5 0.5
```

### Fit HOMMC

`fitHighOrderMultivarMC` method is available to fit HOMMC. Below are the
3 parameters of this method.

1.  seqMat: a character matrix or a data frame, each column represents a
    categorical sequence.
2.  order: order of Multivariate Markov chain. Default is 2.
3.  Norm: Norm to be used. Default is 2.

## A Marketing Example

We tried to replicate the example found in ([Ching et al.
2008](#ref-ching2008higher)) for an application of HOMMC. A soft-drink
company in Hong Kong is facing an in-house problem of production
planning and inventory control. A pressing issue is the storage space of
its central warehouse, which often finds itself in the state of overflow
or near capacity. The company is thus in urgent needs to study the
interplay between the storage space requirement and the overall growing
sales demand. The product can be classified into six possible states (1,
2, 3, 4, 5, 6) according to their sales volumes. All products are
labeled as 1 = no sales volume, 2 = very slow-moving (very low sales
volume), 3 = slow-moving, 4 = standard, 5 = fast-moving or 6 = very
fast-moving (very high sales volume). Such labels are useful from both
marketing and production planning points of view. The data is cointaind
in `sales` object.

``` r
data(sales)
head(sales)
#>      A   B   C   D   E  
#> [1,] "6" "1" "6" "6" "6"
#> [2,] "6" "6" "6" "2" "2"
#> [3,] "6" "6" "6" "2" "2"
#> [4,] "6" "1" "6" "2" "2"
#> [5,] "2" "6" "6" "2" "2"
#> [6,] "6" "1" "6" "3" "3"
```

The company would also like to predict sales demand for an important
customer in order to minimize its inventory build-up. More importantly,
the company can understand the sales pattern of this customer and then
develop a marketing strategy to deal with this customer. Customer’s
sales demand sequences of five important products of the company for a
year. We expect sales demand sequences generated by the same customer to
be correlated to each other. Therefore by exploring these relationships,
one can obtain a better higher-order multivariate Markov model for such
demand sequences, hence obtain better prediction rules.

In ([Ching et al. 2008](#ref-ching2008higher)) application, they choose
the order arbitrarily to be eight, i.e., n = 8. We first estimate all
the transition probability matrices $P_{h}^{ij}$ and we also have the
estimates of the stationary probability distributions of the five
products:.

${\widehat{\mathbf{x}}}^{(1)} = \begin{pmatrix}
0.0818 & 0.4052 & 0.0483 & 0.0335 & 0.0037 & 0.4275
\end{pmatrix}^{\mathbf{T}}$

${\widehat{\mathbf{x}}}^{(2)} = \begin{pmatrix}
0.3680 & 0.1970 & 0.0335 & 0.0000 & 0.0037 & 0.3978
\end{pmatrix}^{\mathbf{T}}$

${\widehat{\mathbf{x}}}^{(3)} = \begin{pmatrix}
0.1450 & 0.2045 & 0.0186 & 0.0000 & 0.0037 & 0.6283
\end{pmatrix}^{\mathbf{T}}$

${\widehat{\mathbf{x}}}^{(4)} = \begin{pmatrix}
0.0000 & 0.3569 & 0.1338 & 0.1896 & 0.0632 & 0.2565
\end{pmatrix}^{\mathbf{T}}$

${\widehat{\mathbf{x}}}^{(5)} = \begin{pmatrix}
0.0000 & 0.3569 & 0.1227 & 0.2268 & 0.0520 & 0.2416
\end{pmatrix}^{\mathbf{T}}$

By solving the corresponding linear programming problems, we obtain the
following higher-order multivariate Markov chain model:

$\mathbf{x}_{r + 1}^{(1)} = \mathbf{P}_{1}^{(12)}\mathbf{x}_{r}^{(2)}$

$\mathbf{x}_{r + 1}^{(2)} = 0.6364\mathbf{P}_{1}^{(22)}\mathbf{x}_{r}^{(2)} + 0.3636\mathbf{P}_{3}^{(22)}\mathbf{x}_{r}^{(2)}$

$\mathbf{x}_{r + 1}^{(3)} = \mathbf{P}_{1}^{(35)}\mathbf{x}_{r}^{(5)}$

$\mathbf{x}_{r + 1}^{(4)} = 0.2994\mathbf{P}_{8}^{(42)}\mathbf{x}_{r}^{(2)} + 0.4324\mathbf{P}_{1}^{(45)}\mathbf{x}_{r}^{(5)} + 0.2681\mathbf{P}_{2}^{(45)}\mathbf{x}_{r}^{(5)}$

$\mathbf{x}_{r + 1}^{(5)} = 0.2718\mathbf{P}_{8}^{(52)}\mathbf{x}_{r}^{(2)} + 0.6738\mathbf{P}_{1}^{(54)}\mathbf{x}_{r}^{(4)} + 0.0544\mathbf{P}_{2}^{(55)}\mathbf{x}_{r}^{(5)}$

According to the constructed 8th order multivariate Markov model,
Products A and B are closely related. In particular, the sales demand of
Product A depends strongly on Product B. The main reason is that the
chemical nature of Products A and B is the same, but they have different
packaging for marketing purposes. Moreover, Products B, C, D and E are
closely related. Similarly, products C and E have the same product
flavor, but different packaging. In this model, it is interesting to
note that both Product D and E quite depend on Product B at order of 8,
this relationship is hardly to be obtained in conventional Markov model
owing to huge amount of parameters. The results show that higher-order
multivariate Markov model is quite significant to analyze the
relationship of sales demand.

``` r

# fit 8th order multivariate markov chain
if (requireNamespace("Rsolnp", quietly = TRUE)) {
object <- fitHighOrderMultivarMC(sales, order = 8, Norm = 2)
}
```

We choose to show only results shown in the paper. We see that $\lambda$
values are quite close, but not equal, to those shown in the original
paper.

### Acknowledgments

We are grateful to Professors Adrian E. Raftery and André Berchtold for
kindly sharing the Koeberg wind-direction and epileptic-seizure series
of ([Berchtold and Raftery 2002](#ref-berchtold2002mixture)), which made
it possible to check `higherOrderLogLik` and `fitMTD` against the
published results.

### References

Anderson, Theodore W, and Leo A Goodman. 1957. “Statistical Inference
about Markov Chains.” *The Annals of Mathematical Statistics*, 89–110.

Berchtold, André. 2001. “Estimation in the Mixture Transition
Distribution Model.” *Journal of Time Series Analysis* 22 (4): 379–97.

Berchtold, André, and Adrian E. Raftery. 2002. “The Mixture Transition
Distribution Model for High-Order Markov Chains and Non-Gaussian Time
Series.” *Statistical Science* 17 (3): 328–56.

Ching, Wai-Ki, Ximin Huang, Michael K Ng, and Tak-Kuen Siu. 2013.
“Higher-Order Markov Chains.” In *Markov Chains*. Springer.

Ching, Wai-Ki, Michael K Ng, and Eric S Fung. 2008. “Higher-Order
Multivariate Markov Chains and Their Applications.” *Linear Algebra and
Its Applications* 428 (2): 492–507.

Csiszár, Imre, and Paul C. Shields. 2000. “The Consistency of the BIC
Markov Order Estimator.” *The Annals of Statistics* 28 (6): 1601–19.

Ghalanos, Alexios, and Stefan Theussl. 2014. *Rsolnp: General Non-Linear
Optimization Using Augmented Lagrange Multiplier Method*.

Katz, Richard W. 1981. “On Some Criteria for Estimating the Order of a
Markov Chain.” *Technometrics* 23 (3): 243–49.

Lèbre, Sophie, and Pierre-Yves Bourguignon. 2008. “An EM Algorithm for
Estimation in the Mixture Transition Distribution Model.” *Journal of
Statistical Computation and Simulation* 78 (1): 1–15.
<https://doi.org/10.1080/00949650701266666>.

MacDonald, Iain L., and Walter Zucchini. 1997. *Hidden Markov and Other
Models for Discrete-Valued Time Series*. Chapman & Hall.

P. J. Avery, and D. A. Henderson. 1999. “Fitting Markov Chain Models to
Discrete State Series.” *Applied Statistics* 48 (1): 53–61.

Raftery, Adrian E. 1985. “A Model for High-Order Markov Chains.”
*Journal of the Royal Statistical Society, Series B* 47 (3): 528–39.

Tong, Howell. 1975. “Determination of the Order of a Markov Chain by
Akaike’s Information Criterion.” *Journal of Applied Probability* 12
(3): 488–97.
