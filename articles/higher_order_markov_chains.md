# Higher order Markov chains

## Higher Order Markov Chains

Continuous time Markov chains are discussed in the CTMC vignette which
is a part of the package.

An experimental `fitHigherOrder` function has been written in order to
fit a higher order Markov chain (Ching et al.
([2008](#ref-ching2008higher))). `fitHigherOrder` takes two inputs

1.  sequence: a categorical data sequence.
2.  order: order of Markov chain to fit with default value 2.

The output will be a `list` which consists of

1.  lambda: model parameter(s).
2.  Q: a list of transition matrices. $Q_{i}$ is the $ith$ step
    transition matrix stored column-wise.
3.  X: frequency probability vector of the given sequence.

Its quadratic programming problem is solved using `solnp` function of
the Rsolnp package ([Ghalanos and Theussl 2014](#ref-pkg:Rsolnp)).

``` r
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  data(rain)
  rain_small <- rain$rain[1:150]
  fitHigherOrder(rain_small, 2)
}
#> $lambda
#> [1] 0.7733301 0.2266699
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

Two caveats apply. First, `fitHigherOrder` chooses the weights $\lambda$
by least squares on the stationary distribution, not by maximum
likelihood, so the value returned is the log-likelihood *of the fitted
model* and not the maximum attainable one. It can therefore be lower for
a higher order than for a lower one, which cannot happen for
maximum-likelihood fits of nested models. Second, the number of
parameters used for the criteria is $k\, r(r - 1) + (k - 1)$, with $r$
the number of states.

The example compares orders one to three on the Alofi Island daily
rainfall and on the preproglucacon DNA sequence, both analysed by ([P.
J. Avery and D. A. Henderson 1999](#ref-averyHenderson)). Both criteria
select the first-order model for both sequences. The values are computed
by the package; they are not claimed to reproduce those of the original
paper.

``` r
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  compareOrders <- function(sequence, orders = 1:3) {
    fits <- lapply(orders, function(k) fitHigherOrder(sequence, k))
    out <- sapply(fits, function(f)
      unlist(higherOrderLogLik(sequence, f, start = max(orders) + 1)[
        c("logLik", "deviance", "AIC", "BIC", "npar")]))
    colnames(out) <- paste("order", orders)
    round(out, 1)
  }
  data(rain)
  print(compareOrders(rain$rain))
  data(preproglucacon)
  print(compareOrders(preproglucacon$preproglucacon))
}
#>          order 1 order 2 order 3
#> logLik   -1038.1 -1047.8 -1047.8
#> deviance  2076.1  2095.5  2095.5
#> AIC       2088.1  2121.5  2135.5
#> BIC       2118.1  2186.5  2235.5
#> npar         6.0    13.0    20.0
#>          order 1 order 2 order 3
#> logLik   -2026.0 -2024.6 -2051.6
#> deviance  4052.0  4049.1  4103.1
#> AIC       4076.0  4099.1  4179.1
#> BIC       4140.3  4233.1  4382.7
#> npar        12.0    25.0    38.0
```

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
squares, whereas the MTD model uses a single matrix $Q$ for all lags and
is estimated by maximum likelihood. The comparison below therefore
evaluates, with `higherOrderLogLik`, the first-order chain estimated on
the observations entering the likelihood and the MTD(2) model with the
weights and the matrix $Q$ printed in the paper (whose rows are the
departure states, hence the transposition).

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
are part of the unit tests of the package.

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

### References

Berchtold, André, and Adrian E. Raftery. 2002. “The Mixture Transition
Distribution Model for High-Order Markov Chains and Non-Gaussian Time
Series.” *Statistical Science* 17 (3): 328–56.

Ching, Wai-Ki, Ximin Huang, Michael K Ng, and Tak-Kuen Siu. 2013.
“Higher-Order Markov Chains.” In *Markov Chains*. Springer.

Ching, Wai-Ki, Michael K Ng, and Eric S Fung. 2008. “Higher-Order
Multivariate Markov Chains and Their Applications.” *Linear Algebra and
Its Applications* 428 (2): 492–507.

Ghalanos, Alexios, and Stefan Theussl. 2014. *Rsolnp: General Non-Linear
Optimization Using Augmented Lagrange Multiplier Method*.

MacDonald, Iain L., and Walter Zucchini. 1997. *Hidden Markov and Other
Models for Discrete-Valued Time Series*. Chapman & Hall.

P. J. Avery, and D. A. Henderson. 1999. “Fitting Markov Chain Models to
Discrete State Series.” *Applied Statistics* 48 (1): 53–61.

Raftery, Adrian E. 1985. “A Model for High-Order Markov Chains.”
*Journal of the Royal Statistical Society, Series B* 47 (3): 528–39.
