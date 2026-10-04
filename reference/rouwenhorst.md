# Discretize an AR(1) process into a Markov chain (Rouwenhorst's method)

Approximates the stationary first-order autoregressive process \$\$y_t =
(1-\rho)\alpha + \rho y\_{t-1} + \varepsilon_t, \qquad \varepsilon_t
\overset{\mathrm{iid}}{\sim} \mathcal N(0,\sigma^2)\$\$ by a
finite-state Markov chain, following Rouwenhorst (1995). Unlike
[`tauchen`](tauchen.md), no grid-width parameter is needed and the
method remains accurate for `rho` close to \\\pm 1\\.

## Usage

``` r
rouwenhorst(alpha, sigma, rho, size)
```

## Arguments

- alpha:

  A single finite number: the unconditional mean of the process.

- sigma:

  A single finite positive number: the standard deviation of the
  innovation \\\varepsilon_t\\.

- rho:

  A single number in \\(-1,1)\\: the autocorrelation (persistence) of
  the process.

- size:

  A single integer of at least `2`: the number of grid points (states)
  of the discretized chain.

## Value

A named list with two elements, `chain` and `states`, in the same form
as returned by [`tauchen`](tauchen.md).

## Details

The grid is \\n=\code{size}\\ evenly spaced points spanning
\\\[\alpha-\psi,\\ \alpha+\psi\]\\ with \\\psi=\sigma_y\sqrt{n-1}\\,
where \\\sigma_y=\sigma/\sqrt{1-\rho^2}\\ is the process's unconditional
standard deviation: this particular width (rather than a fixed multiple
of \\\sigma_y\\ as in [`tauchen`](tauchen.md)) is what the method needs
in order to match the AR(1)'s variance and first-order autocorrelation
exactly at every `size`, including for `rho` near \\\pm1\\. The
transition matrix is built recursively. Let \\\theta=(1+\rho)/2\\ and,
for two states, \$\$\Theta_2 = \begin{pmatrix}\theta & 1-\theta\\
1-\theta & \theta\end{pmatrix}.\$\$ For \\m\\ states (\\2\<m\le n\\),
form the \\m\times m\\ matrix \$\$\Theta_m =
\theta\begin{pmatrix}\Theta\_{m-1} & 0\\ 0 & 0\end{pmatrix} +
(1-\theta)\begin{pmatrix}0 & \Theta\_{m-1}\\ 0 & 0\end{pmatrix} +
(1-\theta)\begin{pmatrix}0 & 0\\ \Theta\_{m-1} & 0\end{pmatrix} +
\theta\begin{pmatrix}0 & 0\\ 0 & \Theta\_{m-1}\end{pmatrix},\$\$ then
divide every interior row (all but the first and last) by \\2\\ to
restore row-stochasticity, since those rows receive contributions from
two of the four corner blocks above.

Rouwenhorst's method reproduces the AR(1)'s unconditional variance and
lag-1 autocorrelation \\\rho\\ exactly, for every `size` (Kopecky and
Suen (2010) find it outperforms [`tauchen`](tauchen.md) and several
other methods across the persistence range typically seen in quarterly
macroeconomic and actuarial time series, e.g. discretized short-rate or
inflation processes for reserving and ALM work).

## References

Rouwenhorst, K. G. (1995). Asset pricing implications of equilibrium
business cycle models. In T. F. Cooley (Ed.), *Frontiers of Business
Cycle Research*, 294-330. Princeton University Press.

Kopecky, K. A. and Suen, R. M. H. (2010). Finite state Markov-chain
approximations to highly persistent processes. *Review of Economic
Dynamics*, 13(3), 701-714.

## See also

[`tauchen`](tauchen.md)

## Examples

``` r
out <- rouwenhorst(alpha = 0, sigma = 1, rho = 0.9, size = 5)
out$states
#>    -4.588    -2.294     0.000     2.294     4.588 
#> -4.588315 -2.294157  0.000000  2.294157  4.588315 
pi <- as.numeric(steadyStates(out$chain))
sum(pi * (out$states - sum(pi * out$states))^2) # matches 1/(1-rho^2) closely
#> [1] 5.263158
1 / (1 - 0.9^2)
#> [1] 5.263158
```
