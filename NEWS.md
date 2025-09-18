# markovchain 0.10.2

# markovchain 0.10.1

# markovchain 0.10.0

-  Uptick pandoc requirements and handling of sparse matrices

# News for version 0.9.6

-   Handling change of pandoc requirements


# News for version 0.9.5

-   Downtick R requirements


# News for version 0.9.4

-   Corrected strange characters


# News for version 0.9.3

-   Generalized application of requireNamespace(..., quietly = TRUE)
-   Other fixes to comply to newer CRAN requirements
-   Move to MIT license

# News for version 0.9.2

-   Add RcppParallel flags to PKG_LIBS

# News for version 0.9.1

## Current changes

-   Uptick Matrics requirements and modified Changelogs

## old changes

-   2022-09-23 0.9.1 Uptick Matrix reqs

-   2022-07-01 0.9.0 Bugfix a state classification error in Rcpp

-   2022-05-21 0.8.9 Removal of Matlab package dependency

-   2021-05-7 0.8.6 Fix a bug in markovchainListFit that made confusion
    between lists and data.frames

-   2020-12-04 0.8.5-2 Fixing unavailable software issues and language
    glitches

-   2020-09-21 0.8.5-1 Coping with etm unavailability

-   2020-05-21 0.8.5 Fixed DoF in verify markov property and supported
    input in compare teorethical

-   2020-05-04 0.8.4.1 2022-09-23 0.9.1

-   Uptick Matrix reqs 2022-07-01 0.9.0

-   Bugfix a state classification error in Rcpp 2022-05-21

-   0.8.9 Removal of Matlab package dependency

-   2021-05-7 0.8.6 Fix a bug in markovchainListFit that made confusion
    between lists and data.frames

-   2020-12-04 0.8.5-2 Fixing unavailable software issues and language
    glitches

-   2020-09-21 0.8.5-1 Coping with etm unavailability

-   2020-05-21 0.8.5 Fixed DoF in verify markov property and supported
    input in compare teorethical

-   2020-05-04 0.8.4.1 Fixed presentation

-   2020-03-16 0.8.4Limiting output lines in vignettes.

-   2020-03-15 0.8.3 Add small changes in code to cope with upcoming R
    4.0.0 (stringsAsFactor=TRUE in data.frame).

-   2019-12-10 0.8.2 Add small changes in code to cope with upcoming R
    4.0.0 (no more check class(x)=='matrix') as well as packages'
    unavailable.

-   2019-08-13 0.7.0 Improves performance and refactors
    `communicatingClasses`, `recurrentClasses`, `transientStates`,
    `is.irreducible`, `canonicForm`, `summary` and `steadyStates`
    methods, porting them to C++ whenever possible and improving the
    algorithmic complexity of the code. Solves a bug with `steadyStates`
    method. Adds the methods `recurrentStates` and `transientClasses`.
    Makes the aforementioned methods work on by column Markov chains.
    Improves tests, adding checking of mathematical structural
    properties and hundreds of random test cases. Fixes documentation
    for `roxygen` and NAMESPACE file for automatic generation using
    `devtools::document()`

-   Bumps Ignacio Cordón as author (ORCID included) 2019-07-01 0.6.9.15
    Fixed confidence interval calculation: true confidence intervals are
    now 1-(1-confidence_interval)/2 Various code refactoring

-   09-12-2018 0.6.9.14 Added plot from MmgraphR Added
    meanFirstPassageTime (thanks to Toni Giorgino) Add orcid Add more
    warning to Statistical Inference Functions

-   12-08-2018 0.6.9.12 Improved Rcpp performance

-   20-04-2018 0.6.9.9 Fixed typo in vignette MAP method now works also
    with lists (issue #141) Fix valgrid error

-   14-08-2017 0.6.9.8-1 Added is.TimeReversible function added
    gm_to_markovchain example

-   10-07-2017 0.6.9.5 Added empirical bayesian estimate Various
    additions from GSOC 2017 (see the new vignette)

-   31-03-2014 Version 0.6.9 Added sort method Revised numeric tolerance
    when creating markovchains Added suggestion for which row to fix

-   16-03-2017 Version 0.6.8 Deep restructuring of statistical tests

-   Add parameter confint to markovchainFit Fixed bug in
    markovchainFitList

-   Handling of NA

-   02-02-2017 Version 0.6.6.2 Add parameter confint to markovchainFit

-   27-01-2017 Version 0.6.6.1 Fixing bug in markovchainListFit

-   22-01-2017 markovchainFit accepts an uneven list now Added confidence intervals
    when markovchainFit is given a matrix 

-   08-12-2016 Added patch to divergence test 

-   20-08-2016 Fully parallelized bootstrapped markovchain fit 

-   08-08-2016 Version 0.6 Added multivariate higher order markov chains Better handlign of steady state analysis on non - recurrent Markov Chains Fixed an error in the igraph conversion 

-   08-07-2016 Fixed C++ 11 variables types 

-   24-06-2016 Version 0.4.5 Speeding up rmarkovchain using parallel and
    RcppParallel library. 

-   14-06-2016 Version 0.4.4.4 Bug fixed for markovchainFit when method = bootstrap

-   09-06-2016 Version 0.4.4.2 added sanitize=false paramter to markovchainFit

-   31-05-2016 Version 0.4.4.1 Improvement of the internal method checkSequence. name method to set and get the names of markovchain object.

-   10-05-2016 Version 0.4.4 rmarkovchain in RCpp (thanks to Deepak and GSOC 2016) Various small fixes

-   05-03-2016 Version 0.4.3.1 fixed a bug in the states classification added options to save output of random sampler in a matrix

-   10-10-2015 Version 0.4.3 fixed an error in plot function

-   08-07-2015 Version 0.3.1 Period to Rcpp (thanks to TAE) communicatingClasses and recurrentClasses (thanks to TAE) Various optimization (thanks to TAE) Initial support for Continuous Time Markov Chains (thanks to SAI) Added new methods: names, !=

-   15-06-2015 Version 0.3 Added a CrashIntro vignette Most probability function rewritten in Rcpp Added standard errors and confidence intervals for MLE (thanks to Tae) Added confidence intervals for bootstap (thanks to Tae) Added bayesian Maximum A Posteriori estimation (thanks to Sai)

-   12-05-2015 Version 0.2.1 Fixed a compatibility issue with R 3 development 

-   12-04-2015 Version 0.2 This is a milestone for markovchain package, since the package project has been selected within the funded GSOC 2015 projects. Thanks to Tae support now the fitting functions have been rewritten in Rcpp.

-   20-03-2015 Version 0.1.3 Fastened the firstpassage time code thanks to Thoralf suggestion

-   01-03-2015 Version 0.1.2 Add GitHub project url 

-   17-02-2015 Version 0.1.1 Fasten markovchain sequence thanks to Mildenberger Thoralf suggestion

-   04-01-2015 Version 0.1.0 It is now possible to fit a markovchain and a markovchainList object from a matrix or data.frame Updated vignettes Added tests

-   21-06-2014 Version 0.0.9.5 Updated vignettes Added a method to convert a square matrix into a markovchain object. 

-   20-04-2014 Version 0.0.9 Updated vignette Added parallel processing for bootstrap estimation

-   09-02-2014 Version 0.0.8 Minor vignette enhancements Added function to find period of a DTMC

-   12-01-2014 Version 0.0.7 Deeply improved vignettes Added predict and summary methods Added function to perform probabilistic analysis

-   31-12-2013 Version 0.0.5 Improved vignettes Add predict methods Add methods for transitory states

-   04-11-2013 Version 0.0.3 Added various method to easily handle markovchain and markovchainList objects Implemented rmarkovchain and bootstrap fit Improved vignettes
