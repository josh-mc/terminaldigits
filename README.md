
<!-- README.md is generated from README.Rmd. Please edit that file -->

# terminaldigits

<!-- badges: start -->
<!-- badges: end -->

The package `terminaldigits` implements simulated tests of uniformity
and independence for terminal digits. For certain parameters,
`terminaldigits` also implements Monte Carlo simulations for type I
errors and power for the test of independence. Simulations are run in
C++ utilizing Rcpp.

## Installation

You can install the development version of `terminaldigits` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("josh-mc/terminaldigits")
```

## Usage

In many cases, terminal digits can be assumed to be uniformly
distributed and independent of preceding digits. A violation of either
of these assumptions may point to a data quality issue.

The following examples are based on a data set taken from the third
round of a decoy experiment involving hand-washing purportedly carried
out in a number of factories in China. For details, see `decoy` and Yu,
Nelson, and Simonsohn (2018).

The `td_uniformity` function tests the assumption of uniformity using
Pearson’s chi-squared statistic for goodness-of-fit.

``` r
library(terminaldigits)

td_uniformity(decoy$weight, decimals = 2, reps = 1000)
#> 
#>  Pearson's chi-squared GOF test for uniformity of terminal digits
#> 
#> data:  decoy$weight
#> Chi-squared = 539.67, p-value = 0.000999
```

The `td_independence` function tests the assumption of independence. The
default statistic is again Pearson’s chi-squared statistic but the
log-likelihood ratio statistic, the Freeman-Tukey statistic, and the
root-mean-square statistic are also available.

``` r
td_independence(decoy$weight, decimals = 2, reps = 1000)
#> 
#>  Chisq test for independence of terminal digits
#> 
#> data:  decoy$weight
#> Chisq = 6422.4, p-value = 0.000999
```

The `td_test` function is a wrapper for the above two functions. For
more details, including a discussion of the `td_simulate` function, see
the package introduction vignette.

## References

Yu, F., Nelson, L., & Simonsohn, U. (2018, December 5). “In Press at
Psychological Science: A New ‘Nudge’ Supported by Implausible Data.”
DataColoda 74. <http://datacolada.org/74>
