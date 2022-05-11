#' Tests of independence and uniformity for terminal digits in a data frame
#'
#' The function `td_tests()` is a wrapper which applies the functions`td_independence()` and
#' `td_uniformity` to a data frame. When a `group` is specified, tests are conducted separated
#' for each group. P-values and p-values adjusted by the false discovery rate (Benjamini
#' and Hochberg, 1995) are reported.
#'
#'@references
#'
#' Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289–300. doi: 10.1111/j.2517-6161.1995.tb02031.x. https://www.jstor.org/stable/2346101.
#'
#' @param data A data frame
#' @param variable A numeric variable. Tests for terminal digits are performed on this variable.
#' @param group A variable used to group the primary variable such that
#'     p-values are calculated separately for each group. The default is set to NULL in which case
#'     p-values are simply calculated for the whole data set.
#' @param decimals an integer specifying the number of decimals. This can be zero if the terminal digit is
#'     not a decimal.
#' @param reps an integer specifying the number of Monte Carlo simulations. The default
#'     is set to 10,000.
#' @param test a string specifying the test of independence. The default is Pearson's
#'     chi-squared statistic ("Chisq"). Also available is the log-likelihood ratio
#'    statistic ("G2"), the Freeman-Tukey statistic ("FT"), and the Root-mean-square
#'    statistic ("RMS").
#' @param tolerance sets an upper bound for rounding errors when evaluating
#'    whether a statistic for a simulation is greater than or equal to the
#'    statistic for the observed data. The default is identical to the tolerance
#'    set for simulations in the `chisq.test` function from the `stats`
#'    package in R.
#'
#'@return A data frame containing the following components:
#'
#' \item{statistic}{the value of the test statistic}
#' \item{p_value_independence}{the simulated p-value for the test of independence}
#' \item{P_value_uniformity}{the simulated p-value for the test of uniformity (chi-squared GOF)}
#' \item{p_value_independence_fdr}{the simulated p-value for the test of independence adjusted via the
#'     false discovery rate (if the `group` variable is specified)}
#' \item{P_value_uniformity}{the simulated p-value for the test of uniformity (chi-squared GOF)
#'     adjusted via the false discovery rate (if the `group` variable is specified)}
#'
#' @export
#'
#' @examples
#'
#'td_tests(decoy, weight, decimals = 2, group = subject, reps = 1000)
#'
#' @useDynLib terminaldigits
#' @importFrom Rcpp sourceCpp
#'
td_tests <- function(data, variable, decimals, group = NULL, reps = 10000, test = "Chisq",
                     tolerance = 64 * .Machine$double.eps) {

  var <- deparse(substitute(variable))
  grp <- deparse(substitute(group))

  if(is.null(data[[grp]])) {

    ind <- td_independence(data[[var]],
                           decimals,
                           reps,
                           test,
                           tolerance)$p.value

    uni <- td_uniformity(data[[var]],
                         decimals,
                         reps,
                         tolerance)$p.value

    val <- list(p_value_independence = ind,
                p_value_uniformity = uni)

  }

  else {

    groups <- unique(data[[grp]])

    results <- data.frame(groups)

    count <- 1

    for(i in groups)  {

      n_data <- data[data[[grp]] == i, ]

      ind <- td_independence(n_data[[var]],
                             decimals,
                             reps,
                             test,
                             tolerance)$p.value

      uni <- td_uniformity(n_data[[var]],
                           decimals,
                           reps,
                           tolerance)$p.value

      results$p_value_independence[count] <- ind
      results$p_value_uniformity[count] <- uni

      count <- count + 1

    }

    results$p_value_independent_fdr <- stats::p.adjust(results$p_value_independence, method = "BH")
    results$p_value_uniformity_fdr <- stats::p.adjust(results$p_value_uniformity, method = "BH")

    val <- results[, c(1, 2, 4, 3, 5)]



  }

  return(val)

}
