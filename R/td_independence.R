#' Test of independence of terminal digits
#'
#' The `td_independence` function tests the independence of terminal digits from
#' preceding digits by constructing a contingency table of counts where rows constitute
#' unique preceding digits and columns constitute unique terminal digits. A test of
#' independence for a contingency tables is then implemented via Monte Carlo
#' simulation.
#'
#' Monte Carlo simulations are implemented for contingency tables with fixed
#' margins using algorithm ASA 144 (Agresti, Wackerly, and Boyett, 1979; Boyett 1979).
#'
#' @param x a numeric vector
#' @param decimals an integer specifying the number of decimals. This can be zero if the terminal digit is
#'     not a decimal.
#' @param reps a positive integer specifying the number of Monte Carlo simulations. The default
#'     is set to 10,000 which may be appropriate for exploratory analysis. A higher
#'     number of simulation should be selected for more precise results.
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
#' @return A list with class "htest" containing the following components:
#'
#' \item{statistic}{the value of the test statistic}
#' \item{p_value}{the simulated p-value for the test}
#' \item{method}{a character string describing the test}
#' \item{data.name}{a character string give the name of the data}
#'
#'
#' @references
#'
#' Agresti, A., Wackerly, D., & Boyett, J. M. (1979). Exact conditional tests
#'     for cross-classifications: approximation of attained significance levels.
#'     Psychometrika, 44(1), 75-83.
#'
#' Boyett, J. M. (1979). Algorithm AS 144: Random r × c tables with
#'     given row and column totals. Journal of the Royal Statistical Society.
#'     Series C (Applied Statistics), 28(3), 329-332.
#'
#' @export
#'
#' @examples
#'
#' td_independence(decoy$weight, decimals = 2)
#'
#' @useDynLib terminaldigits
#' @importFrom Rcpp sourceCpp
#'
#'
#'
td_independence <- function(x,
                            decimals,
                            reps = 10000,
                            test = "Chisq",
                            tolerance = 64 * .Machine$double.eps) {

  if(!class(x) %in% c("numeric", "integer" )) {stop("The vector `x` must be numeric")}

  if (!(test %in% c("Chisq", "G2", "FT", "RMS"))) {
    stop("Specify a valid input for the test, i.e. 'Chisq', 'G2', 'FT', or 'RMS'.")}

  if(reps <= 0) {
    stop("The 'reps' parameter requires a positive integer")
  }


  DNAME <- deparse(substitute(x))

  x <- x[!is.na(x)]
  x <- x[!is.infinite(x)]

  type <- switch(test,
                 "Chisq" = 1,
                 "G2" = 2,
                 "FT" = 3,
                 "RMS" = 4)

  out <- terminal_independence(x,
                               decimals = decimals,
                               reps = reps,
                               tolerance = tolerance,
                               type = type)

  METHOD <- paste(test, "test for independence of terminal digits")

  names(out$statistic) <- test

  val <- list(p.value = out$p_value,
              statistic = out$statistic,
              method = METHOD,
              data.name = DNAME,
              class = "htest")

  class(val) <- "htest"

  return(val)

}
