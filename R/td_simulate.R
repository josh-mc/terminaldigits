#' Monte Carlo simulations for independence of terminal digits
#'
#' The `td_simulate` function performs Monte Carlo simulations to assess
#' type I errors and power for tests of independence of terminal digits for
#' various truncated continuous distributions.
#'
#' Monte Carlo simulations for the null hypothesis are implemented for contingency
#' tables with fixed margins using algorithm ASA 144 (Agresti, Wackerly, and
#' Boyett, 1979; Boyett 1979).
#'
#' @param distribution A string specifying the distribution from which to
#'     draw data for simulations. Options include "normal", "uniform",
#'     and "exponential".
#' @param duplicates A number between 0 and 1 specifying the proportion of
#'     data to be comprised by duplicates. The default value is 0. This is
#'     appropriate for testing type I errors. For testing power, a value
#'     greater than 0 should be entered. For example, entering '0.05' would
#'     ensure that for each simulation, 5% of the data would be comprised by
#'     duplicates.
#' @param n An integer specifying the number of observes to draw from the distribution.
#' @param parameter_1 A numeric value specifying the mean for the normal distribution,
#'     the lower bound of interval for the uniform distribution, or the rate for the
#'     exponential distribution.
#' @param parameter_2 A numeric value specifying the standard deviation for the normal
#'    distribution or the upper bound of the interval for the uniform distribution.
#' @param decimals an integer specifying the number of decimals (including 0)
#'    to which the values drawn from the distribution should be truncated.
#' @param significance a number between 0 and 1 defining the level for
#'     statistical significance. The default is set to 0.05.
#' @param reps an integer specifying the number of Monte Carlo simulations to
#'     implement under the null for each draw. The default is set to 500 but
#'     this is only appropriate for initial exploration.
#' @param simulations an integer specifying the number of Monte Carlo
#'     simulations to perform, i.e. the number of draws from the specified
#'     distribution to be tested. The default is set to 300 but this is only
#'     appropriate for initial exploration.
#' @param tolerance sets an upper bound for rounding errors when evaluating
#'    whether a statistic for a simulation is greater than or equal to the
#'    statistic for the observed data. The default is identical to the tolerance
#'    set for simulations in the `chisq.test` function from the `stats`
#'    package in R.
#'
#' @return A list containing the following components:
#'
#' \item{method}{method employed}
#' \item{distribution}{the distribution}
#' \item{Chisq}{proportion of p-values less than or equal to defined
#'     significance level for Pearson's chi-squared test of independence}
#' \item{G2}{proportion of p-values less than or equal to defined
#'     significance level for log-likelihood ratio test of independence}
#' \item{FT}{proportion of p-values less than or equal to defined
#'     significance level for Freeman-Tukey test of independence}
#' \item{RMS}{proportion of p-values less than or equal to defined
#'     significance level for root-mean-squared test of independence}
#' \item{O}{proportion of p-values less than or equal to defined
#'     significance level for occupancy test of independence}
#' \item{AF}{proportion of p-values less than or equal to defined
#'     significance level for average frequency test of independence}
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
#' @export
#'
#' @examples
#'
#' td_simulate(distribution = "normal",
#' n = 50,
#' parameter_1 = 100,
#' parameter_2 = 1,
#' decimals = 1,
#' reps = 100,
#' simulations = 100)
#'
#' @useDynLib terminaldigits
#' @importFrom Rcpp sourceCpp



td_simulate <- function(distribution,
                        duplicates = 0,
                        n,
                        parameter_1,
                        parameter_2 = NULL,
                        decimals,
                        significance = 0.05,
                        reps = 500,
                        simulations = 300,
                        tolerance = 64 * .Machine$double.eps)  {

  if(!(distribution %in% c("normal", "exponential", "uniform"))) {
    stop("Specify a valid input for the distribution, i.e. 'normal',
         'exponential', or 'uniform'.")}

  if(duplicates < 0 | duplicates > 1) {
    stop("Specify a value for duplicates between 0 and 1")
  }

  if(significance < 0 | significance > 1) {
    stop("Specify a value for duplicates between 0 and 1")
  }

  if(distribution == "normal" & is.null(parameter_2)) {
    stop("Specify parameter_2 (standard deviation)")
  }

  if(distribution == "uniform" & is.null(parameter_2)) {
    stop("Specify parameter_2 (upper bound of interval)")
  }

  if(distribution == "exponential") {
    parameter_2 = 5 #This is a dummy value b/c the C++ function requires a value
    #though the value isn't used.
  }

  dist <- switch(distribution,
                 "normal" = 1,
                 "uniform" = 2,
                 "exponential" = 3)

  out <- perm_basic(distribution = dist,
             duplicates = duplicates * n,
             set_n = n - (duplicates * n),
             set_mean = parameter_1,
             set_sd = parameter_2,
             decimals = decimals,
             reps = reps,
             times = simulations,
             tolerance = tolerance)

  groups <- names(out)

  for(i in groups)  {

    out[[i]][1] <- mean(out[[i]] <= significance)

  }

  val <- list(method = "Monte Carlo simulations for independence of terminal digits",
              distribution = distribution,
              Chisq = out$d_chi_p[1],
              G2 = out$d_g2_p[1],
              FT = out$d_ft_p[1],
              RMS = out$d_rms_p[1],
              O = out$d_perm_p[1],
              AF = out$d_av_fre_p[1])

  return(val)

  }



