#' 3,320 observations from a decoy experiment
#'
#' A data frame containing 3,320 observations (with NA's) from the third round
#' of a decoy experiment involving hand-washing purportedly carried out in a
#' number of factories in China.
#'
#' This series of experiments was published in
#' an article in Psychological Science in 2018. Subsequently, Frank Yu, Leif
#' Nelson, and Uri Simonsohn argued that the data for the experiments could
#' not be \href{http://datacolada.org/74}{trusted},and Simonsohn developed
#' number-bunching in relation to his analysis of the \href{http://datacolada.org/77}{data}.
#' The article was eventually
#' \href{https://journals.sagepub.com/doi/10.1177/0956797618761374}{retracted}. This data
#' frame consists of the data contained in the tab named "Study3-sanitizer usage(grams)".
#'
#' @format A data frame with 3320 rows and 3 variables:
#' \itemize{
#'   \item subject
#'   \item workroom: The room for which the sanitizer weight is recorded.
#'   \item value: The weight in grams for the sanitizer.
#' }
#'
#'@source \url{https://osf.io/wqp7y}

"decoy"

