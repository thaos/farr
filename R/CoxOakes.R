#' Test for exponentiality of Cox and Oakes
#'
#' \code{CoxOakes} returns results of the test of exponentiallity of Cox and Oakes
#'
#' This function performs the test of exponentiallity of Cox and Oakes.
#' For more details, see: Henze and Meintanis (2005, Metrika)
#'
#' @param x a numeric vector of positive values.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{statistic}{the value of the test statistic.}
#'   \item{p.value}{the p-value of the test.}
#' }
#' @examples
#' set.seed(1)
#' x <- rexp(1000)
#' CoxOakes(x)
#' @export
CoxOakes <- function(x){
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1)
    stop("not enough 'x'  data")
  if(sum(x<=0)>0)
    stop(" 'x' has to be positive to be exponentially distributed")
  y <- x / mean(x)
  STATISTIC <- n + sum((1 - y) * log(y))
  NORMAL <- sqrt(6 / n) * (STATISTIC / pi)
  #
  # under the null hyp (exponential)
  # (6/n)^1/2)}*(STATISTIC/pi) is asymtpotically a standard normal distribution
  #
  PVAL <- 2 * (1 - pnorm(abs(NORMAL)))
  #PVAL= (1-pnorm(abs(NORMAL)))
  OUT <- list(statistic = STATISTIC, p.value = PVAL)
  return(OUT)
}
