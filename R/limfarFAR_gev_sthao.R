#' Theoretical far and FAR for three types of changes in a GEV distribution
#'
#' A set of functions to return the theoretical far and the theorecial FAR in 3 specific cases when X and Z
#' follow a GEV distribution.
#'
#' These functions (\code{gumble_farr}, \code{frechet_farr}, \code{weibull_farr},
#' \code{gumble_FAR}, \code{frechet_FAR}, \code{weibull_FAR})
#' return the theoretical far and the FAR in 3 specific cases when X and Z
#' follow a GEV distribution. The three specific cases described in Naveau et al. (2018) are :
#' \enumerate{
#'   \item X ~ Gumbel with cdf G(x) = exp[-exp(-x)]
#'    with x real, mu >= 0  and Z = X + mu
#'   \item X ~ Frechet with cdf G(x) = exp[-x^(-1 / xi)]
#'    with x > 0, sigma > 0  and Z = sigma * X
#'   \item X ~ Weibull with cdf G(x) = exp[-max(0, 1 + xi * x)^(-1 / xi)]
#'    with x real, mu >= 0  and Z = mu + sigma * X
#' }
#' For these three cases, the limits of the far and the FAR coincide when the return period \code{r} grows to
#' infinity and when the threshold \code{u} move towards the upper support bound of the distribution.
#' The functions \code{gumble_lim}, \code{frechet_lim} and \code{weibull_lim} compute those limits.
#'
#' For the full reference, see : \emph{Naveau, P., Ribes, A., Zwiers, F., Hannart, A., Tuel, A., & Yiou, P.
#' Revising return periods for record events in a climate event attribution context.
#' J. Clim., 2018.}, \url{https://doi.org/10.1175/JCLI-D-16-0752.1}
#'
#' @param r the return periods for which the far is to be estimated.
#' @param u the thresholds for which the FAR is to be estimated.
#' @param mu the shift in the GEV location parameter between X and Z
#' @param sigma the multiplying factor in the GEV scale parameter between X and Z
#' @param xi the value of th GEV shape parameter
#'
#' @examples
#'  library(evd)
#'  # sample size
#'  # size <- 100
#'
#'  # Gumbel
#'  # Z = muDelta + X
#'  # x = rgev(size, loc = 0, scale = 1, shape = 0)
#'  # z = rgev(size, loc = muDelta, scale = 1, shape = 0)
#'  muG <-  0; sigmaG <-  1; xiG = 0 ; muDelta = 0.5
#'  r = seq(from = 2, to = 100, length = 100)
#'  u  = qgev(1 - (1 / r), loc = muDelta, scale = 1, shape = 0)
#'
#'  # asymptotic limit for the far and the FAR in this case with a Gumble distributiom
#'  glim <- gumble_lim(mu  = muDelta)
#'  # theoretical far in this case with a Frechet distributiom
#'  gfarr <- gumble_farr(r = r, mu = muDelta)
#'  # theoretical FAR in this case with a Frechet distributiom
#'  gFAR <- gumble_FAR(u = u, mu = muDelta)
#'  par(mfrow = c(1, 2))
#'  ylim = range(gfarr, gFAR, glim)
#'  plot(r, gfarr, ylim = ylim, main = "far theoretical value")
#'  abline(h = glim)
#'  plot(u, gFAR, ylim = ylim, main = "FAR theoretical value")
#'  abline(h = glim)
#'
#'  # Frechet
#'  # Z = sigmaDelta * X
#'  # x = rgev(size, loc = 1, scale = xiF, shape = xiF)
#'  # z = rgev(size, loc = sigmaDelta, scale = xiF * sigmaDelta, shape = xiF)
#'  muF <-  1; sigmaF <- xiF <- .15 ; sigmaDelta <-  1.412538
#'  r = seq(from = 2, to = 100, length = 100)
#'  u  = qgev(1 - (1 / r), loc = muF, scale = sigmaF, shape = xiF)
#'
#'  # asymptotic limit for the far and the FAR in this case with a Frechet distributiom
#'  flim <- frechet_lim(sigma = sigmaDelta, xi = xiF)
#'  # theoretical far in this case with a Frechet distributiom
#'  ffarr <- frechet_farr(r = r, sigma = sigmaDelta, xi = xiF)
#'  # theoretical FAR in this case with a Frechet distributiom
#'  fFAR <- frechet_FAR(u = u, sigma = sigmaDelta, xi = xiF)
#'  par(mfrow = c(1, 2))
#'  ylim = range(ffarr, fFAR, flim)
#'  plot(r, ffarr, ylim = ylim, main = "far theoretical value")
#'  abline(h = flim)
#'  plot(u, fFAR, ylim = ylim, main = "FAR theoretical value")
#'  abline(h = flim)
#'
#'  # Weibull
#'  # Z = muDelta + sigmaDelta * X
#'  # x = rgev(size, loc = 0, scale = 1, shape = xiW)
#'  # z = rgev(size, loc = muDelta, scale =  sigmaDelta, shape = xiW)
#'  muW <-  0 ; sigmaW <- 1 ; xiW <- -0.05;
#'  sigmaDelta <-  1.1 ; muDelta <- (1 - sigmaDelta) / (-xiW)
#'  r = seq(from = 2, to = 100, length = 100)
#'  u  = qgev(1 - (1 / r), loc = muW, scale = sigmaW, shape = xiW)
#'
#'  # asymptotic limit for the far and the FAR in this case with a Weibull distributiom
#'  wlim <- weibull_lim(sigma = sigmaDelta, xi = xiW)
#'  # theoretical far in this case with a Frechet distributiom
#'  wfarr <- weibull_farr(r = r, mu = muDelta, sigma = sigmaDelta, xi = xiW)
#'  # theoretical FAR in this case with a Frechet distributiom
#'  wFAR <- weibull_FAR(u = u, mu = muDelta, sigma = sigmaDelta, xi = xiW)
#'  par(mfrow = c(1, 2))
#'  ylim = range(wfarr, wFAR, wlim)
#'  plot(r, wfarr, ylim = ylim, main = "far theoretical value")
#'  abline(h = wlim)
#'  plot(u, wFAR, ylim = ylim, main = "FAR theoretical value")
#'  abline(h = wlim)
#' @export
gumble_farr <- function(r, mu){
  cst <- exp(-mu)
  1 - (1 + cst * (r-1)) / r
}

#' @rdname gumble_farr
#' @export
frechet_farr <- function(r, sigma, xi){
  cst <- sigma^(-1 / xi)
  1 - (1 + cst * (r - 1)) / r
}

#' @rdname gumble_farr
#' @export
weibull_farr <- function(r, mu, sigma, xi){
  cst <- sigma^(-1 / xi)
  1 - (1 + cst * (r - 1)) / r
}

#' @rdname gumble_farr
#' @export
gumble_FAR <- function(u, mu){
  top <- 1 - exp(-exp(-u))
  bottom <- 1 - exp(-exp(-(u - mu)))
  1 - (top/bottom)
}

#' @rdname gumble_farr
#' @export
frechet_FAR <- function(u, sigma, xi){
  top <- 1 - exp(-(u)^(-1 / xi))
  bottom <- 1 - exp(-(u / sigma)^(-1 / xi))
  1 - (top/bottom)
}

#' @rdname gumble_farr
#' @export
weibull_FAR <- function(u, mu, sigma, xi){
  top <- 1 - exp( -(1 + xi * u)^(-1 / xi))
  cc <- sigma^(1 / xi)
  bottom <- 1 - exp(-cc*(1 + xi * u)^(-1 / xi))
  1 - (top/bottom)
}

#' @rdname gumble_farr
#' @export
gumble_lim <- function(mu){
  1 - exp(-mu)
}

#' @rdname gumble_farr
#' @export
frechet_lim <- function(sigma, xi){
  cst <- sigma^(-1 / xi)
  1 - cst
}

#' @rdname gumble_farr
#' @export
weibull_lim <- function(sigma, xi){
  cst <- sigma^(-1 / xi)
  1 - cst
}

#
# computing far(r) for a GEV
#
gevfar=function(r, mu, sigma, xi){
  if(xi==0) {
    cst= exp(-mu)
  }
  if(xi!=0) {
    cst=(sigma)^(-1/xi)
  }
  out=rep(NA,length(r))
  for(i in 1:length(r)){
    out[i]=1-((1+cst*(r[i]-1))/r[i])
  }
  print(out)
  out <- 1-((1+cst*(r-1))/r)
  print(out)
  return(out)
}
#
# computing FAR(u) for a GEV
#
gevFAR=function(u, mu, sigma, xi){
  m=length(u)
  if(xi==0) {
    top=1-exp(-exp(-u))
    bottom=1-exp(-exp(-(u-mu)))
  }
  if(xi>0) {
    top=1-exp(-(u)^(-1/xi))
    bottom=1-exp(-(u/sigma)^(-1/xi))
  }
  if(xi<0) {
    top=1-exp(-(1+xi*u)^(-1/xi))
    cc=(sigma)^(1/xi)
    bottom=1-exp(-cc*(1+xi*u)^(-1/xi))
  }
  out=1-(top/bottom)
  return(out)
}

