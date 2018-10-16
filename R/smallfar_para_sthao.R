#' @import np
#' @importFrom grDevices rgb
#' @importFrom graphics abline box grid hist lines plot points polygon
#' @importFrom stats confint dexp pexp pnorm qexp quantile sd

estim_GZ.empirical <- function(x, z){
  GZhat <- ecdf(x)(z)
  # GZhat <- GZhat[GZhat > 0]
  # GZhat <- GZhat[GZhat < 1]
  return(GZhat)
}

estim_GZ.kernel <- function(x, z){
  GZhat <- np::npudist(tdat = x, edat = z)$dist  # x = training data & z= evalutaion points of the function npudist (Kernel CDF) in np package
  # GZhat <- GZhat[GZhat > 0]
  # GZhat <- GZhat[GZhat < 1]
  return(GZhat)
}

estim_W <- function(GZhat){
  -log(GZhat)
}

estim_theta <- function(W){
  mean(W)
}

estim_farr_r <- function(theta, rp){
  (1 - theta) * (1 - 1/rp)
}

#
# Asymtpotic confidence intervals for theta hat
#
estim_sigma_theta.wexp <- function(theta){
  #
  # asymptotic stdev of the emprical estimator of theta when W follows an exponential distribution
  # original
  # sigma2_theta <- (1 + theta)^2 / (1 + 2 * theta)+(2 * ( 1 + theta) / ( 2 + theta)) - 2
  # soulivanh
  sigma2_theta <- (1 + theta)^4 / (1 + 2 * theta) + 2 * (1 + theta)^3 / (2 + theta) - 2 * (1 + theta)^2
  sigma_theta <- sqrt(sigma2_theta)
  return(sigma_theta)
}

estim_sigma_farr.wexp <- function(sigma_theta, rp){
  sigma_farr <- (1 - 1/rp) * sigma_theta
  return(sigma_farr)
}

jitter_dat <- function(x, jitter){
  if(length(x) > length(unique(x))){
    warning(paste("making",
                  deparse(substitute(x)),
                  "countinuous to avoid ties"))
    # x <-  x + sample(jitter, length(x)) # adding jittering if data are discrete
    x <- jitter(x)
  }
  return(x)
}

#' Parametric estimation of theta under H0: W ~ exp(theta)
#'
#' \code{estim_theta.wexp} returns an object of class  \code{("thetafit_wexp", "thetafit")} which contains
#' the results of the estimation of theta and of the fit of W ~ exp(theta).
#'
#' This function returns an estimate of theta as defined in Naveau et al (2018) that is used
#' to estimate the far, the fraction of attributable risk for records (with \code{\link{estim_farr.wexp}}).
#' This estimation is made assuming that W = - log(G(Z)) follows an
#' exponentional distribution: W ~ exp(theta). G denotes the Cumulative Distribution Function of
#' the counterfactual variable X.
#'
#' For the full reference, see : \emph{Naveau, P., Ribes, A., Zwiers, F., Hannart, A., Tuel, A., & Yiou, P.
#' Revising return periods for record events in a climate event attribution context.
#' J. Clim., 2018.}, \url{https://doi.org/10.1175/JCLI-D-16-0752.1}
#'
#' @param x the variable of interest in the counterfactual world.
#' @param z the variable of interest in the factual world.
#'
#' @return An object of class \code{("thetafit_wexp", "thetafit")}. It is a list containing the following
#' elements:
#' \describe{
#'   \item{theta_hat}{the estimate of theta}
#'   \item{sigma_theta_hat}{the standard deviation of the estimator of theta
#'   assuming asymptotic gaussianity and obtained via the delta-method}
#'   \item{W}{the estimate of W}
#'   \item{co_test}{the result of the Cox and Oakes test
#'   to test whether W follows an exponential distribution}
#' }
#'
#' @examples
#'  library(evd)
#'
#'  muF <-  1; xiF <- .15; sigmaF <-  1.412538 #  cst^(-xiF) # .05^(-xi1);
#'  # asymptotic limit for the far in this case with a Frechet distributiom
#'  boundFrechet <- frechet_lim(sigma = sigmaF, xi = xiF)
#'  # sample size
#'  size <- 100
#'  # level=.9
#'  set.seed(4)
#'  z = rgev(size, loc = (sigmaF), scale = xiF * sigmaF, shape = xiF)
#'  x = rgev(length(z), loc=(1), scale = xiF, shape=xiF)
#'
#'  rp = seq(from = 2, to = 30, length = 200)
#'  # parametric estimation of far with exponential distribution for theta
#'  theta_fit <- estim_theta.wexp(x = x, z = z)
#'  # Check fit for W ~ exp(theta)
#'  hist(theta_fit)
#'  ecdf(theta_fit)
#'  qqplot(theta_fit)
#'
#'  # Estimate the far
#'  farr_fit.exp <- estim_farr.wexp(theta_hat = theta_fit$theta_hat,
#'   sigma_theta_hat = theta_fit$sigma_theta_hat,
#'    rp = rp)
#'  print(farr_fit.exp)
#'  ylim <- range(boundFrechet, farr_fit.exp$farr_hat)
#'  plot(farr_fit.exp, ylim = ylim, main = "far exponential")
#'  # Theoretical for in this case (Z = sigmaF * X  with X ~ Frechet)
#' lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
estim_theta.wexp <- function(x, z){
  #
  #  x (counterfactual), z (factual), rp = return periods

  # source("CoxOakes.R")
  # rounding numbers to be added to discrete values of x and z
  rounding=seq(from = -0.499, to= 0.499, length = 2 *  max(length(z), length(x))) # rounding numbers to be added to discrete values of x and z
  z <- jitter_dat(z, rounding)
  x <- jitter_dat(x, rounding)
  n <- length(z)


  #
  # Creating W = - log(Ghat(Z)) with Ghat(Z)= average(K(Z-X_i)) for the kernel based on the arctan
  #
  GZhat <- estim_GZ.kernel(x, z)
  W <- estim_W(GZhat)
  # hist(W)
  #
  # Cox and Oakes Test for the exponential distribution
  #
  # KS =  ks.test(W,"pexp",1/theta)
  co_test <- CoxOakes(W[GZhat > 0 & GZhat < 1])
  # print("CoxOakes Test:")
  # print(co_test)

  theta_hat <- 1 / mean(ecdf(x)(z)) - 1
  sigma_theta_hat <- estim_sigma_theta.wexp(theta_hat) / sqrt(n)
  # replace by get_aci()
  # theta_hat_05 <- theta_hat - qnorm(.95) * sigma_theta_hat
  # theta_hat_95 <- theta_hat + qnorm(.95) * sigma_theta_hat

  # cat("theta hat as 1/meanG(Z)-1 with G = empirical cdf  \t\t",
  #     round(theta_hat, 3),
  #     "\t\t 90% CI=", round(theta_hat_05, 3), round(theta_hat_95, 3),
  #     " \n" )
  #
  # cat("Cox and Oakes Test for the exponential distribution (the larger the better) = ",
  #     round(co_test$p.value,4),
  #     " \n")
  theta_fit <- list( theta_hat = theta_hat, sigma_theta_hat = sigma_theta_hat, W = W, co_test = co_test)
  class(theta_fit) <- c("thetafit_wexp", "thetafit")
  return(theta_fit)
}

#' @export
coef.thetafit_wexp <- function(object, ...){
  coef <- object$theta_hat
  names(coef) <- "theta"
  return(coef)
}

#' @export
vcov.thetafit_wexp <- function(object, ...){
  stopifnot( length(object$sigma_theta_hat) == 1)
  vcov <- matrix(object$sigma_theta_hat, ncol = 1, nrow = 1)
  rownames(vcov) <- colnames(vcov) <- "theta"
  return(vcov)
}

#' Parametric estimation of the far under W ~ exp(theta)
#'
#' \code{estim_farr.wexp} returns an object of class  \code{("farrfit_wexp", "farrfit")} which contains
#' the results of the estimation of the far assuming W ~ exp(theta).
#'
#' This function returns an estimate of the far, the fraction of attributable risk for records,
#' as defined in Naveau et al (2018). This estimation is made assuming that W = - log(G(Z)) follows an
#' exponentional distribution: W ~ exp(theta). G denotes the Cumulative Distribution Function of
#' the counterfactual variable X.
#'
#' For the full reference, see : \emph{Naveau, P., Ribes, A., Zwiers, F., Hannart, A., Tuel, A., & Yiou, P.
#' Revising return periods for record events in a climate event attribution context.
#' J. Clim., 2018.}, \url{https://doi.org/10.1175/JCLI-D-16-0752.1}
#'
#' @param theta_hat the value of the estimated theta.
#' The estimation of theta is made assuming W ~ exp(theta).
#' This vector is provided by the function \code{\link{estim_theta.wexp}}
#' @param sigma_theta_hat the value the estimated standard devitation of \code{theta_hat}
#' this vector is provided by the function \code{\link{estim_theta.wexp}}
#' @param rp the return periods for which the far is to be estimated.
#'
#' @return An object of class \code{("farrfit_wexp", "farrfit")}. It is a list containing the following
#' elements:
#' \describe{
#'   \item{rp}{the return periods for which the far is to be estimated.}
#'   \item{farr_hat}{the estimated far for each return period \code{rp}}
#'   \item{sigma__hat}{the standard deviation of the estimator of the far
#'   assuming asymptotic gaussianity and obtained via the delta-method from \code{sigma_theta_hat}}
#' }
#'
#' @examples
#'  library(evd)
#'
#'  muF <-  1; xiF <- .15; sigmaF <-  1.412538 #  cst^(-xiF) # .05^(-xi1);
#'  # asymptotic limit for the far in this case with a Frechet distributiom
#'  boundFrechet <- frechet_lim(sigma = sigmaF, xi = xiF)
#'  # sample size
#'  size <- 100
#'  # level=.9
#'  set.seed(4)
#'  z = rgev(size, loc = (sigmaF), scale = xiF * sigmaF, shape = xiF)
#'  x = rgev(length(z), loc=(1), scale = xiF, shape=xiF)
#'
#'  rp = seq(from = 2, to = 30, length = 200)
#'  # parametric estimation of far with exponential distribution for theta
#'  theta_fit <- estim_theta.wexp(x = x, z = z)
#'  # Check fit for W ~ exp(theta)
#'  hist(theta_fit)
#'  ecdf(theta_fit)
#'  qqplot(theta_fit)
#'
#'  # Estimate the far
#'  farr_fit.exp <- estim_farr.wexp(theta_hat = theta_fit$theta_hat,
#'   sigma_theta_hat = theta_fit$sigma_theta_hat,
#'   rp = rp)
#'  print(farr_fit.exp)
#'  ylim <- range(boundFrechet, farr_fit.exp$farr_hat)
#'  plot(farr_fit.exp, ylim = ylim, main = "far exponential")
#'  # Theoretical for in this case (Z = sigmaF * X  with X ~ Frechet)
#' lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
estim_farr.wexp <- function(theta_hat, sigma_theta_hat, rp){
  farr_hat <- estim_farr_r(rp = rp, theta = theta_hat)
  names(farr_hat) <- paste("farr", rp, sep = "_")
  sigma_farr_hat <- estim_sigma_farr.wexp(sigma_theta =  sigma_theta_hat, rp = rp)
  farr_fit <- list(rp = rp, farr_hat = farr_hat, sigma_farr_hat = sigma_farr_hat)
  class(farr_fit) <- c("farrfit_wexp", "farrfit")
  return(farr_fit)
}

#' @export
coef.farrfit_wexp <- function(object, ...){
  coef <- object$farr_hat
  names(coef) <- names(object$farr_hat)-
  return(coef)
}

#' @export
vcov.farrfit_wexp <- function(object, ...){
  if(length(object$sigma_farr_hat) == 1){
    vcov <- matrix(object$sigma_farr_hat, nrow = 1, ncol = 1)
  } else{
    vcov <- diag(object$sigma_farr_hat)
  }
  rownames(vcov) <- colnames(vcov) <- names(object$farr_hat)
  return(vcov)
}

#' @rdname estim_farr.wexp
#' @export
print.farrfit <- function(x, ...){
  farrdf <- data.frame(rp = get_rp(x),
                      farr_hat = get_farr(x),
                      sigma_farr_hat = get_sigma_farr(x))
  rownames(farrdf) <- NULL
  print("return periods, estimated far and estimated standard deviation for each return period:")
  print(farrdf)
}

get_farr <- function(object, ...){
  UseMethod("get_farr")
}

get_rp <- function(object, ...){
  UseMethod("get_rp")
}

get_sigma_farr <- function(object, ...){
  UseMethod("get_sigma_farr")
}

get_farr.farrfit <- function(object, ...){
  farr <- object$farr_hat
}

get_farr.boot_farrfit <- function(object, ...){
  farr <- apply(object, 1, mean)
}

get_rp.farrfit <- function(object, ...){
  rp <- object$rp
}

get_rp.boot_farrfit <- function(object, ...){
  rp <- dimnames(object)$rp
}

get_sigma_farr.farrfit <- function(object, ...){
  farr <- object$sigma_farr_hat
}

get_sigma_farr.boot_farrfit <- function(object, ...){
  farr <- apply(object, 1, sd)
}

#' @param x an  object of class (\code{"farfit"})
#' @param ... additional arguments for the plot.
#' @rdname estim_farr.wexp
#' @export
plot.farrfit <-function(x, ...){

  farr <- get_farr(x)
  rp <- get_rp(x)
  #
  # quantile for two side 90 confidence interval
  #
  ci90 <- confint(x, names(farr), level = 0.90)
  ellipsis <- list(...)
  if(is.null(ellipsis$ylim)){
    ymin <- min(farr,  ci90, na.rm = TRUE)
    ymax <- max(farr,  ci90, na.rm = TRUE)
    ylim <- c(ymin, ymax)
  } else {
    ylim = ellipsis$ylim
    ellipsis$ylim <- NULL
  }
  if(is.null(ellipsis$col)){
    col = rgb(red = 0.5, green = .1,blue = 0.3,alpha = .5)
  } else {
    col = ellipsis$col
    ellipsis$col <- NULL
  }
  lr <- length(rp)
  args_plot <- list(x =  rp,
                    y = farr,
                    xlab = "Return period r",
                    ylab = "far(r)",
                    type = "n",
                    ylim = ylim)
  args_plot <- c(args_plot, ellipsis)
  do.call(plot, args_plot)
  grid(lwd = 3)
  args_polygon <- list(x = c(rp, rp[lr:1]),
                       y = c(ci90[, 1], ci90[lr:1, 2]),
                       col = col, border=F)
  args_polygon <- c(args_polygon, ellipsis)
  do.call(polygon, args_polygon)
  args_lines <- list(x =  rp,
                    y = farr,
                    col = "black")
  args_lines <- c(args_lines, ellipsis)
  do.call(lines, args_lines)
  #title("far(r)")
}

#' @rdname estim_theta.wexp
#' @export
hist.thetafit_wexp <- function(x, ...){
  #
  # Checking (histogram) if -log(G(Z)) follows an exponentialdistribution
  #
  W <- x$W[is.finite(x$W)]
  hist(W,
       freq = F,
       xlab = "W = -log(G(Z))",
       xlim = c(0, max(W_normalized)),
       main = "Exponential fit for W",
       breaks = max(9, length(W)/10),
       col="lightgray"
  )
  grid(lwd = 3)
  box()
  xx <- seq(from = 0, to = max(W_normalized), length=200)
  lines(xx, dexp(xx, rate = 1/x$theta_hat), col="darkblue", lwd=3, lty=2)
}

#' @export
qqplot.default <- stats::qqplot

#' Generic qqplot function
#' @param x the object to compute the qqplot from.
#' @param ... additional arguments if necessary.
#' @export
#' @export
qqplot <-  function(x, ...){
  UseMethod("qqplot")
}

#' @param ... additional arguments for the plot.
#' @rdname estim_theta.wexp
#' @export
qqplot.thetafit_wexp <- function(x, ...){
  #
  # Checking (qqplot) if -log(G(Z)) follows an exponentialdistribution
  #
  W <- x$W[is.finite(x$W)]
  ll <- length(W)
  pp <- ((1:ll) - 0.5) / ll
  expected <- qexp(pp, rate= 1 / x$theta)
  ci90 <- confint(x)
  expected05 <- qexp(pp, rate = 1/ max(ci90[, 1], 1E-6))
  expected95 <- qexp(pp, rate = 1/ max(ci90[, 2], 1E-6))
  observed = sort(W)
  xylim <- range(expected, observed)
  plot(observed, expected,
       xlim = xylim, ylim = xylim,
       xlab = "Observed", ylab = "Expected",
       main = "Exponential QQ plot for  W=-log(G(Z))",
       col="darkblue", pch = 20, cex = 1)
  grid(lwd = 3)
  box()
  polygon(c(observed, observed[ll:1]), c(expected95, expected05[ll:1]), col = rgb(red=0,green=.0,blue=.5,alpha=.3), border=F)
  abline(a=0, b=1 , col="red", lwd=3)
  points(observed, expected, col="darkblue", pch = 20, cex = 1)
}

#' @export
ecdf.default <- function(x, ...){
  stats::ecdf(x)
}

#' Generic ecdf function
#' @param x the object to compute the ecdf from.
#' @param ... additional arguments if necessary.
#' @export
ecdf <-  function(x, ...){
  UseMethod("ecdf")
}

#' @rdname estim_theta.wexp
#' @export
ecdf.thetafit_wexp <- function(x, ...){
  W <- x$W[is.finite(x$W)]
  plot(ecdf(W),
  xlab = "W = -log(G(Z))",
  ylab = "",
  main = "Theoretical and Empirical CDFs", ...)
  grid(lwd = 3);box()
  xx <- seq(from=0 , to=max(W), length=200)
  lines(xx, pexp(xx, rate = 1 / x$theta_hat), col = "darkgray", lwd = 4, lty = 2)
}

#' Parametric estimation of theta under W ~ exp(theta) in bootstrap samples
#'
#' \code{estim_theta.wexp} returns a list which contains
#' the bootstrap estimates of theta assuming W ~ exp(theta).
#'
#' This function returns bootstrap estimates of theta as defined in Naveau et al (2018) that is used
#' to estimate the far, the fraction of attributable risk for records (with \code{\link{boot_farr_fit.wexp}}).
#' This estimation is made assuming that W = - log(G(Z)) follows an
#' exponentional distribution: W ~ exp(theta). G denotes the Cumulative Distribution Function of
#' the counterfactual variable X. The bootstrap samples of theta are obtained by resampling bootstrap.
#' The first bootrtrap sample corresponds to the original dataset of x and y.
#'
#' For the full reference, see : \emph{Naveau, P., Ribes, A., Zwiers, F., Hannart, A., Tuel, A., & Yiou, P.
#' Revising return periods for record events in a climate event attribution context.
#' J. Clim., 2018.}, \url{https://doi.org/10.1175/JCLI-D-16-0752.1}
#'
#' @param x the variable of interest in the counterfactual world.
#' @param z the variable of interest in the factual world.
#' @param B the number of bootstrap samples to draw.
#'
#' @return An a list containing the following
#' elements:
#' \describe{
#'   \item{theta_boot}{a vector with the estimates of theta for each
#'   bootstrap sample}
#'   \item{xboot}{a list where each element contains the bootstrap sample for x}
#'   \item{zboot}{a list where each element contains the bootstrap sample for z}
#' }
#'
#' @examples
#'  library(evd)
#'
#'  muF <-  1; xiF <- .15; sigmaF <-  1.412538 #  cst^(-xiF) # .05^(-xi1);
#'  # asymptotic limit for the far in this case with a Frechet distributiom
#'  boundFrechet <- frechet_lim(sigma = sigmaF, xi = xiF)
#'  # sample size
#'  size <- 100
#'  # level=.9
#'  set.seed(4)
#'  z = rgev(size, loc = (sigmaF), scale = xiF * sigmaF, shape = xiF)
#'  x = rgev(length(z), loc=(1), scale = xiF, shape=xiF)
#'
#'  rp = seq(from = 2, to = 30, length = 200)
#'  # Resampling bootstrap for the estimation of the far assuming an exponential distribution for theta
#'  theta_boot.exp <- boot_theta_fit.wexp(x = x, z = z, B = 10)
#'
#'  # Estimate the far from the bootstrap samples of theta
#'  boot_farr.exp <- boot_farr_fit.wexp(theta_boot = theta_boot.exp$theta_boot , rp = rp)
#'  confint(boot_farr.exp)
#'  print(boot_farr.exp)
#'  ylim <- range(boundFrechet, boot_farr.exp)
#'  plot(boot_farr.exp, ylim = ylim, main = "boot far exponential")
#'  # Theoretical for in this case (Z = sigmaF * X  with X ~ Frechet)
#' lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
boot_theta_fit.wexp <- function(x, z , B = 100){
  xboot <- lapply(seq.int(B), function(b) sample(x, length(x), replace = TRUE))
  zboot <- lapply(seq.int(B), function(b) sample(z, length(z), replace = TRUE))
  xboot[[1]] <- x
  zboot[[1]] <- z
  theta_boot <- mapply(FUN = function(x, z){
    theta_fit <- estim_theta.wexp(x = x, z = z)
    return(theta_fit$theta_hat)
    }, xboot, zboot, SIMPLIFY = TRUE)
  list(theta_boot = theta_boot, xboot = xboot, zboot = zboot)
}

#' Parametric estimation of the far from bootstrap samples of theta assuming W ~ exp(theta) s
#'
#' \code{boot_farr_fit.wexp} returns an object of class \code{("boot_farrfit.wexp", "boot_farrfit", "farrfit")}
#' which contains the estimates of the far for each bootstrap sample and
#' for different return periods \code{rp}
#'
#' This function returns bootstrap estimates of far, the fraction of attributable risk for records,
#' as defined in Naveau et al (2018). This estimation is made assuming that W = - log(G(Z)) follows an
#' exponentional distribution: W ~ exp(theta). G denotes the Cumulative Distribution Function of
#' the counterfactual variable X. The far is estimate from the bootstrap samples of theta
#' that are obtained by resampling bootstrap.
#' The first bootrtrap sample corresponds to the original dataset of x and y.
#'
#' For the full reference, see : \emph{Naveau, P., Ribes, A., Zwiers, F., Hannart, A., Tuel, A., & Yiou, P.
#' Revising return periods for record events in a climate event attribution context.
#' J. Clim., 2018.}, \url{https://doi.org/10.1175/JCLI-D-16-0752.1}
#'
#' @param theta_boot the results obtained from the function \code{boot_theta_fit.wexp}.
#' @param rp the return periods for which theta is to be estimated.
#'
#' @return An object of class \code{("boot_farrfit.wexp", "boot_farrfit", "farrfit")}.
#' It is a matrix where each columm contains the far estimated for the returns periods \code{rp}
#' for a given bootstrap sample.
#'
#' @examples
#'  library(evd)
#'
#'  muF <-  1; xiF <- .15; sigmaF <-  1.412538 #  cst^(-xiF) # .05^(-xi1);
#'  # asymptotic limit for the far in this case with a Frechet distributiom
#'  boundFrechet <- frechet_lim(sigma = sigmaF, xi = xiF)
#'  # sample size
#'  size <- 100
#'  # level=.9
#'  set.seed(4)
#'  z = rgev(size, loc = (sigmaF), scale = xiF * sigmaF, shape = xiF)
#'  x = rgev(length(z), loc=(1), scale = xiF, shape=xiF)
#'
#'  rp = seq(from = 2, to = 30, length = 200)
#'  # Resampling bootstrap for the estimation of far assuming an exponential distribution for theta
#'  theta_boot.exp <- boot_theta_fit.wexp(x = x, z = z, B = 10)
#'
#'  # Estimate the far from the bootstrap samples of theta
#'  boot_farr.exp <- boot_farr_fit.wexp(theta_boot = theta_boot.exp$theta_boot , rp = rp)
#'  confint(boot_farr.exp)
#'  print(boot_farr.exp)
#'  ylim <- range(boundFrechet, boot_farr.exp)
#'  plot(boot_farr.exp, ylim = ylim, main = "boot far exponential")
#'  # Theoretical for in this case (Z = sigmaF * X  with X ~ Frechet)
#' lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
boot_farr_fit.wexp <- function(theta_boot , rp){
  farr_boot <- vapply(theta_boot,
                     FUN = function(theta) estim_farr_r(rp = rp, theta = theta),
                     FUN.VALUE = rp)
  dimnames(farr_boot) <- list(rp = rp, B = seq.int(length(theta_boot)))
  class(farr_boot) <- c("boot_farrfit.wexp", "boot_farrfit", "farrfit")
  return(farr_boot)
}

#' Compute confidence intervals from bootstrap samples of the far
#'
#' \code{confint.boot_farrfit} returns an matrix that contains the confidence interval
#' for the far computed empirically from the bootstrap estimates of the far
#' for different return periods \code{rp}
#'
#' This function returns a two-column matrix that contains the confidence intervals
#' for the far computed empirically from the bootstrap estimates of the far.
#'
#' @param object an object with the class \code{boot_farrfit}. It is a matrix containing
#' the bootstrap samples of the far. Each line of this matrix represent the estimated far for
#' different return periods \code{rp} and colum corresponds to a differents bootstrap samples of the
#' data.
#' @param parm a vector of return levels for which to compute the confidence interval for the far.
#' The return levels have to be selected from the one present in \code{rownames(object)}.
#'
#' @param level a numerical value between 0 and 1 corresponding to the confidence level
#' of the confidence intervals.
#'
#' @param ... not used.
#'
#' @return A two-column matrix that contains the confidence interval
#' for the far computed empirically from the bootstrap estimates of the far.
#' The first colum is for the lower bound of the confidence interval and the second one
#' for the upper bound. Each line of the matrix represents the condidence interval for a
#' different return period \code{rp}
#'
#' @examples
#'  library(evd)
#'
#'  muF <-  1; xiF <- .15; sigmaF <-  1.412538 #  cst^(-xiF) # .05^(-xi1);
#'  # asymptotic limit for the far in this case with a Frechet distributiom
#'  boundFrechet <- frechet_lim(sigma = sigmaF, xi = xiF)
#'  # sample size
#'  size <- 100
#'  # level=.9
#'  set.seed(4)
#'  z = rgev(size, loc = (sigmaF), scale = xiF * sigmaF, shape = xiF)
#'  x = rgev(length(z), loc=(1), scale = xiF, shape=xiF)
#'
#'  rp = seq(from = 2, to = 30, length = 200)
#'  # Resampling bootstrap for estimation of far assuming an exponential distribution for theta
#'  theta_boot.exp <- boot_theta_fit.wexp(x = x, z = z, B = 10)
#'
#'  # Estimate the far from the bootstrap samples of theta
#'  boot_farr.exp <- boot_farr_fit.wexp(theta_boot = theta_boot.exp$theta_boot , rp = rp)
#'  confint(boot_farr.exp)
#'  print(boot_farr.exp)
#'  ylim <- range(boundFrechet, boot_farr.exp)
#'  plot(boot_farr.exp, ylim = ylim, main = "boot far exponential")
#'  # Theoretical for in this case (Z = sigmaF * X  with X ~ Frechet)
#' lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
confint.boot_farrfit <- function(object, parm = dimnames(object)$rp, level = 0.95, ...){
  irow <- dimnames(object)$rp %in% parm
  alpha <- (1 - level) / 2
  ic <- t(apply(object[irow, ], 1, quantile, probs = c(alpha, 1 - alpha)))
  return(ic)
}

