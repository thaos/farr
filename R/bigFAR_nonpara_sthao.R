estim_FAR.emp <- function(x, z, u){
    p0 <- 1 - ecdf(x)(u)
    p1 <- 1 - ecdf(z)(u) # + 0.5 / length(z)
    FAR <- 1 - p0 / p1
    class(FAR) <- c("FARfit_emp", "FARfit")
    return(FAR)
}

#' Empirical estimation of the FAR from bootstrap samples
#'
#' \code{boot_FAR_fit.emp} returns an object of class \code{("boot_FARfit.emp", "boot_FARfit", "FARfit")}
#' which contains the estimates of the FAR for each bootstrap sample and
#' for different return periods \code{rp}
#'
#' This function returns bootstrap empirical estimates of FAR, the fraction of attributable risk
#' where events are defined in terms of a threshold exceendance.
#' The FAR is estimate from the bootstrap samples of the data x and z
#' that are obtained by resampling bootstrap.
#' The first bootstrap sample is obtained from the original dataset of x and y.
#'
#' @param x the variable of interest in the counterfactual world.
#' @param z the variable of interest in the factual world.
#' @param u the thresholds used to define the events.
#' @param B the number of bootstrap samples to draw.
#'
#' @return An object of class \code{("boot_FARfit.emp", "boot_FARfit", "FARfit")}.
#' It is a matrix where each columm contains the FAR estimated for the returns periods \code{rp}
#' for a given bootstrap sample.
#'
#' @examples
#'  library(evd)
#'
#'  muF <-  1; xiF <- .15; sigmaF <-  1.412538 #  cst^(-xiF) # .05^(-xi1);
#'  # asymptotic limit for the far in this case with a Frechet distributiom
#'  boundFrechet <- 1 - sigmaF^(-1/xiF)
#'  # sample size
#'  size <- 100
#'  # level=.9
#'  set.seed(4)
#'  z = rgev(size, loc = (sigmaF), scale = xiF * sigmaF, shape = xiF)
#'  x = rgev(length(z), loc=(1), scale = xiF, shape=xiF)
#'
#'  rp = seq(from = 2, to = 30, length = 200)
#'  u  = qgev(1 - (1 / rp),loc = muF, scale = xiF, shape = xiF)
#'
#'  # Resampling bootstrap for the empirical estimation of the FAR
#'  boot_FAR.emp <- boot_FAR_fit.emp(x = x, z = z, u = u, B = 10)
#'  print(boot_FAR.emp)
#'  confint(boot_FAR.emp)
#'
#'  ylim <- range(boundFrechet, boot_FAR.emp)
#'  plot(boot_FAR.emp, ylim = ylim, main = "boot FAR empirical")
#'  # Theoretical FAR for in this case (Z = sigmaF * X  with X ~ Frechet)
#'  lines(rp, frechet_FAR(u = u, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
boot_FAR_fit.emp <- function(x, z , u, B = 100){
  xboot <- lapply(seq.int(B), function(b) sample(x, length(x), replace = TRUE))
  zboot <- lapply(seq.int(B), function(b) sample(z, length(z), replace = TRUE))
  xboot[[1]] <- x
  zboot[[1]] <- z
  FAR_boot <- mapply(FUN = function(x, z){
    FAR_fit <- estim_FAR.emp(x = x, z = z, u = u)
    return(FAR_fit)
    }, xboot, zboot, SIMPLIFY = TRUE)
  dimnames(FAR_boot) <- list(u = u, B = seq.int(B))
  class(FAR_boot) <- c("boot_FARfit.emp", "boot_FARfit", "FARfit")
  return(FAR_boot)
}

#' Compute confidence intervals from bootstrap samples of the FAR
#'
#' \code{confint.boot_FARfit} returns an matrix that contains the confidence interval
#' for the far computed empirically from the bootstrap estimates of the far
#' for different exceedance thresholds \code{u}
#'
#' This function returns a two-column matrix that contains the confidence intervals
#' for the FAR computed empirically from the bootstrap estimates of the FAR.
#'
#' @param object an object with the class \code{boot_FARfit}. It is a matrix containing
#' the bootstrap samples of the far. Each line of this matrix represent the estimated far for
#' different thresholds \code{u} and colum corresponds to a differents bootstrap samples of the
#' data.
#' @param parm a vector of thresholds for which to compute the confidence interval for the FAR.
#' The thresholds have to be selected from the one present in \code{rownames(object)}.
#'
#' @param level a numerical value between 0 and 1 corresponding to the confidence level
#' of the confidence intervals.
#'
#' @param ... not used.
#'
#' @return A two-column matrix that contains the confidence interval
#' for the FAR computed empirically from the bootstrap estimates of the FAR.
#' The first colum is for the lower bound of the confidence interval and the second one
#' for the upper bound. Each line of the matrix represents the condidence interval for a
#' different threshold \code{u}
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
#'  u  = qgev(1 - (1 / rp),loc = muF, scale = xiF, shape = xiF)
#'
#'  # Resampling bootstrap for the empirical estimation of the FAR
#'  boot_FAR.emp <- boot_FAR_fit.emp(x = x, z = z, u = u, B = 10)
#'  print(boot_FAR.emp)
#'  confint(boot_FAR.emp)
#'
#'  ylim <- range(boundFrechet, boot_FAR.emp)
#'  plot(boot_FAR.emp, ylim = ylim, main = "boot FAR empirical")
#'  # Theoretical FAR for in this case (Z = sigmaF * X  with X ~ Frechet)
#'  lines(rp, frechet_FAR(u = u, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
confint.boot_FARfit <- function(object, parm = dimnames(object)$u, level = 0.95, ...){
  irow <- dimnames(object)$u %in% parm
  alpha <- (1 - level) / 2
  ic <- t(apply(object[irow, ], 1, quantile, probs = c(alpha, 1 - alpha)))
  return(ic)
}

#' @rdname boot_FAR_fit.emp
#' @export
print.FARfit <- function(x, ...){
  FARdf <- data.frame(u = get_u(x),
                      FAR_hat = get_FAR(x),
                      sigma_FAR_hat = get_sigma_FAR(x))
  rownames(FARdf) <- NULL
  print("return periods, estimated FAR and estimated standard deviation for each return period:")
  print(FARdf)
}


#' @param ... additional arguments for the plot.
#' @rdname boot_FAR_fit.emp
#' @export
plot.FARfit <-function(x, ...){
  FAR <- get_FAR(x)
  u <- get_u(x)
  #
  # quantile for two side 90 confidence interval
  #
  ci90 <- confint(x, names(FAR), level = 0.90)
  ellipsis <- list(...)
  if(is.null(ellipsis$ylim)){
    ymin <- min(FAR,  ci90, na.rm = TRUE)
    ymax <- max(FAR,  ci90, na.rm = TRUE)
    ylim <- c(ymin, ymax)
  } else {
    ylim = ellipsis$ylim
    ellipsis$ylim <- NULL
  }
  if(is.null(ellipsis$col)){
    col = rgb(red = 0.5, green = .1, blue = 0.3, alpha = .5)
  } else {
    col = ellipsis$col
    ellipsis$col <- NULL
  }
  lu <- length(u)
  args_plot <- list(x =  u,
                    y = FAR,
                    xlab = "Return level u",
                    ylab = "FAR(u)",
                    type = "n",
                    ylim = ylim)
  args_plot <- c(args_plot, ellipsis)
  do.call(plot, args_plot)
  grid(lwd = 3)
  args_polygon <- list(x = c(u, u[lu:1]),
                       y = c(ci90[, 1], ci90[lu:1, 2]),
                       col = col, border=F)
  args_polygon <- c(args_polygon, ellipsis)
  do.call(polygon, args_polygon)
  args_lines <- list(x =  rp,
                     y = farr,
                     col = "black")
  args_lines <- c(args_lines, ellipsis)
  do.call(lines, args_lines)
  # title("Exponential FAR(r)")
}

get_FAR.boot_FARfit <- function(object, ...){
  FAR <- apply(object, 1, mean)
}

get_u.boot_FARfit <- function(object, ...){
  u <- dimnames(object)$u
}

get_sigma_FAR.boot_FARfit <- function(object, ...){
  FAR <- apply(object, 1, sd)
}

get_FAR <- function(object, ...){
  UseMethod("get_FAR")
}

get_u <- function(object, ...){
  UseMethod("get_u")
}

get_sigma_FAR <- function(object, ...){
  UseMethod("get_sigma_FAR")
}
