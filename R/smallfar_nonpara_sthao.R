################################################################
###
###        non parametric estimation of p_{1,r}
###       with its stdev lambda_r
###
################################################################
min_FUN<-function(u,v){
   #  min(u,v)
  ((u + v) - abs(u - v)) / 2
}

outer_sthao <- function(u, v, FUN){
  umat <- matrix(u, nrow = length(u), ncol = length(v))
  vmat <- matrix(v, nrow = length(u), ncol = length(v), byrow = TRUE)
  FUN(umat, vmat)
}

mGZ_hat <-function(u, rp){
  #
  # Arguments u= a vector of reals on [0,1] and RP = a vector of integers greater than 2
  #
  # Value = empirical mean of E(U^(r-2)V^(r-2) min(U,V))
  #
  if(min(rp) < 2){
    stop("the integer r has to be greater than 2")
  }
  l <- length(u)
  # ll <- (l * (l - 1)) / 2
  triangle_MAT <- matrix(0, l, l)
  triangle_MAT[lower.tri(triangle_MAT)] <- 1
  ll <- sum(triangle_MAT)

  min_MAT <- outer(u, u ,FUN = min_FUN)
  prod_MAT<-outer(u, u, FUN = function(u,v) u*v)

  # min.MAT <-  outer_sthao(u, u, FUN = min.FUN)
  # prod.MAT <- outer_sthao(u, u, FUN=function(u,v){u*v})

  out <- vapply(rp,
                FUN = function(r){
                  r_MAT <- ((prod_MAT)^(r-2))*min_MAT*triangle_MAT
                  return(sum(r_MAT)/ll)
                },
                FUN.VALUE = numeric(1)
  )

  return(out)
}

non_para_pr <- function(x, z, rp){

  if(min(rp) < 2){
    stop("the length of r has to be greater than 2")
  }
  n <- length(z)
  GZ <- ecdf(x)(z)
  # estimates of p_{1, r} and p_{1, 2r-1}
  p1_r <- p1_2rm1 <- rep(NA, length(rp))
  for(i in seq_along(rp)){
    p1_r[i] <- mean(GZ^(rp[i] - 1))
    p1_2rm1[i] <- mean(GZ^(2 * rp[i] - 2))
  }

  # Empirical mean of E(U^(r-2)V^(r-2) min(U,V))
  m_vec <- mGZ_hat(GZ, rp)

  # stdev of p_{1,r}
  lambda_r <- p1_2rm1 - (p1_r)^2 + ((rp - 1)^2) * (m_vec - (p1_r)^2)
  if(any(is.na(lambda_r))){
    warning("Negative values in lambda2_r \n
            = > may require more data for a more precise estimation")
  }
  lambda_r <- sqrt(lambda_r)
  sigma_p1_r <- lambda_r / sqrt(n)
  p1_r_df <- data.frame(rp = rp, p1_r = p1_r, lambda_r = lambda_r, sigma_p1_r = sigma_p1_r)
  return(p1_r_df)
}

non_para_farr <- function(p1_r_df){
  #
  # Argument: a matrix obtained from the function non.para.pr
  #
  #
  # Value: a data.frame with three columns c("rp", "farr_hat", "sigma_farr_hat")
  #
  rp <- p1_r_df$rp
  p1_r <-  p1_r_df$p1_r
  lambda_r <- p1_r_df$lambda_r
  sigma_p1_r <- p1_r_df$sigma_p1_r

  farr_hat <- 1 - 1 / (rp * p1_r)
  names(farr_hat) <- paste("farr", rp, sep = "_")
  sigma_farr_hat <- sigma_p1_r / (rp * p1_r^2)
  farr_r_df <-  list(rp = rp, farr_hat = farr_hat, sigma_farr_hat = sigma_farr_hat)
  return(farr_r_df)
}

#' Non-parametric estimation of the far
#'
#' \code{estim_farr.np} returns an object of class  \code{("farrfit_np", "farrfit")} which contains
#' the results of the non-parametric estimation of the far.
#'
#' This function returns an non-parametric estimate of the far, the fraction of attributable risk for records,
#' as defined in Naveau et al (2018).
#'
#' For the full reference, see : \emph{Naveau, P., Ribes, A., Zwiers, F., Hannart, A., Tuel, A., & Yiou, P.
#' Revising return periods for record events in a climate event attribution context.
#' J. Clim., 2018.}, \url{https://doi.org/10.1175/JCLI-D-16-0752.1}
#'
#' @param x the variable of interest in the counterfactual world.
#' @param z the variable of interest in the factual world.
#' @param rp the return periods for which the far is to be estimated.
#'
#' @return An object of class \code{("farrfit_np", "farrfit")}. It is a list containing the following
#' elements:
#' \describe{
#'   \item{rp}{the return periods for which the far is estimated.}
#'   \item{farr_hat}{the estimate of the far for each return period \code{rp}}
#'   \item{sigma_farr_hat}{the standard deviation of the estimator of the far
#'   assuming asymptotic gaussianity}.
#' }
#'
#' @examples
#'  library(evd)
#'
#'  muF <-  1; xiF <- .15; sigmaF <-  1.412538 #  cst^(-xiF) # .05^(-xi1);
#'  # asymptotic limit for the farr in this case with a Frechet distributiom
#'  boundFrechet <- frechet_lim(sigma = sigmaF, xi = xiF)
#'  # sample size
#'  size <- 100
#'  # level=.9
#'  set.seed(4)
#'  z = rgev(size, loc = (sigmaF), scale = xiF * sigmaF, shape = xiF)
#'  x = rgev(length(z), loc=(1), scale = xiF, shape=xiF)
#'
#'  rp = seq(from = 2, to = 30, length = 200)
#'  # non-parametric estimation of far
#'  farr_fit.np <- estim_farr.np(x = x, z = z, rp = rp)
#'  print(farr_fit.np)
#'  ylim <- range(boundFrechet, farr_fit.np$farr_hat)
#'  plot(farr_fit.np, ylim = ylim, main = "far empirical")
#'  # Theoretical for in this case (Z = sigmaF * X  with X ~ Frechet)
#'  lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
estim_farr.np  <- function(x, z, rp){
  p1_r_df <- non_para_pr(x = x, z = z, rp = rp)
  farr_r_df <- non_para_farr(p1_r_df = p1_r_df)
  class(farr_r_df) <- c("farrfit_np", "farrfit")
  farr_r_df
}

#' @export
coef.farrfit_np <- function(object, ...){
  coef <- object$farr_hat
  names(coef) <- names(object$farr_hat)
  return(coef)
}

#' @export
vcov.farrfit_np <- function(object, ...){
  if(length(object$sigma_farr_hat) == 1){
    vcov <- matrix(object$sigma_farr_hat, nrow = 1, ncol = 1)
  } else{
    vcov <- diag(object$sigma_farr_hat)
  }
  rownames(vcov) <- colnames(vcov) <- names(object$farr_hat)
  return(vcov)
}

#' Non-parametric estimation of the far from bootstrap samples
#'
#' \code{boot_farr_fit.np} returns an object of class \code{("boot_farrfit.np", "boot_farrfit", "farrfit")}
#' which contains the estimates of the far for each bootstrap sample and
#' for different return periods \code{rp}
#'
#' This function returns bootstrap non-parametric estimates of far, the fraction of attributable risk for records,
#' as defined in Naveau et al (2018).The far is estimated from the bootstrap samples of the data x and z
#' that are obtained by resampling bootstrap.
#' The first bootstrap sample corresponds to the original dataset of x and z.
#'
#' For the full reference, see : \emph{Naveau, P., Ribes, A., Zwiers, F., Hannart, A., Tuel, A., & Yiou, P.
#' Revising return periods for record events in a climate event attribution context.
#' J. Clim., 2018.}, \url{https://doi.org/10.1175/JCLI-D-16-0752.1}
#'
#' @param x the variable of interest in the counterfactual world.
#' @param z the variable of interest in the factual world.
#' @param rp the return periods for which the far is to be estimated.
#' @param B the number of bootstrap samples to draw.
#'
#' @return An object of class \code{("boot_farrfit.np", "boot_farrfit", "farrfit")}.
#' It is a matrix where each columm contains the far estimated for the returns periods \code{rp}
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
#'
#'  # Resampling bootstrap for the non-parametrc estimation of the far
#'  boot_farr.np <- boot_farr_fit.np(x = x, z = z, rp = rp, B = 10)
#'  print(boot_farr.np)
#'  confint(boot_farr.np)
#'
#'  ylim <- range(boundFrechet, boot_farr.np)
#'  plot(boot_farr.np, ylim = ylim, main = "boot far non-parametric")
#'  # Theoretical far for in this case (Z = sigmaF * X  with X ~ Frechet)
#'  lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
#'  abline(h = boundFrechet, col = "red", lty = 2)
#' @export
boot_farr_fit.np <- function(x, z , rp, B = 100){
  xboot <- lapply(seq.int(B), function(b) sample(x, length(x), replace = TRUE))
  zboot <- lapply(seq.int(B), function(b) sample(z, length(z), replace = TRUE))
  xboot[[1]] <- x
  zboot[[1]] <- z
  farr_boot <- mapply(FUN = function(x, z){
    farr_fit <- estim_farr.np(x = x, z = z, rp = rp)
    return(farr_fit$farr_hat)
    }, xboot, zboot, SIMPLIFY = TRUE)
  dimnames(farr_boot) <- list(rp = rp, B = seq.int(B))
  class(farr_boot) <- c("boot_farrfit.np", "boot_farrfit", "farrfit")
  return(farr_boot)
}

