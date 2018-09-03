### Gumble case
library(evd)
library(magrittr)
library(devtools)
devtools::load_all()
# RNG
set.seed(1)

# sample size
size <-  250 * 40


# Gumbel
# Simulations
rp <- 50
t <- seq.int(size)/size
mu = (seq(0, 2, length.out = size))^2
x = rgev(size * 4, loc = 0, scale = 1, shape = 0)
z = rgev(size, loc = mu, scale = 1, shape = 0)
theta_theo <- 1 / exp(mu)
p12_theo <- 1 / (1 + theta_theo)
G_theo <- function(z) pgev(z, loc = 0, scale = 1, shape = 0)
far_theo_rp <- (1 - theta_theo) * (1 - 1/rp)

kernel_gauss <- function(x, h = 1){
  dnorm(x / h) / h
}
attr(kernel_gauss, "K2_integrated") <-  1 / (2 * sqrt(pi))

kernel_epanechnikov <- function(x, h = 1){
  3/7 * (1 - (x/h)^2) / h
}
attr(kernel_epanechnikov, "K2_integrated") <-  3/5


dmat <- t %>% matrix(ncol = 1) %>%
  dist(method = "euclidian") %>%
  as.matrix()

kmat <- apply(dmat, 2, kernel_gauss, h = 25)


estim_p12 <- function(x, t, z, kernel, h = length(t)^(-1/5)){
  n <- length(t)
  m <- length(x)
  dmat <- t %>% matrix(ncol = 1) %>%
    dist(method = "euclidian") %>%
    as.matrix()
  kmat <- apply(dmat, 2, kernel_gauss, h = h)
  Gm <- ecdf(x)
  # Gm <- G_theo
  GmZ <- Gm(z)
  p12 <- apply(kmat, 2, function(weight){
    weighted.mean(GmZ, weight)
  })
  sig2 <- apply(kmat, 2, function(weight){
    weighted.mean((p12 - GmZ)^2, weight)
  })
  print(summary((1/m * p12 + (1 - 1/m) * sig2 ) / sig2))
  sig2 <- 1/m * p12 + (1 - 1/m) * sig2
  # K22 <- integrate(function(x) kernel(x)^2, lower = -Inf, upper = Inf)$value
  if( !is.null(attr(kernel, "K2_integrated"))){
    K22 <- attr(kernel, "K2_integrated")
  } else{
    K22 <- integrate(function(x) kernel(x)^2, lower = -Inf, upper = Inf)$value
  }
  ft <- apply(kmat, 2, mean)
  p12_var <- 1/(n * h) * sig2/ft * K22
  list(p12hat = p12, p12hat_sd = sqrt(p12_var), GmZ = GmZ)
}
p12 <- estim_p12(x = x, t = t, z = z, kernel = kernel_gauss, h = length(t)^(-1/5)/5)
p12 <- estim_p12(x = x, t = t, z = z, kernel = kernel_epanechnikov, h = length(t)^(-1/5)/5)
with(p12, plot(t, GmZ, pch = 20, col = "lightgrey"))
with(p12, lines(t, p12hat, lty = 2))
with(p12, lines(t, p12hat + 1.96 * p12hat_sd))
with(p12, lines(t, p12hat - 1.96 * p12hat_sd))
lines(t, p12_theo, col = "red")

xylim <- with(p12,
              range(p12_theo,
                    p12hat + 1.96 * p12hat_sd,
                    p12hat - 1.96 * p12hat_sd)
)
with(p12, plot(p12_theo, p12hat,
               type = "l", lty = 2,
               xlim = xylim, ylim = xylim))
with(p12, lines(p12_theo, p12hat + 1.96 * p12hat_sd))
with(p12, lines(p12_theo, p12hat - 1.96 * p12hat_sd))
abline(a = 0, b = 1, col = "red")

estim_theta <-function(p12hat, p12hat_sd) {
  thetahat <- 1/p12hat - 1
  thetahat_var <- p12hat_sd^2 / p12hat^4
  list(thetahat = thetahat, thetahat_sd = sqrt(thetahat_var))
}
theta <- with(p12, estim_theta(p12hat = p12hat, p12hat_sd = p12hat_sd))
with(theta, plot(t, 1/p12$GmZ - 1))
with(theta, lines(t, thetahat, lty = 2))
with(theta, lines(t, thetahat + 1.96 * thetahat_sd))
with(theta, lines(t, thetahat - 1.96 * thetahat_sd))
lines(t, theta_theo, col = "red")

xylim <- with(theta,
             range(theta_theo,
                   thetahat + 1.96 * thetahat_sd,
                   thetahat - 1.96 * thetahat_sd)
             )
with(theta, plot(theta_theo, thetahat,
                 type = "l", lty = 2,
                 xlim = xylim, ylim = xylim))
with(theta, lines(theta_theo, thetahat + 1.96 * thetahat_sd))
with(theta, lines(theta_theo, thetahat - 1.96 * thetahat_sd))
abline(a = 0, b = 1, col = "red")

estim_far <-function(thetahat, thetahat_sd, rp) {
  farhat <- (1 - thetahat) * (1 - 1/rp)
  farhat_var <- (1 - 1/rp)^2 * thetahat_sd^2
  list(farhat = farhat, farhat_sd = sqrt(farhat_var))
}
far_rp <- with(theta, estim_far(thetahat = thetahat,
                                thetahat_sd = thetahat_sd,
                                rp = rp))
with(theta, plot(t, (1 - (1/p12$GmZ - 1)) * (1 - 1/rp)))
with(far_rp, lines(t, farhat, lty = 2))
with(far_rp, lines(t, farhat + 1.96 * farhat_sd))
with(far_rp, lines(t, farhat - 1.96 * farhat_sd))
lines(t, far_theo_rp, col = "red")

xylim <- with(far_rp,
              range(far_theo_rp,
                    farhat + 1.96 * farhat_sd,
                    farhat - 1.96 * farhat_sd)
)
with(far_rp, plot(far_theo_rp, farhat,
                 type = "l", lty = 2,
                 xlim = xylim, ylim = xylim))
with(far_rp, lines(far_theo_rp, farhat + 1.96 * farhat_sd))
with(far_rp, lines(far_theo_rp, farhat - 1.96 * farhat_sd))
abline(a = 0, b = 1, col = "red")


################################################################################
# Test de distribution exponentielle de W

normalize_W <- function(W, theta){
  exp(log(W) - log(theta))
}

# In theory, OK !
W_theo <- - log(pgev(z, loc = 0, scale = 1, shape = 0))
W_norm_theo <- normalize_W(W_theo, theta_theo)
U_norm_theo <- 1 - exp(-W_norm_theo)
hist(U_norm_theo, breaks = 100, freq = FALSE)
farr::CoxOakes(x = W_norm_theo)
ks.test(x = W_norm_theo, y = "pexp", rate = 1)

W <- -log(farr::estim_GZ.kernel(x, z))
plot(W_theo, W)
abline(a = 0, b = 1, col = "red")

replace_zeros <- function(W_norm){
  izero <- which(W_norm == 0)
  if (length(izero) == 0) return(W_norm)
  minW <- min(W_norm[-izero])
  W_norm[izero] <- runif(length(izero), min = 0, max = minW)
  return(W_norm)
}
replace_ones <- function(Gz){
  ione <- which(Gz >= 1)
  if (length(ione) == 0) return(Gz)
  maxW <- max(Gz[-ione])
  Gz[ione] <- runif(length(ione), min = maxW, max = 1 - 1E-16)
  return(Gz)
}

W_norm_with_theta_theo <- normalize_W(W, theta_theo)
W_norm_with_theta_theo <- replace_zeros(W_norm_with_theta_theo)
hist(W_norm_with_theta_theo, breaks = 100, freq = FALSE)
lines(seq(0,10, 0.01), dexp(seq(0, 10, 0.01)))
plot(W_norm_theo, W_norm_with_theta_theo )
abline(a = 0, b = 1, col = "red")
farr::CoxOakes(x = W_norm_with_theta_theo)
ks.test(x = W_norm_with_theta_theo, y = "pexp", rate = 1)

W_normalized <- normalize_W(W, theta$thetahat)
W_normalized <- replace_zeros(W_normalized)
hist(W_normalized, breaks = 100, freq = FALSE)
lines(seq(0,10, 0.01), dexp(seq(0, 10, 0.01)))
U_normalized <- 1 - exp(-W_normalized)
hist(U_normalized, breaks = 100, freq = FALSE)
plot(W_norm_theo, W_norm_with_theta_theo )
abline(a = 0, b = 1, col = "red")
farr::CoxOakes(x = W_normalized)
ks.test(x = W_normalized, y = "pexp", rate = 1)
ks.test(x = U_normalized, y = "punif")

compute_Ustat <- function(Gz, theta){
  Gz <- replace_ones(Gz)
  hist(Gz, breaks = 1000)
  print(summary(Gz))
  print(sum(Gz == 1))
  1 - Gz^(1 / theta)
}
Ustat_theo <- compute_Ustat(Gz = pgev(z, loc = 0, scale = 1, shape = 0),
                            theta = theta_theo)
hist(Ustat_theo, breaks = 100)
ks.test(x = Ustat_theo, y = "punif")

Ustat_norm_with_theta_theo <- compute_Ustat(Gz = farr::estim_GZ.kernel(x, z),
                                            theta = theta_theo)
hist(Ustat_norm_with_theta_theo, breaks = 100)
ks.test(x = Ustat_norm_with_theta_theo, y = "punif")

Ustat_normalized <- compute_Ustat(Gz = farr::estim_GZ.kernel(x, z),
                                            theta = theta)
hist(Ustat_normalized, breaks = 100)
ks.test(x = Ustat_normalized, y = "punif")

estim_farr.ns <- function(theta, rp){
  rp_mat <- matrix(rp, ncol = length(rp), nrow = length(theta), byrow = TRUE)
  theta_mat <- matrix(theta, ncol = length(rp), nrow = length(theta))
  farr_mat <- farr::estim_farr_r(theta = theta_mat, rp = rp_mat)
  dimnames(farr_mat) <- list(theta = theta, rp = rp)
  return(farr_mat)
}
nsfar <- estim_farr.ns(theta = theta, rp = 2:10)
nsfar_theo <- estim_farr.ns(theta = theta_theo, rp = 2:10)
par(mfrow = c(3, 1))
ylim <- range(nsfar, nsfar_theo)
matplot(nsfar, pch = paste(2:10), col = rainbow(9), ylim = ylim)
matplot(nsfar_theo, pch = paste(2:10), col = rainbow(9), ylim = ylim)
matplot(nsfar - nsfar_theo, pch = paste(2:10), col = rainbow(9))


choose_h_for_wexp <- function(x = x, t = t, z = z,
                              kernel = kernel_epanechnikov){
  ftomin <- function(h){
    p12 <- estim_p12(x = x, t = t, z = z, kernel = kernel, h = h)
    theta <- with(p12, estim_theta(p12hat = p12hat, p12hat_sd = p12hat_sd))
    W_normalized <- normalize_W(W, theta$thetahat)
    W_normalized <- replace_zeros(W_normalized)
    pvalue <- farr::CoxOakes(x = W_normalized)$p.value
    return(1 - pvalue)
  }
  roughly <-  length(t)^(-1/5)
  optimize(ftomin, lower = 0, upper = max(t),
           tol = .Machine$double.eps^(1/10))
}
h_chosen <- choose_h_for_wexp(x = x, t = t, z = z,
                  kernel = kernel_epanechnikov)
p12 <- estim_p12(x = x, t = t, z = z, kernel = kernel_epanechnikov, h = h_chosen$minimum)
with(p12, plot(t, GmZ))
with(p12, lines(t, p12hat, lty = 2))
with(p12, lines(t, p12hat + 1.96 * p12hat_sd))
with(p12, lines(t, p12hat - 1.96 * p12hat_sd))
lines(t, p12_theo, col = "red")

cv_loo_h <- function(x = x, t = t, z = z,
                   kernel = kernel_epanechnikov, h = h){
  dmat <- t %>% matrix(ncol = 1) %>%
    dist(method = "euclidian") %>%
    as.matrix()
  kmat <- apply(dmat, 2, kernel_gauss, h = h)
  Gm <- ecdf(x)
  GmZ <- Gm(z)
  ft <- apply(kmat, 2, mean)
  compute_res2_i <- function(i){
    p12i <- weighted.mean(GmZ[-i], w = kmat[-i, i])
    e2i <- (GmZ[i] - p12i)^2
  }
  e2 <- vapply(seq_along(t),
               FUN = compute_res2_i,
               FUN.VALUE = numeric(1))
  cvh <- weighted.mean(e2, ft)
  return(cvh)
}
htotest <- seq(.Machine$double.eps, max(t), by = length(t)^(-1/5)/10)
cv <- sapply(htotest,
        cv_loo_h, x = x, t = t, z = z,
       kernel = kernel_epanechnikov)
plot(htotest, cv)
hchosen <- htotest[which.min(cv)]
hchosen <- optimize(cv_loo_h, interval = range(htotest),
                    x = x, t = t, z = z,
                    kernel = kernel_epanechnikov)$minimum

p12 <- estim_p12(x = x, t = t, z = z, kernel = kernel_epanechnikov, h = hchosen)
with(p12, plot(t, GmZ, pch = 20, col = "lightgrey"))
with(p12, lines(t, p12hat, lty = 2))
with(p12, lines(t, p12hat + 1.96 * p12hat_sd))
with(p12, lines(t, p12hat - 1.96 * p12hat_sd))
lines(t, p12_theo, col = "red")

