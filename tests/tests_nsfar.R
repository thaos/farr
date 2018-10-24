### Gumble case
library(evd)
library(magrittr)
library(devtools)
library(truncdist)
devtools::load_all()
source("tests/tests_nsfar_algo.R")
# RNG
set.seed(1)


##########################################################################
# Example 2
# Frechet
# Simulations
size <-  250 * 4
rp <- 50
t <- seq.int(size)/size
sigma <- seq(1, 2, length.out = size)
xi = 0.5
x = rgev(size * 5, loc = 1, scale = xi, shape = xi)
z = rgev(size, loc = sigma, scale = xi * sigma, shape = xi)
theta_theo <- sigma^(-1 / xi)
p12_theo <- 1 / (1 + theta_theo)
G_theo <- function(z) pgev(z, loc = 1, scale = xi, shape = xi)
far_theo_rp <- (1 - theta_theo) * (1 - 1/rp)


##########################################################################
# Example 1
# Gumbel
# Simulations
size <-  250 * 4
rp <- 50
t <- seq.int(size)/size
mu = 2 + exp(-seq(0, 2, length.out = size))
mu = 2 + seq(0, 5, length.out = size)
x = rgev(size * 5, loc = 0, scale = 1, shape = 0)
z = rgev(size, loc = mu, scale = 1, shape = 0)
theta_theo <- 1 / exp(mu)
p12_theo <- 1 / (1 + theta_theo)
G_theo <- function(z) pgev(z, loc = 0, scale = 1, shape = 0)
far_theo_rp <- (1 - theta_theo) * (1 - 1/rp)


# kernel choice:
# - kernel_epanechnikov
# - kernel_gauss
krnl <- kernel_epanechnikov

p12 <- estim_p12_ns(x = x, t = t, z = z, kernel = krnl, h = length(t)^(-1/5)/5)
with(p12, plot(t, GmZ, pch = 20, col = "lightgrey"))
with(p12, lines(t, p12_hat, lty = 2))
with(p12, lines(t, p12_hat + 1.96 * sigma_p12_hat))
with(p12, lines(t, p12_hat - 1.96 * sigma_p12_hat))
lines(t, p12_theo, col = "red")

xylim <- with(p12,
              range(p12_theo,
                    p12_hat + 1.96 * sigma_p12_hat,
                    p12_hat - 1.96 * sigma_p12_hat)
)
with(p12, plot(p12_theo, p12_hat,
               type = "l", lty = 2,
               xlim = xylim, ylim = xylim))
with(p12, lines(p12_theo, p12_hat + 1.96 * sigma_p12_hat))
with(p12, lines(p12_theo, p12_hat - 1.96 * sigma_p12_hat))
abline(a = 0, b = 1, col = "red")


theta <- with(p12, estim_theta_ns(p12_hat = p12_hat, sigma_p12_hat = sigma_p12_hat))
with(theta, plot(t, 1/p12$GmZ - 1))
with(theta, lines(t, theta_hat, lty = 2))
with(theta, lines(t, theta_hat + 1.96 * sigma_theta_hat))
with(theta, lines(t, theta_hat - 1.96 * sigma_theta_hat))
lines(t, theta_theo, col = "red")

xylim <- with(theta,
             range(theta_theo,
                   theta_hat + 1.96 * sigma_theta_hat,
                   theta_hat - 1.96 * sigma_theta_hat)
             )
with(theta, plot(theta_theo, theta_hat,
                 type = "l", lty = 2,
                 xlim = xylim, ylim = xylim))
with(theta, lines(theta_theo, theta_hat + 1.96 * sigma_theta_hat))
with(theta, lines(theta_theo, theta_hat - 1.96 * sigma_theta_hat))
abline(a = 0, b = 1, col = "red")

theta <- estim_theta.nswexp(x = x, t = t, z = z, kernel = krnl, h = NULL)
theta$utest_pvalue
coef(theta)
vcov(theta)
theta_ci95 <- confint(theta, level = c(0.95))
plot(theta)
lines(t, theta_theo, col = "red")
hist(theta)
qqplot(theta)
ecdf(theta)


far_rp <- with(theta, estim_farr_ns(theta_hat = theta_hat,
                                sigma_theta_hat = sigma_theta_hat,
                                rp = rp))
with(theta, plot(t, (1 - (1/p12$GmZ - 1)) * (1 - 1/rp),
                 pch = 20, col = "lightgrey"))
lines(t, far_rp[, 1], lty = 2)
lines(t, far_rp[, 1] + 1.96 * far_rp[, 2])
lines(t, far_rp[, 1] - 1.96 * far_rp[, 2])
lines(t, far_theo_rp, col = "red")

xylim <- range(far_theo_rp, far_rp[, 1] + 1.96 * far_rp[, 2], far_rp[, 1] - 1.96 * far_rp[, 2])
plot(far_theo_rp, far_rp[, 1],
     type = "l", lty = 2,
     xlim = xylim, ylim = xylim)
lines(far_theo_rp, far_rp[, 1] + 1.96 * far_rp[, 2])
lines(far_theo_rp, far_rp[, 1] - 1.96 * far_rp[, 2])
abline(a = 0, b = 1, col = "red")

rp_list <- seq(2, 10, length.out = 30)
far_rp <- with(theta, estim_farr.nswexp(theta, rp = rp_list))
coef(far_rp)
vcov(far_rp)
farr_ci <- confint(far_rp, level = 0.95)
plot(far_rp)
plot(extract_rp(far_rp, rp = head(rp_list, n = 1)))





################################################################################
# Test de distribution exponentielle de W


# In theory, OK !
W_theo <- - log(pgev(z, loc = 0, scale = 1, shape = 0))
W_norm_theo <- normalize_W(W_theo, theta_theo)
hist(W_norm_theo, breaks = 100, freq = FALSE)
lines(seq(0,10, 0.01), dexp(seq(0, 10, 0.01)))
farr::CoxOakes(x = W_norm_theo)
ks.test(x = W_norm_theo, y = "pexp", rate = 1)

W <- -log(farr::estim_GZ.kernel(x, z))
plot(W_theo, W)
abline(a = 0, b = 1, col = "red")

W_norm_with_theta_theo <- normalize_W(W, theta_theo)
W_norm_with_theta_theo <- replace_zeros(W_norm_with_theta_theo)
hist(W_norm_with_theta_theo, breaks = 100, freq = FALSE)
lines(seq(0,10, 0.01), dexp(seq(0, 10, 0.01)))
plot(W_norm_theo, W_norm_with_theta_theo )
abline(a = 0, b = 1, col = "red")
farr::CoxOakes(x = W_norm_with_theta_theo)
ks.test(x = W_norm_with_theta_theo, y = "pexp", rate = 1)

W_normalized <- normalize_W(W, theta$theta_hat)
W_normalized <- replace_zeros(W_normalized)
hist(W_normalized, breaks = 100, freq = FALSE)
lines(seq(0,10, 0.01), dexp(seq(0, 10, 0.01)))
plot(W_norm_theo, W_norm_with_theta_theo )
abline(a = 0, b = 1, col = "red")
farr::CoxOakes(x = W_normalized)
ks.test(x = W_normalized, y = "pexp", rate = 1)

hchosen <- choose_h_for_wexp(x = x, t = t, z = z,
                  kernel = krnl)$minimum
p12 <- estim_p12_ns(x = x, t = t, z = z, kernel = krnl, h = hchosen)
with(p12, plot(t, GmZ, pch = 20, col = "lightgrey"))
with(p12, lines(t, p12_hat, lty = 2))
with(p12, lines(t, p12_hat + 1.96 * sigma_p12_hat))
with(p12, lines(t, p12_hat - 1.96 * sigma_p12_hat))
lines(t, p12_theo, col = "red")

htotest <- seq(.Machine$double.eps, max(t)/3, by = length(t)^(-1/5)/30)[-1]
cv <- sapply(htotest,
             cv_loo_h, x = x, t = t, z = z,
             kernel = krnl)
plot(htotest, cv)
hchosen <- htotest[which.min(cv)]
# hchosen <- optimize(cv_loo_h, interval = range(htotest),
                    # x = x, t = t, z = z,
                    # kernel = krnl)$minimum

p12 <- estim_p12_ns(x = x, t = t, z = z, kernel = krnl, h = hchosen)
with(p12, plot(t, GmZ, pch = 20, col = "lightgrey"))
with(p12, lines(t, p12_hat, lty = 2))
with(p12, lines(t, p12_hat + 1.96 * sigma_p12_hat))
with(p12, lines(t, p12_hat - 1.96 * sigma_p12_hat))
lines(t, p12_theo, col = "red")

p12boot <- boot_p12(x = x, t = t, z = z,
                    kernel = krnl, h = hchosen,
                    B = 100)
plot(t, p12$GmZ, pch = 20, col = "lightgrey")
lines(t, apply(p12boot$p12_boot, 1, mean), lwd = 1, lty = 2)
lines(t, apply(p12boot$p12_boot, 1, quantile, probs = 0.975), lwd = 1, lty = 1)
lines(t, apply(p12boot$p12_boot, 1, quantile, probs = 0.025), lwd = 1, lty = 1)
lines(t, p12_theo, col = "red", lwd = 1)

thetaboot <- boot_theta(p12_boot = p12boot$p12_boot)
thetaboot_ci025 <- apply(thetaboot, 1, quantile, probs = 0.025)
thetaboot_ci975 <- apply(thetaboot, 1, quantile, probs = 0.975)
xylim <- range(theta_theo,
               thetaboot_ci025,
               thetaboot_ci975)
plot(theta_theo, apply(thetaboot, 1, mean),
                 type = "l", lty = 2,
                 xlim = xylim, ylim = xylim)
lines(theta_theo, thetaboot_ci025)
lines(theta_theo, thetaboot_ci975)
abline(a = 0, b = 1, col = "red")

thetaboot <- boot_theta_fit.nswexp(x = x, t = t, z = z,
                      kernel = krnl, h = hchosen,
                      B = 100)
thetaboot_ci025 <- apply(thetaboot$theta_boot, 1, quantile, probs = 0.025)
thetaboot_ci975 <- apply(thetaboot$theta_boot, 1, quantile, probs = 0.975)
xylim <-  range(theta_theo,
                thetaboot_ci025,
                thetaboot_ci975)
plot(theta_theo, apply(thetaboot$theta_boot, 1, mean),
     type = "l", lty = 2,
     xlim = xylim, ylim = xylim)
lines(theta_theo, thetaboot_ci025)
lines(theta_theo, thetaboot_ci975)
abline(a = 0, b = 1, col = "red")

bootfarr <- boot_farr_fit.nswexp(thetaboot, rp = rp)
bootfarr <- boot_farr_fit.nswexp(thetaboot, rp = seq(2, 10, length.out = 20))
confint(bootfarr)
plot(bootfarr)


tref <- runif(1)
size <- 100000
h = size^(-1/5)/100

samples_BBNW <- sapply(1:1000, function(i) drawBBNW(tref = tref, size = size, h = h))
mean(samples_BBNW^2)
var(samples_BBNW)

vsize <- seq(1000, 20000, by = 1000)
vvar <- sapply(vsize, varempirical)

# pairs(t(samples_BBNW))
mBBNW <- apply(samples_BBNW, 2, mean)
hist(mBBNW, breaks = 100)
plot(density(mBBNW))
summary(mBBNW)
ks.test(mBBNW, "pnorm", mean = 0, sd = sd(mBBNW))

mean(samples_BBNW^2)
mu_theo <- 2 + exp(tref)
theta_theo <- 1 / exp(mu_theo)
p12_theo <- 1 / (1 + theta_theo)
p12var_theo <- 1/(size * h) * p12_theo * (1 - p12_theo) * attr(kernel_gauss,"K2_integrated")

# non-stationnary exponential test
x_h0 <- rgev(length(x), loc = 0, scale = 1, shape = 0)
z_h0 <- rgev(length(z), loc = -log(theta_theo), scale = 1, shape = 0)
W_h0 <- -log(ecdf(x_h0)(z_h0))/theta_theo
W <- -log(ecdf(x)(z))/theta_theo

hist(rexp(length(z)), breaks = 100)
hist(W_h0, breaks = 100)
hist(W, breaks = 100)
qqplot(W_h0, W)
ks.test(W_h0, W)
Matching::ks.boot(W_h0, W, nboots = 10000)
ks.test(W, "dexp")

##############################################################################
#########################################################################
# Example 1
# Gumbel
# Simulations
size <-  250 * 4
rp <- 50
t <- seq.int(size)/size
mu = 2 + exp(-seq(0, 2, length.out = size))
mu = 2 + seq(0, 5, length.out = size)
x = rgev(size * 5, loc = 0, scale = 1, shape = 0)
z = rgev(size, loc = mu, scale = 1, shape = 0)
theta_theo <- 1 / exp(mu)
p12_theo <- 1 / (1 + theta_theo)
G_theo <- function(z) pgev(z, loc = 0, scale = 1, shape = 0)
far_theo_rp <- (1 - theta_theo) * (1 - 1/rp)

# Example 2
# Frechet
# Simulations
size <-  250 * 4
rp <- 50
t <- seq.int(size)/size
sigma <- seq(1, 2, length.out = size)
xi = 0.5
x = rgev(size * 5, loc = 1, scale = xi, shape = xi)
z = rgev(size, loc = sigma, scale = xi * sigma, shape = xi)
theta_theo <- sigma^(-1 / xi)
p12_theo <- 1 / (1 + theta_theo)
G_theo <- function(z) pgev(z, loc = 1, scale = xi, shape = xi)
far_theo_rp <- (1 - theta_theo) * (1 - 1/rp)




kstest_betatrunc(x, z)

ks.test(GZhat, pbeta, 1/theta_theo, 1)
ks.test(GZhat_trunc, pbeta, 1/theta_theo, 1)
W <- estim_W(GZhat)
W_normalized <- normalize_W(W, theta_theo)
dist_ks <- ks.test(W_normalized, "pexp", rate = 1)$statistic
dist_ks_h0 <- sapply(1:1000, function(i) gen_h0(x, z, theta_theo))
hist(dist_ks_h0, breaks = 50)
abline(v = dist_ks)
pvalue <- mean(dist_ks < dist_ks_h0)




check_pvalues <- function(){
  size <-  250 * 4
  t <- seq.int(size)/size
  mu = 2 + seq(0, 5, length.out = size)
  x = rgev(size * 5, loc = 0, scale = 1, shape = 0)
  z = rgev(size, loc = mu, scale = 1, shape = 0)
  theta_theo <- 1 / exp(mu)
  compute_exptest_pvalue(x, z, theta, n_h0 = 200)
}

pvalues <- sapply(1:200, function(i) check_pvalues())
hist(pvalues)
ks.test(pvalues, punif)

qqplot(W_h0_normalized, W_normalized)
qqplot(GZhat, GZhat_h0)
ks.test(W_h0_normalized, W_normalized)
ks.test(GZhat, GZhat_h0)
Matching::ks.boot(W_h0_normalized, W_normalized)
ks.test(W_h0_normalized, "pexp", rate = 1)
ks.test(W_normalized, "pexp", rate = 1)
hist(W_normalized, breaks = 100, freq = FALSE)
lines(seq(0,10, 0.01), dexp(seq(0, 10, 0.01)))
hist(W_h0_normalized, breaks = 100, freq = FALSE)
lines(seq(0,10, 0.01), dexp(seq(0, 10, 0.01)))

hist(Gd, breaks = 100, freq = FALSE)
lines(seq(0,10, 0.01), dexp(seq(0, 10, 0.01)))

