library(evd)
library(farr)

muF <-  1; xiF <- .15; sigmaF <-  1.412538 #  cst^(-xiF) # .05^(-xi1);
# asymptotic limit
boundFrechet <- 1 - sigmaF^(-1/xiF)
# sample size
size=100
# level=.9
set.seed(4)
z = rgev(size, loc = (sigmaF), scale = xiF * sigmaF, shape = xiF)
x = rgev(length(z), loc=(1), scale = xiF, shape=xiF)
rp = seq(from = 2, to = 30, length = 200)
u  = qgev(1 - (1 / rp),loc = muF, scale = xiF, shape = xiF)

# parametric estimation of far with exponential distribution for theta
theta_fit <- estim_theta.wexp(x = x, z = z) 
hist(theta_fit)
ecdf(theta_fit)
qqplot(theta_fit)
farr_fit.exp <- estim_farr.wexp(theta_hat = theta_fit$theta_hat, sigma_theta_hat = theta_fit$sigma_theta_hat, rp = rp) 
plot(farr_fit.exp)
print(farr_fit.exp)
theta_boot.exp <- boot_theta_fit.wexp(x = x, z = z, B = 10)
boot_farr.exp <- boot_farr_fit.wexp(theta_boot = theta_boot.exp$theta_boot , rp = rp)
confint(boot_farr.exp)
plot(boot_farr.exp)
print(boot_farr.exp)

# non-parametric estimation of far
farr_fit.emp <- estim_farr.emp(x = x, z = z, rp = rp) 
plot(farr_fit.emp)
print(farr_fit.emp)
boot_farr.emp <- boot_farr_fit.emp(x = x, z = z, rp = rp, B = 10)
confint(boot_farr.emp)
plot(boot_farr.emp)
print(boot_farr.emp)


# non-parametric estimation of FAR
# FAR_fit.emp <- estim_FAR.emp(x = x, z = z, u = u) 
boot_FAR.emp <- boot_FAR_fit.emp(x = x, z = z, u = u, B = 10)
confint(boot_FAR.emp)
plot(boot_FAR.emp)
print(boot_FAR.emp)

ylim <- range(boundFrechet, farr_fit.exp$farr_hat, farr_fit.emp$farr_hat)
par(mfrow = c(1, 5))
plot(farr_fit.exp, ylim = ylim, main = "far exponential")
lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
abline(h = boundFrechet, col = "red", lty = 2)

plot(boot_farr.exp, ylim = ylim, main = "boot far exponential",
     col = rgb(red = 0.3, green = 0.2, blue = 0.3, alpha = .5))
lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
abline(h = boundFrechet, col = "red", lty = 2)

plot(farr_fit.emp, ylim = ylim, main = "far non-parametric", 
     col = rgb(red = 0.5, green = 0.1, blue = 0.5, alpha = .5))
lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
abline(h = boundFrechet, col = "red", lty = 2)

plot(boot_farr.emp, ylim = ylim, main = "boot far non-parametric",
     col = rgb(red = 0.1, green = 0.5, blue = 0.5, alpha = .5))
lines(rp, frechet_farr(r = rp, sigma = sigmaF, xi = xiF), col = "red", lty = 2)


plot(boot_FAR.emp, ylim = ylim, main = "boot FAR non-parametric",
     col = rgb(red = 0.1, green = 0.3, blue = 0.7, alpha = .5))
lines(rp, frechet_FAR(u = u, sigma = sigmaF, xi = xiF), col = "red", lty = 2)
abline(h = boundFrechet, col = "red", lty = 2)

