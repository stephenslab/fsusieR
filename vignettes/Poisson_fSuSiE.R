## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(fsusieR)
library(susieR)
library(ebnm)

set.seed(1)
"%!in%" <- function(x, y) !(x %in% y)

data(N3finemapping)
X_raw <- N3finemapping$X
N <- 200
mysd <- 1.1
genotype <- X_raw[sample(1:nrow(X_raw), size = N), 1:100]
idx_zero_var <- which(apply(genotype, 2, var) < 1e-15)
if (length(idx_zero_var) > 0) genotype <- genotype[, -idx_zero_var]
X <- (genotype - 0.99 * min(genotype)) / (0.5 * max(genotype))
G <- X

lev_res <- 6
L <- 2
lf <- list()
lf[[1]] <- rep(0, 2^lev_res); lf[[1]][20:30] <- 1.5
lf[[2]] <- rep(0, 2^lev_res); lf[[2]][50:60] <- 1.5

true_pos <- sample(1:ncol(G), replace = FALSE, size = L)

plot(lf[[1]], type = "l", main = paste("Effect of SNP", true_pos[1]), ylab = "Effect", xlab = "Position")
plot(lf[[2]], type = "l", main = paste("Effect of SNP", true_pos[2]), ylab = "Effect", xlab = "Position")

count.data <- vector("list", N)
True_intensity <- list()
for (i in 1:N) {
  predictor <- rep(0, 2^lev_res)
  for (l in 1:L) {
    G[, true_pos[l]] <- G[, true_pos[l]] - min(G[, true_pos[l]])
    predictor <- predictor + G[i, true_pos[l]] * lf[[l]] + 0.3
  }
  lambda <- exp(predictor + rnorm(length(predictor), sd = mysd))
  count.data[[i]] <- rpois(n = length(predictor), lambda = lambda)
  True_intensity[[i]] <- log(lambda)
}
Y <- do.call(rbind, count.data)
True_intensity <- do.call(rbind, True_intensity)

## -----------------------------------------------------------------------------
res0 <- susiF(Y = log1p(Y), X = X, L = 3)
res1 <- susiF(Y = HFT(Y), X = X, L = 3)
res <- Pois_fSuSiE(Y = Y, X = X, L = 3 , post_processing = "smash")

## -----------------------------------------------------------------------------
true_pos
res0$cs
res1$cs
res$susiF.obj$cs

## -----------------------------------------------------------------------------
plot(res0$ind_fitted_func, True_intensity, pch = 19, xlab = "Estimated Intensity", ylab = "True Intensity", main = "Comparison of Estimated Intensities")
points(res1$ind_fitted_func, True_intensity, col = "green", pch = 19)
points(res$susiF.obj$ind_fitted_func, True_intensity, col = "blue", pch = 19)
abline(a = 0, b = 1, lty = 2)
legend("bottomright", legend = c("log1p", "Haar-Fisz", "Poisson-fSuSiE"), col = c("black", "green", "blue"), pch = 19)

## -----------------------------------------------------------------------------
mean((c(res1$ind_fitted_func) - c(True_intensity))^2)
mean((c(res0$ind_fitted_func) - c(True_intensity))^2)
mean((c(res$susiF.obj$ind_fitted_func) - c(True_intensity))^2)

## -----------------------------------------------------------------------------

plot(lf[[1]], lwd = 2, lty = 2, type = "l", main = "Recovered Function 1", ylab = "Effect", xlab = "Position")
lines(res0$fitted_func[[1]], col = "black")
lines(res1$fitted_func[[1]], col = "green")
lines(res$susiF.obj$fitted_func[[1]], col = "blue")
legend("topright", legend = c("True", "log1p", "Haar-Fisz", "Poisson-fSuSiE"), col = c("black", "black", "green", "blue"), lty = c(2, 1, 1, 1))

plot(lf[[2]], lwd = 2, lty = 2, type = "l", main = "Recovered Function 2", ylab = "Effect", xlab = "Position")
lines(res0$fitted_func[[2]], col = "black")
lines(res1$fitted_func[[2]], col = "green")
lines(res$susiF.obj$fitted_func[[2]], col = "blue")
legend("topleft", legend = c("True", "log1p", "Haar-Fisz", "Poisson-fSuSiE"), col = c("black", "black", "green", "blue"), lty = c(2, 1, 1, 1))

## -----------------------------------------------------------------------------
res_2iter <- Pois_fSuSiE(Y = Y, X = X, L = 3, max.iter = 2, post_processing = "smash")
res_2iter$susiF.obj$cs

plot(res_2iter$susiF.obj$ind_fitted_func, True_intensity, col = "lightblue", pch = 19, xlab = "Estimated Intensity", ylab = "True Intensity")
points(res$susiF.obj$ind_fitted_func, True_intensity, col = "blue", pch = 19)
abline(a = 0, b = 1)

## -----------------------------------------------------------------------------
mean ( (res_2iter$susiF.obj$fitted_func[[2]] -  lf[[2]])^2)
mean ( (res $susiF.obj$fitted_func[[2]] -  lf[[2]])^2)

## -----------------------------------------------------------------------------

plot(lf[[1]], lwd = 2, lty = 2, type = "l", main = "Effect Function 1 (2 vs 20 iterations)")
lines(res_2iter$susiF.obj$fitted_func[[1]], col = "lightblue")
lines(res$susiF.obj$fitted_func[[1]], col = "blue")
legend("topright", legend = c("True", "2-iter", "10-iter"), col = c("black", "lightblue", "blue"), lty = c(2, 1, 1))

plot(lf[[2]], lwd = 2, lty = 2, type = "l", main = "Effect Function 2 (2 vs 20 iterations)")
lines(res_2iter$susiF.obj$fitted_func[[2]], col = "lightblue")
lines(res$susiF.obj$fitted_func[[2]], col = "blue")
legend("topright", legend = c("True", "2-iter", "10-iter"), col = c("black", "lightblue", "blue"), lty = c(2, 1, 1))


## -----------------------------------------------------------------------------

res <- Pois_fSuSiE(Y = Y, X = X, L = 3 ,max.iter=20, post_processing = "HMM", print=TRUE )

