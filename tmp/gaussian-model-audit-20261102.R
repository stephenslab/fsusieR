# Independent numerical audit; does not alter the package or installed library.
pkgload::load_all("C:/Document/Serieux/Travail/Package/fsusieR", quiet = TRUE)
options(digits = 10)
lse <- function(x) { a <- max(x); a + log(sum(exp(x - a))) }
say <- function(label, x) cat(label, ": ", paste(signif(x, 10), collapse = ", "), "\n", sep = "")
log_mvn <- function(y, V) {
  U <- chol(V)
  z <- forwardsolve(t(U), y)
  -length(y) / 2 * log(2*pi) - sum(log(diag(U))) - sum(z*z) / 2
}

# Reference uses n-dimensional Gaussian marginal densities, independently
# of the package's regression-summary Bayes factor implementation.
reference_ser <- function(Y, X, sigma2, prior, index) {
  p <- ncol(X); T <- ncol(Y)
  lbf <- numeric(p); mu <- variance <- matrix(0, p, T)
  conditional_kl <- numeric(p)
  for (j in seq_len(p)) {
    x <- X[, j]; d <- sum(x*x)
    for (k in seq_len(T)) {
      group <- if (inherits(prior, "mixture_normal")) 1L else
        which(vapply(index, function(z) k %in% z, logical(1)))
      g <- prior[[group]]$fitted_g
      v0 <- g$sd^2
      logs <- vapply(v0, function(v)
        log_mvn(Y[, k], diag(sigma2, nrow(X)) + v * tcrossprod(x)), 0)
      mixlog <- lse(log(g$pi) + logs)
      lbf[j] <- lbf[j] + mixlog - log_mvn(Y[, k], diag(sigma2, nrow(X)))
      w <- exp(log(g$pi) + logs - mixlog)
      v <- sigma2 * v0 / (sigma2 + d*v0)
      m <- v0 * sum(x*Y[, k]) / (sigma2 + d*v0)
      mu[j, k] <- sum(w*m)
      variance[j, k] <- sum(w*(v+m*m)) - mu[j, k]^2
      nk <- numeric(length(v0))
      active <- v0 > 0
      nk[active] <- (log(v0[active]/v[active]) - 1 +
                      (v[active]+m[active]^2)/v0[active])/2
      conditional_kl[j] <- conditional_kl[j] +
        sum(w * (log(w) - log(g$pi) + nk))
    }
  }
  alpha <- exp(lbf - lse(lbf))
  kl <- sum(alpha * (log(alpha) + log(p) + conditional_kl))
  list(lbf=lbf, alpha=alpha, mu=mu, variance=variance, kl=kl)
}
as_obj <- function(posteriors, prior, sigma2) {
  structure(list(L=length(posteriors), G_prior=prior, sigma2=sigma2,
    alpha=lapply(posteriors, `[[`, "alpha"),
    fitted_wc=lapply(posteriors, `[[`, "mu"),
    fitted_wc2=lapply(posteriors, `[[`, "variance"),
    lBF=lapply(posteriors, `[[`, "lbf"),
    KL=vapply(posteriors, `[[`, 0, "kl")), class="susiF")
}
calc_kl <- function(obj, Y, X, index, l=1L)
  cal_KL_l(obj, l, X, Y[, -ncol(Y), drop=FALSE], Y[, ncol(Y)], index)
objective <- function(obj, Y, X, index)
  get_objective(obj, Y, X, Y[, -ncol(Y), drop=FALSE], Y[, ncol(Y)], index)

set.seed(20261102)
X <- colScale(matrix(rnorm(44), 11, 4))
Y <- colScale(outer(X[, 1], c(1.2, -.5, .7, .2, -.9, .4, 0, .3)) +
                matrix(rnorm(88), 11, 8), scale=FALSE)
s2 <- 1.3
index <- gen_wavelet_indx(3)
g <- ashr::normalmix(c(.55, .3, .15), rep(0, 3), c(0, .4, 1.2))
pooled <- structure(list(list(fitted_g=g)), class="mixture_normal")
per_scale <- structure(lapply(seq_along(index), function(i) {
  gi <- g; gi$sd <- gi$sd * (0.6 + .25*i)
  list(fitted_g=gi)
}), class="mixture_normal_per_scale")

for (prior in list(pooled, per_scale)) {
  cat("\nPRIOR: ", class(prior), "\n", sep="")
  ref <- reference_ser(Y, X, s2, prior, index)
  ss <- cal_Bhat_Shat(Y, X, sigma2=s2)
  bf <- log_BF(prior, ss$Bhat, ss$Shat, NULL, index)
  pm <- post_mat_mean(prior, ss$Bhat, ss$Shat, bf, NULL, index, e=0)
  pv <- post_mat_sd(prior, ss$Bhat, ss$Shat, bf, NULL, index, e=0)^2
  obj <- as_obj(list(ref), prior, s2)
  evidence <- sum(dnorm(Y, sd=sqrt(s2), log=TRUE)) + lse(ref$lbf) - log(ncol(X))
  say("Max log-BF error vs dense Gaussian integration", max(abs(bf-ref$lbf)))
  say("Max posterior mean/variance error", c(max(abs(pm-ref$mu)), max(abs(pv-ref$variance))))
  say("Immediate KL identity error", calc_kl(obj,Y,X,index)-ref$kl)
  say("Single-effect ELBO minus exact log evidence", objective(obj,Y,X,index)-evidence)
  stopifnot(max(abs(bf-ref$lbf)) < 1e-8,
            max(abs(pm-ref$mu)) < 1e-8, max(abs(pv-ref$variance)) < 1e-8,
            abs(calc_kl(obj,Y,X,index)-ref$kl) < 1e-8,
            abs(objective(obj,Y,X,index)-evidence) < 1e-8)

  # q and its prior are unchanged: their KL must not depend on sigma2.
  changed <- obj; changed$sigma2 <- 2*s2
  say("True KL, recomputed KL after changing sigma2 only",
      c(ref$kl, calc_kl(changed,Y,X,index)))
  wrong <- changed; wrong$KL <- calc_kl(changed,Y,X,index)
  say("ELBO error caused by refreshing KL after variance change",
      objective(wrong,Y,X,index)-objective(changed,Y,X,index))

  # After the next coordinate update, q1 and prior1 remain unchanged.
  residual2 <- Y-X %*% (ref$alpha*ref$mu)
  ref2 <- reference_ser(residual2, X, s2, prior, index)
  two <- as_obj(list(ref,ref2), prior,s2)
  say("Cached KL1, incorrectly refreshed KL1 after effect2 update",
      c(ref$kl,calc_kl(two,Y,X,index,1L)))
  say("KL2 identity error at its own update",calc_kl(two,Y,X,index,2L)-ref2$kl)

  # Global unit changes with the entire model transformed must preserve BF.
  factor <- 7
  scaled_prior <- prior
  for (i in seq_along(scaled_prior)) scaled_prior[[i]]$fitted_g$sd <-
    scaled_prior[[i]]$fitted_g$sd*factor
  scaled_ss <- cal_Bhat_Shat(Y*factor, X, sigma2=s2*factor^2)
  scaled_bf <- log_BF(scaled_prior,scaled_ss$Bhat,scaled_ss$Shat,NULL,index)
  say("Global unit-change log-BF error with matching prior/noise",max(abs(scaled_bf-bf)))
}

cat("\nTRANSFORMS\n")
wave <- function(y) { z <- DWT2(y); cbind(z$D,z$C) }
W <- wave(diag(8))
say("Orthonormal DWT error",max(abs(crossprod(W)-diag(8))))
position_scaled <- colScale(Y)
wave_scaled <- colScale(wave(position_scaled))
transform <- diag(1/attr(position_scaled,"scaled:scale")) %*% W %*%
  diag(1/attr(wave_scaled,"scaled:scale"))
V <- s2*crossprod(transform)
say("Noise variances after both observed-SD scalings",range(diag(V)))
say("Maximum absolute off-diagonal noise covariance",max(abs(V-diag(diag(V)))))
A <- remap_data(diag(7),1:7,verbose=FALSE)$Y
V_interp <- crossprod(A %*% W)
say("Interpolation: input positions, output positions, covariance rank",
    c(nrow(A),ncol(A),qr(V_interp)$rank))
say("Interpolation: max departure from identity noise covariance",
    max(abs(V_interp-diag(ncol(V_interp)))))

cat("\nAPPROXIMATION / MASKING / NUMERICS\n")
B <- matrix(.1,1001,8); S <- matrix(.3,1001,8); flat <- rep(0,1001)
means <- post_mat_mean(pooled,B,S,flat,NULL,index)
sds <- post_mat_sd(pooled,B,S,flat,NULL,index)
say("Diffuse alpha: max abs mean, max posterior SD",c(max(abs(means)),max(sds)))
say("Exact conditional SD for the same row",ashr::postsd(g,ashr::set_data(B[1,],S[1,]))[1])

# The likelihood ignores masked coefficients, but the M-step still sees them.
mask <- 1:4
masked_B <- matrix(rep(c(0,0,0,0,3,-2,.5,1),each=4),4,8)
masked_S <- matrix(1,4,8)
em_full <- EM_pi(pooled,masked_B,masked_S,index,max_step=5,
                control_mixsqp=list(verbose=FALSE),lowc_wc=mask,
                nullweight=1,max_SNP_EM=ncol(X),tol_null_prior=0)
em_active <- EM_pi(pooled,masked_B[,-mask],masked_S[,-mask],list(1:4),max_step=5,
                  control_mixsqp=list(verbose=FALSE),lowc_wc=NULL,
                  nullweight=1,max_SNP_EM=ncol(X),tol_null_prior=0)
say("Masked-vs-omitted coefficient mixture-weight difference",
    max(abs(em_full$tpi_k-em_active$tpi_k)))

extremeB <- matrix(c(100,0,rep(0,14)),2,8)
extremeS <- matrix(.1,2,8)
badbf <- suppressWarnings(log_BF(pooled,extremeB,extremeS,NULL,index))
stable <- vapply(1:2,function(j) sum(vapply(1:8,function(k)
  lse(log(g$pi)+dnorm(extremeB[j,k],sd=sqrt(g$sd^2+extremeS[j,k]^2),log=TRUE))-
    dnorm(extremeB[j,k],sd=extremeS[j,k],log=TRUE),0)),0)
say("Extreme data log-BF: current",badbf)
say("Extreme data log-BF: stable reference",stable)

cat("\nSINGLE-EFFECT CONVERGENCE\n")
set.seed(123)
X1 <- matrix(rnorm(30*6),30,6)
Y1 <- matrix(rnorm(30*16,sd=2),30,16)+outer(X1[,1],sin(1:16))
fit1 <- susiF(Y1,X1,L=1,standardize_Y=FALSE,post_processing="none",
              filter_cs=FALSE,cal_obj=TRUE,nullweight=1,verbose=FALSE)
say("L=1 number of objective evaluations",length(fit1$ELBO))

# Public e=0 must reach every posterior update, including L > 1.
seen_e <- numeric()
record_e <- function(e) seen_e <<- c(seen_e,e)
trace("update_susiF_obj.susiF",where=asNamespace("fsusieR"),
      tracer=quote(record_e(e)),print=FALSE)
fit2 <- susiF(Y1,X1,L=2,L_start=2,standardize_Y=FALSE,post_processing="none",
              filter_cs=FALSE,greedy=FALSE,backfit=FALSE,cal_obj=TRUE,
              nullweight=1,verbose=FALSE,e=0,maxit=4)
untrace("update_susiF_obj.susiF",where=asNamespace("fsusieR"))
say("Internal e values with public e=0, L=2",unique(seen_e))

# Independent check of the expected residual sum of squares.
ref1 <- reference_ser(Y,X,s2,pooled,index)
ref2 <- reference_ser(Y-X%*%(ref1$alpha*ref1$mu),X,s2,pooled,index)
ob <- as_obj(list(ref1,ref2),pooled,s2)
M <- lapply(list(ref1,ref2),function(z) z$alpha*z$mu)
Q <- lapply(list(ref1,ref2),function(z) z$alpha*(z$variance+z$mu^2))
erss <- sum((Y-X%*%Reduce(`+`,M))^2)+
  sum(vapply(1:2,function(l) sum(colSums(X^2)*Q[[l]])-sum((X%*%M[[l]])^2),0))
say("ERSS error on standardized X",get_ER2(ob,Y,X)-erss)
say("Residual variance M-step error",estimate_residual_variance(ob,Y,X)-erss/length(Y))

cat("Audit complete. Source code was not modified.\n")
