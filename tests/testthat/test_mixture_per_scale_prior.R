# tests/test_mixture_per_scale_prior.R

suppressPackageStartupMessages({
  library(testthat)
  library(ashr)
  library(wavethresh)
})
# mixsqp sometimes prints a build-version warning; silence it for clean tests
suppressWarnings(suppressPackageStartupMessages(library(mixsqp)))
suppressPackageStartupMessages(library(fsusieR))

skip_if_not_installed("ashr")
skip_if_not_installed("wavethresh")
skip_if_not_installed("mixsqp")
skip_if_not_installed("fsusieR")

set.seed(2)
control_mixsqp <- list(verbose = FALSE)

# --- Helper: safe column scaling (if not provided by your pkg) ---
if (!exists("colScale", mode = "function")) {
  colScale <- function(M, scale = TRUE, center = TRUE) {
    sc <- if (scale) apply(M, 2, sd) else rep(1, ncol(M))
    sc[sc == 0 | is.na(sc)] <- 1
    ct <- if (center) colMeans(M) else rep(0, ncol(M))
    sweep(sweep(M, 2, ct, "-"), 2, sc, "/")
  }
}

test_that("mixture-per-scale prior end-to-end sanity", {
  # --------- Data gen ----------
  f1 <- simu_IBSS_per_level(lev_res = 9, alpha = 1, prop_decay = 1.5)
  if (interactive()) {
    plot(f1$sim_func, type = "l", ylab = "y")
  }
  N <- 500
  P <- 10
  lfsr_curve <- 0.05
  nullweight <- 10 / sqrt(N)

  set.seed(23)
  G <- matrix(sample(c(0,1,2), size = N * P, replace = TRUE), nrow = N, ncol = P)
  beta0 <- 0
  beta1 <- 1
  pos1  <- 5
  rsnr  <- 10

  noisy.data <- lapply(seq_len(N), function(i) {
    f1_obs <- f1$sim_func
    beta1 * G[i, pos1] * f1_obs + rnorm(length(f1$sim_func), sd = (1 / rsnr) * sd(f1$sim_func))
  })
  noisy.data <- do.call(rbind, noisy.data)

  Y <- noisy.data
  X <- G
  Y <- colScale(Y, scale = FALSE)
  X <- colScale(X)

  W <- DWT2(Y)
  update_D <- W
  Y_f <- cbind(W$D, W$C) # Using a column-like phenotype
  update_Y <- Y_f

  v1 <- rep(1, ncol(X))
  tt <- cal_Bhat_Shat(Y_f, X, v1, lowc_wc = NULL)
  indx_lst <- gen_wavelet_indx(9)
  Bhat <- tt$Bhat
  Shat <- tt$Shat
  init_pi0_w <- 1

  # ---------- Prior ----------
  G_prior <- init_prior(
    Y = Y_f, X = X,
    prior = "mixture_normal_per_scale",
    v1 = v1, indx_lst = indx_lst,
    lowc_wc = NULL,
    control_mixsqp = control_mixsqp,
    nullweight = nullweight,
    max_SNP_EM = 100,
    tol_null_prior = 0
  )$G_prior

  # ---------- lBF localization ----------
  lBF <- log_BF(G_prior, Bhat, Shat, indx_lst, lowc_wc = NULL)
  expect_equal(which.max(lBF), pos1)

  # ---------- susiF object init ----------
  greedy  <- TRUE
  backfit <- TRUE

  susiF_obj <- init_susiF_obj(
    L_max = 2, G_prior, Y, X,
    L_start = 2, greedy = greedy, backfit = backfit
  )

  expect_equal(get_pi(susiF_obj, 1), get_pi_G_prior(G_prior))
  expect_equal(get_pi(susiF_obj, 2), get_pi_G_prior(G_prior))

  susiF_obj <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                              L_start = 4, greedy = greedy, backfit = backfit)
  susiF_obj <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                              L_start = 1, greedy = greedy, backfit = backfit)
  expect_equal(get_pi(susiF_obj, 1), get_pi_G_prior(G_prior))

  # ---------- Expansion ----------
  obj <- init_susiF_obj(L_max = 10, G_prior, Y, X,
                        L_start = 3, greedy = greedy, backfit = backfit)
  expect_equal(obj$L_max, 10)
  expect_equal(obj$L, 3)
  obj <- expand_susiF_obj(obj, L_extra = 7)
  expect_equal(obj$L_max, 10)
  expect_equal(obj$L, length(obj$fitted_wc))
  expect_equal(obj$L, length(obj$alpha))
  expect_equal(obj$L, length(obj$fitted_wc2))
  expect_equal(obj$L, length(obj$G_prior))
  expect_equal(obj$L, length(obj$cs))
  expect_equal(obj$L, length(obj$est_pi))
  expect_equal(obj$L, length(obj$est_sd))
  expect_equal(obj$L, length(obj$lBF))
  expect_equal(obj$L, length(obj$cred_band))

  # ---------- Internal prior ----------
  susiF_obj <- init_susiF_obj(L_max = 2, G_prior, Y, X,
                              L_start = 2, greedy = greedy, backfit = backfit)
  expect_equal(get_G_prior(susiF_obj), G_prior)

  # ---------- Class checks ----------
  expect_equal(
    class(init_prior(
      Y = Y_f, X = X,
      prior = "mixture_normal_per_scale",
      v1 = v1, indx_lst = indx_lst,
      lowc_wc = NULL,
      control_mixsqp = control_mixsqp,
      nullweight = nullweight,
      max_SNP_EM = 100,
      tol_null_prior = 0
    )$G_prior)[1],
    "mixture_normal_per_scale"
  )

  if (interactive()) {
    plot(Bhat, post_mat_mean(G_prior, Bhat, Shat, indx_lst = indx_lst, lowc_w = NULL))
    plot(Shat, post_mat_sd(G_prior, Bhat, Shat, indx_lst = indx_lst, lowc_w = NULL))
  }

  expect_equal(class(get_pi_G_prior(G_prior))[1], "pi_mixture_normal_per_scale")
  expect_equal(class(get_sd_G_prior(G_prior))[1], "sd_mixture_normal_per_scale")

  # ---------- Likelihood ----------
  L <- L_mixsq(G_prior, Bhat, Shat, indx_lst)
  expect_equal(class(L)[1], "lik_mixture_normal_per_scale")

  # ---------- Responsibilities ----------
  zeta <- cal_zeta(lBF)
  expect_equal(which.max(zeta), pos1)

  # ---------- M-step ----------
  tpi <- m_step(L, zeta, indx_lst,
                init_pi0_w = init_pi0_w,
                control_mixsqp = control_mixsqp,
                nullweight = nullweight,
                tol_null_prior = 0)
  expect_equal(class(tpi)[1], "pi_mixture_normal_per_scale")

  for (i in 1:8) {
    expect_gte(fsusieR::get_pi0(tpi = tpi)[i], c(0, 0.5, rep(1, 8))[i])
  }

  G_update <- update_prior(G_prior, tpi)
  expect_true(identical(get_pi_G_prior(G_update), tpi))

  # ---------- EM over pi ----------
  outEM <- EM_pi(G_prior, Bhat, Shat, indx_lst,
                 init_pi0_w = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 lowc_wc = NULL,
                 nullweight = nullweight,
                 tol_null_prior = 0)

  expect_equal(class(outEM$tpi_k)[1], "pi_mixture_normal_per_scale")
  expect_type(outEM$lBF, "double")
  expect_equal(length(outEM$lBF), nrow(Bhat))
  for (i in 1:8) {
    expect_gte(fsusieR::get_pi0(tpi = outEM$tpi_k)[i], c(0, 0.5, rep(1, 8))[i])
  }

  susiF_obj <- update_pi(susiF_obj, 1, outEM$tpi_k)
  expect_equal(get_pi(susiF_obj, 1), outEM$tpi_k)

  new_alpha <- cal_zeta(outEM$lBF)
  susiF_obj <- update_alpha(susiF_obj, 1, new_alpha)
  expect_equal(get_alpha(susiF_obj, 1), new_alpha)

  G_prior <- update_prior(G_prior, tpi = outEM$tpi_k)
  susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst)

  expect_equal(
    susiF_obj$fitted_wc[[1]],
    post_mat_mean(G_prior, Bhat, Shat, lBF = outEM$lBF, indx_lst = indx_lst, lowc_w = NULL)
  )
  expect_equal(
    susiF_obj$fitted_wc2[[1]],
    post_mat_sd(G_prior, Bhat, Shat, lBF = outEM$lBF, indx_lst = indx_lst, lowc_w = NULL)^2
  )
  expect_equal(get_alpha(susiF_obj, 1), cal_zeta(outEM$lBF))
  expect_equal(get_G_prior(susiF_obj), G_prior)
  expect_equal(susiF_obj$lBF[[1]], outEM$lBF)

  # ---------- Partial residual ----------
  update_T <- cal_partial_resid(
    obj = susiF_obj, l = 1, X = X, D = W$D, C = W$C, indx_lst = indx_lst
  )
  Lcur <- 2
  l    <- 1
  id_L <- (1:Lcur)[-((l %% Lcur) + 1)]
  manual_update_D <- W$D - Reduce("+", lapply(id_L, function(lid)
    (X * rep(susiF_obj$alpha[[lid]], rep.int(N, P))) %*%
      (susiF_obj$fitted_wc[[lid]][, -indx_lst[[length(indx_lst)]]])
  ))
  manual_update_C <- W$C - Reduce("+", lapply(id_L, function(lid)
    (X * rep(susiF_obj$alpha[[lid]], rep.int(N, P))) %*%
      susiF_obj$fitted_wc[[lid]][, indx_lst[[length(indx_lst)]]]
  ))
  manual_update <- cbind(manual_update_D, manual_update_C)
  expect_equal(update_T, manual_update)

  # ---------- Fit curves precision ----------
  susiF_obj <- update_susiF_obj(susiF_obj, 2, outEM, Bhat, Shat, indx_lst)
 # susiF_obj <- update_susiF_obj(susiF_obj, 3, outEM, Bhat, Shat, indx_lst)
  #susiF_obj <- update_susiF_obj(susiF_obj, 4, outEM, Bhat, Shat, indx_lst)
  susiF_obj$cs=susiF_obj$cs[1]
  susiF_obj$L=1
  # NOTE: fix the original typo: fitted_funcfitted_func -> fitted_func
  fit1 <- unlist(update_cal_fit_func(obj = susiF_obj, Y = Y, X = X,
                                     indx_lst = indx_lst, TI = FALSE)$fitted_func[[1]])
  expect_equal(sqrt(mean( (fit1 - f1$sim_func)^2)), 0, tol = 0.01)

  # ---------- cal_Bhat_Shat lowc_wc masking ----------
  tt1 <- cal_Bhat_Shat(update_Y, X, v1, lowc_wc = NULL)
  tt2 <- cal_Bhat_Shat(update_Y, X, v1, lowc_wc = 1:10)
  expect_equal(c(tt1$Bhat[, -c(1:10)]), c(tt2$Bhat[, -c(1:10)]))
  expect_equal(c(tt1$Shat[, -c(1:10)]), c(tt2$Shat[, -c(1:10)]))
  expect_equal(c(tt2$Bhat[, 1:10]), rep(0, length(c(tt2$Bhat[, 1:10]))))
  expect_equal(c(tt2$Shat[, 1:10]), rep(1, length(c(tt2$Shat[, 1:10]))))
})

test_that("susiF single-effect recovers one curve", {
  set.seed(1)
  sim <- simu_test_function(rsnr = 2, is.plot = FALSE)
  Y <- sim$noisy.data
  X <- sim$G
  out <- susiF(Y, X, L = 1, prior = "mixture_normal_per_scale", nullweight = 0)
  expect_equal(unlist(out$alpha), c(1, rep(0, 9)), tol = 1e-5)
  expect_equal(sum(abs(unlist(out$fitted_func) - sim$f1)), 0, tol = 0.2 * length(sim$f1))
})

test_that("susiF two-effects recovers two curves", {

  set.seed(1)
  sim <- simu_test_function(rsnr = 1, pos2 = 2, is.plot = FALSE)
  Y <- sim$noisy.data
  X <- sim$G
  out <- susiF(Y, X, L = 2, prior = "mixture_normal_per_scale",
               nullweight = 0, init_pi0_w = 1)
  expect_equal(Reduce("+", out$alpha), c(1, 1, rep(0, 8)), tol = 1e-5)

  d1 <- min(sqrt(mean( unlist(out$fitted_func[[1]]) - sim$f1)^2),
            sqrt(mean( unlist(out$fitted_func[[1]]) - sim$f2)^2))
  d2 <- min(sqrt(mean( unlist(out$fitted_func[[2]]) - sim$f1)^2),
            sqrt(mean( unlist(out$fitted_func[[2]]) - sim$f2)^2))
  expect_lte(d1, 0.25  )
  expect_lte(d2, 0.25)

})

test_that("Objective, variance, and KL update hook up", {
  # Reuse fixtures from the first test by regenerating quickly
  set.seed(2)
  f1 <- simu_IBSS_per_level(lev_res = 9, alpha = 1, prop_decay = 1.5)
  N <- 200; P <- 6; rsnr <- 10; pos1 <- 3
  G <- matrix(sample(c(0,1,2), size = N * P, replace = TRUE), nrow = N)
  Y <- t(replicate(N, f1$sim_func)) + 0 * G[, pos1] # simpler for speed
  Y <- Y + matrix(rnorm(length(Y), sd = (1/rsnr) * sd(f1$sim_func)), nrow = N)

  Y <- colScale(Y, scale = FALSE)
  X <- colScale(G)

  W <- DWT2(Y); Y_f <- cbind(W$D, W$C)
  v1 <- rep(1, ncol(X))
  indx_lst <- gen_wavelet_indx(9)
  tt <- cal_Bhat_Shat(Y_f, X, v1, lowc_wc = NULL)
  Bhat <- tt$Bhat; Shat <- tt$Shat
  nullweight <- 10 / sqrt(N)

  G_prior <- init_prior(
    Y = Y_f, X = X,
    prior = "mixture_normal_per_scale",
    v1 = v1, indx_lst = indx_lst,
    lowc_wc = NULL,
    control_mixsqp = control_mixsqp,
    nullweight = nullweight,
    max_SNP_EM = 100,
    tol_null_prior = 0
  )$G_prior

  susiF_obj <- init_susiF_obj(L_max = 2, G_prior, Y, X, L_start = 2,
                              greedy = TRUE, backfit = TRUE)

  outEM <- EM_pi(G_prior, Bhat, Shat, indx_lst,
                 init_pi0_w = 1, control_mixsqp = control_mixsqp,
                 lowc_wc = NULL, nullweight = nullweight, tol_null_prior = 0)
  G_prior <- update_prior(G_prior, tpi = outEM$tpi_k)
  susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst)

  sigma2 <- estimate_residual_variance(susiF_obj, Y_f, X)
  susiF_obj <- update_residual_variance(susiF_obj, sigma2 = sigma2)
  expect_equal(susiF_obj$sigma2, sigma2)

  KL_l <- cal_KL_l(susiF_obj, l = 1, Y = Y_f, X = X, D = W$D, C = W$C, indx_lst = indx_lst)
  susiF_obj <- update_KL(susiF_obj, Y = Y_f, X = X, D = W$D, C = W$C, indx_lst = indx_lst)
  expect_equal(susiF_obj$KL[1], KL_l)

  # Just ensure objective computes without error
  expect_silent(get_objective(susiF_obj, Y_f, X, D = W$D, C = W$C, indx_lst = indx_lst))
})

