## ----------------------------------------------------------------------------
## Tests for the mixture_normal prior path in fsusieR.
##
## Layout:
##   1. Shared simulation + IBSS scaffolding (run once at top, no test_that)
##   2. Tests grouped by topic
## ----------------------------------------------------------------------------

library(testthat)
library(ashr)
library(wavethresh)
library(mixsqp)
library(fsusieR)

# ---- 1. shared simulation -------------------------------------------------

set.seed(2)
f1   <- simu_IBSS_per_level(lev_res = 9, alpha = 1, prop_decay = 1.5)

N    <- 500
P    <- 10
pos1 <- 5
beta1      <- 1
rsnr       <- 1
nullweight <- 0.1
greedy     <- TRUE
backfit    <- TRUE
control_mixsqp <- list(verbose = FALSE)
init_pi0_w     <- 1
lowc_wc        <- NULL

set.seed(23)
G <- matrix(sample(c(0, 1, 2), size = N * P, replace = TRUE), nrow = N, ncol = P)

noisy.data <- matrix(NA, N, length(f1$sim_func))
for (i in seq_len(N)) {
  noisy.data[i, ] <- beta1 * G[i, pos1] * f1$sim_func +
                     rnorm(length(f1$sim_func),
                           sd = (1 / rsnr) * sd(f1$sim_func))
}

Y <- colScale(noisy.data, scale = FALSE)
X <- colScale(G)

W        <- DWT2(Y)
Y_f      <- cbind(W$D, W$C)
indx_lst <- gen_wavelet_indx(9)
v1       <- rep(1, dim(X)[2])

## cal_Bhat_Shat in marginal mode (sigma2 omitted) is the natural fit for
## these single-modality tests.
tt   <- cal_Bhat_Shat(Y_f, X, lowc_wc = NULL)
Bhat <- tt$Bhat
Shat <- tt$Shat

G_prior <- init_prior(Y              = Y_f,
                      X              = X,
                      prior          = "mixture_normal",
                      v1             = v1,
                      indx_lst       = indx_lst,
                      lowc_wc        = NULL,
                      control_mixsqp = control_mixsqp,
                      nullweight     = nullweight,
                      max_SNP_EM     = 100,
                      tol_null_prior = 0)$G_prior

lBF <- log_BF(G_prior, Bhat, Shat, indx_lst = indx_lst, lowc_wc = NULL)


# ---- 2a. prior class ------------------------------------------------------

test_that("init_prior returns a mixture_normal G_prior", {
  expect_s3_class(
    init_prior(Y              = Y_f,
               X              = X,
               prior          = "mixture_normal",
               v1             = v1,
               indx_lst       = indx_lst,
               lowc_wc        = NULL,
               control_mixsqp = control_mixsqp,
               nullweight     = nullweight,
               max_SNP_EM     = 100,
               tol_null_prior = 0)$G_prior,
    "mixture_normal"
  )
})

test_that("get_pi_G_prior / get_sd_G_prior have expected classes", {
  expect_s3_class(get_pi_G_prior(G_prior), "pi_mixture_normal")
  expect_s3_class(get_sd_G_prior(G_prior), "sd_mixture_normal")
})


# ---- 2b. lBF locates the causal SNP --------------------------------------

test_that("Max lBF is at the causal SNP position", {
  expect_equal(which.max(lBF), pos1)
})


# ---- 2c. EM step machinery -----------------------------------------------

test_that("L_mixsq returns lik_mixture_normal", {
  L <- L_mixsq(G_prior, Bhat, Shat, indx_lst)
  expect_s3_class(L, "lik_mixture_normal")
})

test_that("cal_zeta puts the highest assignment at the causal SNP", {
  zeta <- cal_zeta(lBF)
  expect_equal(which.max(zeta), pos1)
})

test_that("m_step output has class pi_mixture_normal", {
  L    <- L_mixsq(G_prior, Bhat, Shat, indx_lst)
  zeta <- cal_zeta(lBF)
  tpi  <- m_step(L, zeta, indx_lst,
                 init_pi0_w     = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  expect_s3_class(tpi, "pi_mixture_normal")
})

test_that("m_step with init_pi0_w = 1 returns a near-null mixture", {
  L    <- L_mixsq(G_prior, Bhat, Shat, indx_lst)
  zeta <- cal_zeta(lBF)
  tpi  <- m_step(L, zeta, indx_lst,
                 init_pi0_w     = 1,
                 control_mixsqp = control_mixsqp,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  expect_gt(fsusieR::get_pi0(tpi = tpi), (1 - 1 / 2^9) - 0.01)
})

test_that("update_prior writes tpi back into G_prior unchanged", {
  L    <- L_mixsq(G_prior, Bhat, Shat, indx_lst)
  zeta <- cal_zeta(lBF)
  tpi  <- m_step(L, zeta, indx_lst,
                 init_pi0_w     = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  G_update <- update_prior(G_prior, tpi)
  expect_identical(get_pi_G_prior(G_update), tpi)
})

test_that("EM_pi outputs have expected classes / shapes", {
  outEM <- EM_pi(G_prior, Bhat, Shat, indx_lst,
                 init_pi0_w     = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 lowc_wc        = NULL,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  expect_s3_class(outEM$tpi_k, "pi_mixture_normal")
  expect_type(outEM$lBF, "double")
  expect_length(outEM$lBF, dim(Bhat)[1])
  expect_gt(fsusieR::get_pi0(tpi = outEM$tpi_k), (1 - 1 / 2^9) - 0.01)
})


# ---- 2d. susiF object construction ---------------------------------------

test_that("init_susiF_obj seeds pi from the supplied G_prior (L = 1)", {
  obj <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                        L_start = 1, greedy = greedy, backfit = backfit)
  expect_equal(get_pi(obj, 1), get_pi_G_prior(G_prior))
})

test_that("init_susiF_obj seeds pi from the supplied G_prior (L = 2)", {
  obj <- init_susiF_obj(L_max = 2, G_prior, Y, X,
                        L_start = 2, greedy = greedy, backfit = backfit)
  expect_equal(get_pi(obj, 1), get_pi_G_prior(G_prior))
  expect_equal(get_pi(obj, 2), get_pi_G_prior(G_prior))
})

test_that("get_G_prior round-trips through init_susiF_obj", {
  obj <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                        L_start = 1, greedy = greedy, backfit = backfit)
  expect_equal(get_G_prior(obj), G_prior)
})


# ---- 2e. update_alpha / update_pi / update_susiF_obj ----------------------

test_that("update_alpha.susiF writes alpha back into the susiF object", {
  obj   <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                          L_start = 1, greedy = greedy, backfit = backfit)
  outEM <- EM_pi(G_prior, Bhat, Shat, indx_lst,
                 init_pi0_w     = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 lowc_wc        = NULL,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  new_alpha <- cal_zeta(outEM$lBF)
  obj <- update_alpha.susiF(obj, 1, new_alpha)
  expect_equal(get_alpha(obj, 1), new_alpha)
})

test_that("update_pi writes tpi back into the susiF object", {
  obj   <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                          L_start = 1, greedy = greedy, backfit = backfit)
  outEM <- EM_pi(G_prior, Bhat, Shat, indx_lst,
                 init_pi0_w     = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 lowc_wc        = NULL,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  obj <- update_pi(obj, 1, outEM$tpi_k)
  expect_equal(get_pi(obj, 1), outEM$tpi_k)
})

test_that("update_susiF_obj sets fitted_wc / fitted_wc2 / alpha / G_prior", {
  obj   <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                          L_start = 1, greedy = greedy, backfit = backfit)
  outEM <- EM_pi(G_prior, Bhat, Shat, indx_lst,
                 init_pi0_w     = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 lowc_wc        = NULL,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  G_post <- update_prior(G_prior, tpi = outEM$tpi_k)
  obj    <- update_susiF_obj(obj, 1, outEM, Bhat, Shat, indx_lst)

  expect_equal(obj$fitted_wc[[1]],
               post_mat_mean(G_post, Bhat, Shat,
                             lBF      = outEM$lBF,
                             indx_lst = indx_lst,
                             lowc_wc  = lowc_wc))
  expect_equal(obj$fitted_wc2[[1]],
               post_mat_sd(G_post, Bhat, Shat,
                           lBF      = outEM$lBF,
                           indx_lst = indx_lst,
                           lowc_wc  = lowc_wc)^2)
  expect_equal(get_alpha(obj, 1), cal_zeta(outEM$lBF))
  expect_equal(get_G_prior(obj),  G_post)
})


# ---- 2f. partial residual ------------------------------------------------

test_that("cal_partial_resid matches the manual residualisation", {
  obj   <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                          L_start = 1, greedy = greedy, backfit = backfit)
  outEM <- EM_pi(G_prior, Bhat, Shat, indx_lst,
                 init_pi0_w     = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 lowc_wc        = NULL,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  obj <- update_susiF_obj(obj, 1, outEM, Bhat, Shat, indx_lst)

  update_T <- cal_partial_resid(
    obj      = obj,
    l        = 1,
    X        = X,
    D        = W$D,
    C        = W$C,
    indx_lst = indx_lst
  )

  id_L <- 1
  update_D <- W$D - Reduce("+", lapply(id_L, function(l)
    (X * rep(obj$alpha[[l]], rep.int(N, P))) %*%
      obj$fitted_wc[[l]][, -dim(obj$fitted_wc[[l]])[2]]))
  update_C <- W$C - Reduce("+", lapply(id_L, function(l)
    (X * rep(obj$alpha[[l]], rep.int(N, P))) %*%
      obj$fitted_wc[[l]][,  dim(obj$fitted_wc[[l]])[2]]))
  manual_update <- cbind(update_D, update_C)

  expect_equal(update_T, manual_update)
})


# ---- 2g. fitted curve recovery -------------------------------------------

test_that("update_cal_fit_func (post_processing='none') recovers f1", {
  obj   <- init_susiF_obj(L_max = 1, G_prior, Y, X,
                          L_start = 1, greedy = greedy, backfit = backfit)
  outEM <- EM_pi(G_prior, Bhat, Shat, indx_lst,
                 init_pi0_w     = init_pi0_w,
                 control_mixsqp = control_mixsqp,
                 lowc_wc        = NULL,
                 nullweight     = nullweight,
                 tol_null_prior = 0)
  obj <- update_susiF_obj(obj, 1, outEM, Bhat, Shat, indx_lst)
  obj <- update_cal_fit_func(obj      = obj,
                             indx_lst = indx_lst,
                             Y        = Y,
                             X        = X,
                             post_processing = "none")
  ## Mean absolute error per position should be small relative to sd(f1).
  mae <- mean(abs(obj$fitted_func[[1]] - f1$sim_func))
  expect_lt(mae, 0.1 * sd(f1$sim_func))
})


# ---- 2h. end-to-end susiF performance ------------------------------------

test_that("susiF (L = 1) recovers a single causal effect with mixture_normal", {
  set.seed(1)
  sim <- simu_test_function(rsnr = 0.5, is.plot = FALSE)
  out <- susiF(sim$noisy.data, sim$G, L = 1, prior = "mixture_normal")

  expect_equal(unlist(out$alpha), c(1, rep(0, 9)), tolerance = 1e-5)
  expect_lt(mean((unlist(out$fitted_func) - sim$f1)^2),
            0.2 * length(sim$f1))
})

test_that("susiF (L = 2) recovers two causal effects with mixture_normal", {
  set.seed(1)
  sim <- simu_test_function(N = 500, rsnr = 1.5, pos2 = 2, is.plot = FALSE)
  out <- susiF(sim$noisy.data, sim$G, L = 2, prior = "mixture_normal")



  expect_equal(Reduce("+", out$alpha), c(1, 1, rep(0, 8)), tolerance = 1e-5)

  ## fitted_func[[1]] and fitted_func[[2]] each correspond to one of the two
  ## true effects (we don't know the ordering); so for each fitted curve, the
  ## min absolute error against {f1, f2} should be small.
  err1<- min(sqrt(mean((out$fitted_func[[1]] - sim$f1)^2)),
             sqrt(mean((out$fitted_func[[1]] - sim$f2)^2)))
  err2 <- min(sqrt(mean((out$fitted_func[[2]] - sim$f1)^2)),
              sqrt(mean((out$fitted_func[[2]] - sim$f2)^2)))
  expect_lt(err1, 0.2 * length(sim$f1))
  expect_lt(err2, 0.2 * length(sim$f1))

})

