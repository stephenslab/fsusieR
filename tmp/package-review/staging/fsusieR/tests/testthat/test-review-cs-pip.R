review_prior <- function(per_scale = FALSE, null = FALSE) {
  g <- list(fitted_g = ashr::normalmix(if (null) c(1, 0) else c(.5, .5),
                                      0, c(0, 1)))
  structure(rep(list(g), if (per_scale) 4L else 1L),
            class = if (per_scale) "mixture_normal_per_scale" else "mixture_normal")
}

review_data <- function() {
  set.seed(901)
  X <- matrix(rnorm(40 * 5), 40, 5)
  Y <- outer(X[, 1], c(0, 0, 1, 2, 3, 2, 1, 0)) +
    matrix(rnorm(40 * 8), 40, 8)
  list(X = X, Y = Y)
}

review_fit <- function(...) {
  args <- modifyList(c(review_data(), list(
    L = 2L, L_start = 2L, maxit = 5L, tol = 0,
    greedy = FALSE, backfit = FALSE, filter_cs = FALSE,
    prior = "mixture_normal_per_scale", post_processing = "none",
    verbose = FALSE, nullweight = 1, e = 0
  )), list(...))
  suppressWarnings(do.call(susiF, args))
}

test_that("overlapping CSs retain every plotted membership", {
  obj <- structure(list(L = 2L, P = 3L,
                        alpha = list(c(.94, .06, 0), c(.94, 0, .06)),
                        cs = vector("list", 2)), class = "susiF")
  obj <- update_cal_pip(update_cal_cs(obj, .95))
  expect_equal(obj$cs, list(c(1L, 2L), c(1L, 3L)))
  p <- plot_susiF_pip(obj)
  memberships <- p$layers[[2]]$data
  expect_equal(memberships$snp, c(1L, 2L, 1L, 3L))
  expect_equal(as.integer(memberships$CS), c(1L, 1L, 2L, 2L))
  expect_equal(memberships$y, obj$pip[memberships$snp])
  expect_equal(length(unique(memberships$ring_size[memberships$snp == 1])), 2)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[2]]), 4)
  # Requesting only the PIP plot must not require fitted effect curves.
  expect_s3_class(plot_susiF(obj, which_plot = "pip"), "ggplot")
  expect_error(plot_susiF_pip(obj, pos_SNP = 1:2), "original SNP")
})

test_that("null effects contribute neither PIPs nor credible sets", {
  for (per_scale in c(FALSE, TRUE)) {
    g <- review_prior(per_scale, null = TRUE)
    obj <- init_susiF_obj(2, g, matrix(0, 5, 8), matrix(1, 5, 3),
                          L_start = 2, greedy = FALSE, backfit = FALSE)
    obj <- update_cal_pip(update_cal_cs(obj))
    expect_equal(obj$pip, rep(0, 3))
    expect_equal(lengths(obj$cs), c(0L, 0L))
    expect_equal(obj$n_cs, 0)
    obj <- check_cs(obj, X = matrix(rnorm(15), 5, 3))
    expect_equal(obj$L, 1)
    expect_equal(obj$n_cs, 0)
    expect_equal(obj$cs[[1]], integer(0))
    for (post in c("none", "smash", "TI", "HMM")) {
      out <- update_cal_fit_func(obj, indx_lst = gen_wavelet_indx(3),
                                 Y = matrix(rnorm(40), 5, 8),
                                 X = matrix(rnorm(15), 5, 3),
                                 post_processing = post)
      expect_equal(out$fitted_func[[1]], rep(0, 8))
    }
  }
})

test_that("a nonnull scale is not hidden by averaging null weights", {
  obj <- list(alpha = list(rep(1/3, 3)), G_prior = review_prior(TRUE),
              est_pi = list(c(rep(list(c(1, 0)), 3), list(c(.998, .002)))),
              tol_null_prior = .001)
  expect_false(fsusieR:::.fs_null_effects(obj))
})

test_that("CS coverage uses the requested threshold and excludes zero mass", {
  obj <- structure(list(L = 1L, P = 3L, alpha = list(c(.8, .2, 0)),
                        cs = list(integer())), class = "susiF")
  expect_equal(update_cal_cs(obj, .8)$cs[[1]], 1L)
  expect_equal(update_cal_cs(obj, 1)$cs[[1]], c(1L, 2L))
  fit <- review_fit(cov_lev = .8)
  for (l in seq_along(fit$cs)) {
    if (!length(fit$cs[[l]])) next
    mass <- sort(fit$alpha[[l]], decreasing = TRUE)
    expected <- which(cumsum(mass) >= .8 - 1e-12)[1L]
    expect_length(fit$cs[[l]], expected)
  }
})

test_that("short fits perform complete sweeps without unfitted greedy effects", {
  for (iterations in 1:2) {
    fit <- review_fit(maxit = iterations, L = 4, L_start = 1,
                      greedy = TRUE, backfit = FALSE)
    expect_equal(fit$niter, iterations)
    expect_true(all(is.finite(unlist(fit$lBF))))
    expect_equal(vapply(fit$alpha, sum, numeric(1)), rep(1, fit$L))
    expect_length(fit$lbf, fit$L)
  }
  one <- review_fit(L = 1L, L_start = 1L, maxit = 4)
  expect_equal(one$niter, 4)
  expect_false(one$converged)
})

test_that("effect log BFs and all SNP fields retain original column mapping", {
  dat <- review_data()
  X <- cbind(constant = 1, dat$X[, 1:2], constant2 = 0, dat$X[, 3:5])
  fit <- review_fit(X = X)
  expect_equal(fit$variable_index, c(2L, 3L, 5L, 6L, 7L))
  expect_equal(unname(fit$pip[c(1, 4)]), c(0, 0))
  expect_equal(fit$lbf, vapply(fit$lBF, function(x)
    fsusieR:::.fs_log_mean_exp(x[fit$variable_index]), numeric(1)))
  for (l in seq_along(fit$cs)) {
    if (length(fit$cs[[l]])) {
      expect_gte(sum(fit$alpha[[l]][fit$cs[[l]]]), .95 - 1e-12)
      expect_true(all(fit$pip >= fit$alpha[[l]] - 1e-12))
    }
  }
  expect_equal(fit$n_cs, sum(lengths(fit$cs) > 0))
})

test_that("mixture BFs remain finite at extreme z scores and under the null", {
  b <- matrix(c(100, 80, 0, 1), 2, 2)
  for (df in list(NULL, 20)) {
    for (null in c(FALSE, TRUE)) {
      g <- review_prior(null = null)
      bf <- log_BF(g, b, matrix(1, 2, 2), NULL, list(1:2), df = df)
      expect_true(all(is.finite(bf)))
      if (null) expect_equal(bf, c(0, 0))
      per <- structure(rep(g, 2), class = "mixture_normal_per_scale")
      expect_equal(log_BF(per, b, matrix(1, 2, 2), NULL,
                          list(1L, 2L), df = df), bf)
    }
  }
  g <- review_prior()
  expected <- rowSums(matrix(log(.5 + .5 * exp(-log(2)/2 + c(b/100)^2/4)), 2))
  expect_equal(log_BF(g, b/100, matrix(1, 2, 2), NULL, list(1:2)), expected)
  L <- L_mixsq(g, b, matrix(1, 2, 2), list(1:2))
  expect_true(all(rowSums(unclass(L)) > 0))
})

test_that("a single SNP still learns its prior and respects nullweight", {
  b <- matrix(c(3, 0, -2, .2, 1, 4, -4, .1), 1)
  for (per_scale in c(FALSE, TRUE)) {
    g <- review_prior(per_scale)
    args <- list(G_prior = g, indx_lst = gen_wavelet_indx(3),
                 control_mixsqp = list(verbose = FALSE, eps = 1e-10),
                 lowc_wc = NULL, nullweight = 10, tol_null_prior = 0)
    single <- do.call(EM_pi, c(args, list(Bhat = b, Shat = matrix(1, 1, 8))))
    double <- do.call(EM_pi, c(args, list(Bhat = b[c(1, 1), ],
                                          Shat = matrix(1, 2, 8))))
    expect_equal(single$tpi_k, double$tpi_k, tolerance = 1e-7)
    expect_equal(single$lBF, double$lBF[1], tolerance = 1e-7)
  }
})

test_that("skipped conditional means and variances agree", {
  b <- matrix(1, 2, 8)
  for (per_scale in c(FALSE, TRUE)) {
    g <- review_prior(per_scale)
    args <- list(G_prior = g, Bhat = b, Shat = b, lBF = c(0, 0),
                 lowc_wc = NULL, indx_lst = gen_wavelet_indx(3), e = .9)
    expect_equal(do.call(post_mat_mean, args), 0 * b)
    expect_equal(do.call(post_mat_sd, args), 0 * b)
  }
})

test_that("raw reconstruction respects the wavelet filter", {
  dat <- review_data()
  fit <- review_fit(L = 1, L_start = 1, filter.number = 1,
                    family = "DaubExPhase")
  wc <- as.numeric((fit$alpha[[1]] / fit$csd_X) %*% fit$fitted_wc[[1]])
  template <- wavethresh::wd(rep(0, 8), filter.number = 1, family = "DaubExPhase")
  template$D <- wc[1:7]
  template$C[length(template$C)] <- wc[8]
  expect_equal(unname(fit$fitted_func[[1]]), unname(wavethresh::wr(template)))
  changed <- change_fit(fit, dat$Y, dat$X, to = "none", verbose = FALSE)
  expect_equal(changed$fitted_func, fit$fitted_func)
})

test_that("Gaussian input errors are caught before fitting", {
  dat <- review_data()
  expect_error(review_fit(X = dat$X[-1, ]), "matching rows")
  bad <- dat$X; bad[1, 1] <- NA_real_
  expect_error(review_fit(X = bad), "finite numeric")
  expect_error(review_fit(maxit = 0), "positive integer")
  expect_error(review_fit(cov_lev = 0), "cov_lev")
  expect_error(review_fit(pos = rep(1, 8)), "strictly increasing")
})

test_that("HMM refits handle overlapping CSs with duplicate lead SNPs", {
  calls <- 0L
  local_mocked_bindings(fit_hmm = function(x, sd, ...) {
    calls <<- calls + 1L
    list(posterior = list(mean = x, lfsr = rep(.1, length(x))))
  }, .package = "fsusieR")
  dat <- review_data()
  obj <- structure(list(L = 2L, cs = list(c(1L, 2L), c(1L, 3L)),
                        alpha = list(c(.94, .06, 0, 0, 0),
                                     c(.94, 0, .06, 0, 0))), class = "susiF")
  out <- HMM_regression(obj, dat$Y, dat$X)
  expect_equal(calls, 1L)
  expect_length(out$fitted_func, 2L)
  expect_equal(out$fitted_func[[2]], rep(0, 8))
  expect_equal(out$lfsr_func[[2]], rep(1, 8))
})

test_that("complete Gaussian coordinate sweeps increase the penalized ELBO", {
  for (prior in c("mixture_normal", "mixture_normal_per_scale")) {
    for (penalty in c(1, 10)) {
      fit <- review_fit(prior = prior, nullweight = penalty, cal_obj = TRUE,
                        maxit = 8, e = 0)
      expect_true(all(is.finite(fit$ELBO)))
      expect_true(all(diff(fit$ELBO) >= -1e-7), info = paste(prior, penalty))
    }
  }
})

test_that("interpolated response grids retain their original coordinates", {
  dat <- review_data()
  pos <- seq(100, 160, length.out = 7)
  eb <- EBmvFR(dat$Y[, 1:7], dat$X, pos = pos, maxit = 1, verbose = FALSE)
  expect_equal(range(eb$outing_grid), range(pos))
  for (post in c("none", "TI")) {
    fit <- review_fit(Y = dat$Y[, 1:7], pos = pos, post_processing = post)
    expect_length(fit$fitted_func[[1]], length(fit$outing_grid))
    expect_equal(ncol(fit$ind_fitted_func), length(fit$outing_grid))
    expect_equal(range(fit$outing_grid), range(pos))
  }
})

test_that("masked coefficients do not alter the learned mixture", {
  g <- review_prior()
  args <- list(G_prior = g, control_mixsqp = list(verbose = FALSE),
               nullweight = 1, tol_null_prior = 0, max_step = 5)
  single <- do.call(EM_pi, c(args, list(Bhat = matrix(3, 1, 1),
    Shat = matrix(1, 1, 1), indx_lst = list(1L), lowc_wc = NULL)))
  masked <- do.call(EM_pi, c(args, list(Bhat = matrix(c(3, rep(0, 7)), 1),
    Shat = matrix(1, 1, 8), indx_lst = gen_wavelet_indx(3), lowc_wc = 2:8)))
  expect_equal(single$tpi_k, masked$tpi_k, tolerance = 1e-6)
  expect_equal(single$lBF, masked$lBF, tolerance = 1e-6)
})

test_that("changing post-processing preserves original SNP indices and units", {
  dat <- review_data()
  X <- cbind(1, dat$X[, 1] * 3, 0, dat$X[, -1])
  fit <- review_fit(X = X)
  changed <- change_fit(fit, dat$Y, X, to = "smash", verbose = FALSE)
  direct <- review_fit(X = X, post_processing = "smash")
  expect_equal(changed$cs, fit$cs)
  expect_equal(changed$pip, fit$pip)
  expect_equal(unname(changed$fitted_func), unname(direct$fitted_func))
  expect_equal(changed$ind_fitted_func, direct$ind_fitted_func)
})

test_that("Poisson optimization gradients match the actual objective", {
  f <- fsusieR:::pois_mean_GP_opt_obj
  g <- fsusieR:::pois_mean_GP_opt_obj_gradient
  args <- list(x = c(2, 1e9), s = c(1, 2), beta = c(0, 3),
               sigma2 = c(.5, 2), n = 2)
  theta <- c(.3, 21, log(.2), -16)
  observed <- do.call(g, c(list(theta = theta), args))
  numerical <- vapply(seq_along(theta), function(k) {
    plus <- minus <- theta
    plus[k] <- plus[k] + 1e-4
    minus[k] <- minus[k] - 1e-4
    (do.call(f, c(list(theta = plus), args)) -
       do.call(f, c(list(theta = minus), args))) / 2e-4
  }, numeric(1))
  expect_equal(observed, numerical, tolerance = 1e-5)
  fit <- pois_mean_GP(c(0, 10, 1000), prior_mean = 1, prior_var = 2)
  expected <- vga_pois_solver(log(c(.1, 10, 1000)), c(0, 10, 1000),
                              rep(1, 3), 1, 2, tol = 1e-10)
  expect_equal(fit$posterior$posteriorMean_latent, expected$m, tolerance = 1e-4)
  expect_equal(fit$posterior$posteriorVar_latent, expected$v, tolerance = 1e-4)
})

test_that("Poisson smoothing returns the lognormal posterior variance", {
  fit <- pois_smooth_split(c(2, 3, 2, 4, 6, 5, 4, 3), maxiter = 2)
  expect_equal(fit$Vmean, expm1(fit$Vmu) * fit$Emean^2)
})

test_that("null simulation priors terminate and report the correct null fraction", {
  single <- simu_IBSS_ash_vanilla(lev_res = 3, pi0 = 1)
  scales <- simu_IBSS_per_level(lev_res = 3, pi0 = 1)
  expect_equal(single$sim_func, rep(0, 8))
  expect_equal(single$emp_pi0, 1)
  expect_equal(scales$sim_func, rep(0, 8))
  expect_equal(scales$emp_pi0, rep(1, 3))
})

test_that("simulation coefficients use their own sampled component scale", {
  local_mocked_bindings(rnorm = function(n, mean, sd) sd,
                        rchisq = function(n, df) rep(1, n),
                        runif = function(n) rep(1, n), .package = "fsusieR")
  set.seed(902)
  simulated <- simu_IBSS_ash_vanilla(lev_res = 4, length_grid = 3, pi0 = 0)
  expect_equal(sort(unique(simulated$true_coef)), c(1, 2))
})
