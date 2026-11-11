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
