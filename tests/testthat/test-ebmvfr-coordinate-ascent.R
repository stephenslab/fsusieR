ebmvfr_test_prior <- function(per_scale = FALSE, T = 8L) {
  g <- list(fitted_g = ashr::normalmix(c(0.5, 0.3, 0.2), rep(0, 3), c(0, 0.4, 1.5)))
  if (per_scale) {
    structure(rep(list(g), log2(T) + 1L), class = "mixture_normal_per_scale")
  } else {
    structure(list(g), class = "mixture_normal")
  }
}

ebmvfr_test_data <- function(p = 3L, orthogonal = FALSE) {
  set.seed(953)
  X <- scale(matrix(rnorm(24L * p), 24L, p))
  if (orthogonal) X <- qr.Q(qr(X)) %*% diag(seq_len(p), p)
  if (!orthogonal && p > 1L) X[, 2] <- 0.9 * X[, 1] + 0.3 * X[, 2]
  Z <- scale(matrix(rnorm(24L * 8L), 24L, 8L), scale = FALSE)
  Z <- Z + tcrossprod(X[, 1], seq(-1, 1, length.out = 8))
  list(X = X, Z = Z, idx = gen_wavelet_indx(3))
}

test_that("EBmvFR normal means agree with ashr, including degenerate components", {
  bhat <- c(0, 0.01, -0.2, 1, -5, 20)
  shat <- c(1, 0.01, 0.3, 0.5, 1, 0.1)
  for (weights in list(c(0.5, 0.3, 0.2), c(0, 0.5, 0.5), c(1, 0, 0), c(0, 0, 1))) {
    g <- ashr::normalmix(weights, rep(0, 3), c(0, 0.4, 1.5))
    data <- ashr::set_data(bhat, shat)
    post <- .ebmvfr_normal_means(bhat, shat, g)
    expect_equal(post$mean, ashr::postmean(g, data), tolerance = 1e-10)
    expect_equal(post$variance, ashr::postsd(g, data)^2, tolerance = 1e-9)
    expect_equal(post$counts, rowSums(ashr::comp_postprob(g, data)), tolerance = 1e-10)
    expect_equal(sum(post$counts), length(bhat), tolerance = 1e-12)
    expect_true(is.finite(post$KL) && post$KL >= -1e-10)

    # Independent KL identity: posterior expected log likelihood minus
    # normal-means log evidence. Includes the categorical mixture entropy.
    likelihood <- vapply(seq_along(g$pi), function(k) {
      dnorm(bhat, 0, sqrt(shat^2 + g$sd[k]^2))
    }, numeric(length(bhat)))
    log_evidence <- log(drop(likelihood %*% g$pi))
    usable <- is.finite(log_evidence)
    expected_loglik <- -0.5 * (log(2 * pi * shat^2) +
      ((bhat - post$mean)^2 + post$variance) / shat^2)
    if (all(usable)) {
      expect_equal(post$KL, sum(expected_loglik - log_evidence), tolerance = 1e-9)
    }
  }
  one <- .ebmvfr_normal_means(0.5, 0.2, ashr::normalmix(1, 0, 1))
  expect_equal(one$mean, 0.5 / 1.04)
  expect_equal(one$variance, 0.04 / 1.04)
})

test_that("EBmvFR uses current sigma2 and retains standard errors without another root", {
  dat <- ebmvfr_test_data()
  obj <- init_EBmvFR_obj(ebmvfr_test_prior(), dat$Z, dat$X)
  obj$sigma2 <- 0.25
  low <- fit_effect.EBmvFR(obj, 1, dat$X, dat$Z[, -8], dat$Z[, 8], dat$idx, NULL)
  obj$sigma2 <- 25
  high <- fit_effect.EBmvFR(obj, 1, dat$X, dat$Z[, -8], dat$Z[, 8], dat$idx, NULL)
  expected_se <- sqrt(0.25 / sum(dat$X[, 1]^2))
  expect_equal(low$MLE_wc2[[1]][1, ], rep(expected_se, 8))
  expect_equal(high$MLE_wc2[[1]][1, ], rep(10 * expected_se, 8))
  expect_gt(max(abs(low$fitted_wc[[1]][1, ] - high$fitted_wc[[1]][1, ])), 0.01)
})

test_that("EBmvFR partial residuals work for one and two predictors and match the cache", {
  for (p in 1:3) {
    dat <- ebmvfr_test_data(p)
    obj <- init_EBmvFR_obj(ebmvfr_test_prior(), dat$Z, dat$X)
    obj$fitted_wc[[1]] <- matrix(seq_len(p * 8) / 20, p, 8)
    for (j in seq_len(p)) {
      direct <- cal_partial_resid.EBmvFR(obj, j, dat$X, dat$Z[, -8], dat$Z[, 8], dat$idx)
      expected <- dat$Z
      for (k in setdiff(seq_len(p), j)) {
        expected <- expected - tcrossprod(dat$X[, k], obj$fitted_wc[[1]][k, ])
      }
      expect_equal(dim(direct), dim(expected))
      expect_equal(as.vector(direct), as.vector(expected), tolerance = 1e-12)
      obj$residual <- dat$Z - dat$X %*% obj$fitted_wc[[1]]
      cached <- cal_partial_resid.EBmvFR(obj, j, dat$X, dat$Z[, -8], dat$Z[, 8], dat$idx)
      expect_equal(as.vector(cached), as.vector(direct), tolerance = 1e-12)
      obj <- fit_effect.EBmvFR(obj, j, dat$X, dat$Z[, -8], dat$Z[, 8], dat$idx, NULL)
      expect_equal(obj$residual, dat$Z - dat$X %*% obj$fitted_wc[[1]], tolerance = 1e-12)
      obj$residual <- NULL
    }
  }
})

test_that("EBmvFR uses a Frobenius residual norm and predictor-weighted posterior variances", {
  X <- matrix(c(-1, 0, 1), 3, 1)
  Z <- cbind(c(-1, 0, 1), c(1, 0, -1))
  obj <- structure(list(fitted_wc = list(matrix(0, 1, 2)),
                        fitted_wc2 = list(matrix(0, 1, 2))), class = "EBmvFR")
  expect_equal(get_ER2(obj, Z, X), 4)
  obj$fitted_wc[[1]][] <- 1
  obj$fitted_wc2[[1]][] <- 0.01
  Z <- X %*% obj$fitted_wc[[1]]
  expect_equal(get_ER2(obj, Z, X), 0.04)
  expect_equal(estimate_residual_variance(obj, Z, X), 0.04 / 6)

  dat <- ebmvfr_test_data()
  obj <- init_EBmvFR_obj(ebmvfr_test_prior(), dat$Z, dat$X)
  obj$fitted_wc[[1]] <- matrix(seq_len(24) / 30, 3, 8)
  obj$fitted_wc2[[1]] <- matrix(seq_len(24) / 100, 3, 8)
  expected <- sum((dat$Z - dat$X %*% obj$fitted_wc[[1]])^2)
  for (k in seq_len(8)) {
    cov_fit <- dat$X %*% diag(obj$fitted_wc2[[1]][, k]) %*% t(dat$X)
    expected <- expected + sum(diag(cov_fit))
  }
  expect_equal(get_ER2(obj, dat$Z, dat$X), expected)
})

test_that("EBmvFR prior update is Kim's fixed-q responsibility average", {
  dat <- ebmvfr_test_data()
  prior <- structure(list(list(fitted_g = ashr::normalmix(c(0.5, 0.5), c(0, 0), c(0, 1)))),
                     class = "mixture_normal")
  obj <- init_EBmvFR_obj(prior, dat$Z, dat$X)
  mle <- list(Bhat = matrix(0, 1, 8), Shat = matrix(1, 1, 8))
  for (j in 1:3) obj <- update_effect.EBmvFR(obj, j, mle, NULL, dat$idx)
  before <- .ebmvfr_objective(obj, dat$Z, dat$X)
  result <- update_prior_EBmvFR(obj)
  expected <- c(sqrt(2), 1) / (sqrt(2) + 1)
  expect_equal(as.numeric(get_pi_G_prior(result$G_prior)), expected)
  expect_equal(result$fitted_wc, obj$fitted_wc)
  expect_equal(result$fitted_wc2, obj$fitted_wc2)
  expect_gt(.ebmvfr_objective(result, dat$Z, dat$X), before)
  expect_equal(as.numeric(result$est_pi[[1]]), expected)

  # Repeating an M-step with the SAME q must leave pi unchanged.
  repeated <- update_prior_EBmvFR(result)
  expect_equal(repeated$G_prior, result$G_prior)
  expect_equal(repeated$KL, result$KL)
})

test_that("EBmvFR ELBO equals exact evidence for orthogonal X and fixed hyperparameters", {
  dat <- ebmvfr_test_data(orthogonal = TRUE)
  for (per_scale in c(FALSE, TRUE)) {
    prior <- ebmvfr_test_prior(per_scale)
    obj <- init_EBmvFR_obj(prior, dat$Z, dat$X)
    obj$sigma2 <- 0.8
    for (j in seq_len(ncol(dat$X))) {
      obj <- fit_effect.EBmvFR(obj, j, dat$X, dat$Z[, -8], dat$Z[, 8], dat$idx, NULL)
    }
    exact <- sum(dnorm(dat$Z, 0, sqrt(obj$sigma2), log = TRUE))
    for (j in seq_len(ncol(dat$X))) {
      se <- sqrt(obj$sigma2 / sum(dat$X[, j]^2))
      bhat <- drop(crossprod(dat$X[, j], dat$Z)) / sum(dat$X[, j]^2)
      for (s in seq_along(prior)) {
        idx <- obj$prior_groups[[s]]
        g <- prior[[s]]$fitted_g
        density <- vapply(seq_along(g$pi), function(k) {
          dnorm(bhat[idx], 0, sqrt(se^2 + g$sd[k]^2))
        }, numeric(length(idx)))
        if (length(idx) == 1) density <- matrix(density, 1)
        exact <- exact + sum(log(drop(density %*% g$pi)) - dnorm(bhat[idx], 0, se, log = TRUE))
      }
    }
    expect_equal(.ebmvfr_objective(obj, dat$Z, dat$X), exact, tolerance = 1e-10)
  }
})

test_that("every EBmvFR coordinate block increases the same ELBO for correlated X", {
  dat <- ebmvfr_test_data()
  for (per_scale in c(FALSE, TRUE)) {
    obj <- init_EBmvFR_obj(ebmvfr_test_prior(per_scale), dat$Z, dat$X)
    previous <- .ebmvfr_objective(obj, dat$Z, dat$X)
    for (iter in 1:8) {
      for (j in 1:3) {
        obj <- fit_effect.EBmvFR(obj, j, dat$X, dat$Z[, -8], dat$Z[, 8], dat$idx, NULL)
        value <- .ebmvfr_objective(obj, dat$Z, dat$X)
        expect_gte(value - previous, -1e-9)
        previous <- value
      }
      obj <- update_prior_EBmvFR(obj)
      value <- .ebmvfr_objective(obj, dat$Z, dat$X)
      expect_gte(value - previous, -1e-9)
      previous <- value
      obj <- update_residual_variance(obj, estimate_residual_variance(obj, dat$Z, dat$X))
      value <- .ebmvfr_objective(obj, dat$Z, dat$X)
      expect_gte(value - previous, -1e-9)
      previous <- value
    }
  }
})

test_that("masked wavelet columns do not bias EBmvFR priors or residual variance", {
  dat <- ebmvfr_test_data()
  mask <- c(dat$idx[[1]], dat$idx[[2]])
  prior <- ebmvfr_test_prior(TRUE)
  obj <- init_EBmvFR_obj(prior, dat$Z, dat$X, lowc_wc = mask)
  for (j in 1:3) obj <- fit_effect.EBmvFR(obj, j, dat$X, dat$Z[, -8], dat$Z[, 8], dat$idx, mask)
  obj <- update_prior_EBmvFR(obj)
  expect_true(all(obj$fitted_wc[[1]][, mask] == 0))
  expect_true(all(obj$fitted_wc2[[1]][, mask] == 0))
  expect_equal(obj$G_prior[1:2], prior[1:2])
  changed <- dat$Z
  changed[, mask] <- 1e6
  expect_equal(get_ER2(obj, changed, dat$X), get_ER2(obj, dat$Z, dat$X))
  expect_equal(estimate_residual_variance(obj, changed, dat$X),
               get_ER2(obj, dat$Z, dat$X) / (nrow(dat$Z) * (8 - length(mask))))
  expect_true(is.finite(.ebmvfr_objective(obj, dat$Z, dat$X)))
})

test_that("EBmvFR public fits support small p, ELBO and parameter stopping", {
  for (p in 1:2) {
    dat <- ebmvfr_test_data(p)
    for (prior in c("mixture_normal", "mixture_normal_per_scale")) {
      fit <- EBmvFR(dat$Z, dat$X, prior = prior, cal_obj = TRUE,
                    verbose = FALSE, maxit = 12, tol = 0)
      expect_equal(fit$niter, 12)
      expect_false(fit$converged)
      expect_length(fit$ELBO, fit$niter + 1)
      expect_true(all(is.finite(fit$ELBO)))
      expect_true(all(diff(fit$ELBO) >= -1e-8))
      expect_gt(fit$sigma2, 0)
      expect_null(fit$residual)
      expect_equal(dim(fit$fitted_func), c(p, 8L))
      expect_equal(fit$ind_fitted_func,
                   unname(scale(dat$X, scale = FALSE)) %*% fit$fitted_func,
                   ignore_attr = TRUE)
    }
  }
  dat <- ebmvfr_test_data()
  loose <- EBmvFR(dat$Z, dat$X, verbose = FALSE, cal_obj = FALSE, tol = 1e6, maxit = 8)
  strict <- EBmvFR(dat$Z, dat$X, verbose = FALSE, cal_obj = FALSE, tol = 0, maxit = 2)
  one <- EBmvFR(dat$Z, dat$X, verbose = FALSE, maxit = 1)
  expect_equal(loose$niter, 1)
  expect_true(loose$converged)
  expect_equal(strict$niter, 2)
  expect_false(strict$converged)
  expect_equal(one$niter, 1)
  expect_length(loose$ELBO, 0)
  expect_error(EBmvFR(dat$Z, cbind(1, dat$X), verbose = FALSE), "constant")
})
