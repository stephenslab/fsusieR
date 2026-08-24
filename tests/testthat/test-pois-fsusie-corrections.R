test_that("Poisson inputs are validated", {
  set.seed(1)
  Y <- matrix(rpois(24, 3), nrow = 4)
  X <- cbind(c(0, 1, 2, 1), c(2, 0, 1, 2))

  checked <- fsusieR:::validate_pois_fsusie_inputs(Y, X, 2, rep(1, 4))
  expect_identical(dim(checked$Y), c(4L, 6L))
  expect_error(
    fsusieR:::validate_pois_fsusie_inputs(replace(Y, 1, -1), X, 2, rep(1, 4)),
    "non-negative integer"
  )
  expect_error(
    fsusieR:::validate_pois_fsusie_inputs(replace(Y, 1, 0.5), X, 2, rep(1, 4)),
    "non-negative integer"
  )
  expect_error(
    fsusieR:::validate_pois_fsusie_inputs(Y, X, 2, c(1, 1, 1, 0)),
    "strictly positive"
  )
  expect_error(
    fsusieR:::validate_pois_fsusie_inputs(Y, cbind(X, 1), 2, rep(1, 4)),
    "constant column"
  )
})


test_that("B posterior variance is centered and includes curve variance", {
  mock <- list(
    cs = list(c(1L, 2L)),
    alpha = list(c(0.25, 0.75)),
    fitted_func = list(c(2, 3)),
    fitted_var = list(c(0.5, 0.25))
  )
  X <- cbind(c(0, 1, 2), c(2, 0, 1))
  X_centered <- sweep(X, 2, colMeans(X), "-")
  alpha <- mock$alpha[[1]]
  mean_curve <- mock$fitted_func[[1]]
  curve_variance <- mock$fitted_var[[1]]
  expected <- tcrossprod(
    as.numeric((X_centered^2) %*% alpha),
    mean_curve^2 + curve_variance
  ) - tcrossprod(
    as.numeric(X_centered %*% alpha)^2,
    mean_curve^2
  )

  observed <- fsusieR:::compute_B_pv(mock, X, n_grid = 2)
  expect_equal(observed, pmax(expected, 0), tolerance = 1e-12)
  expect_equal(
    observed,
    fsusieR:::compute_B_pv(
      mock, sweep(X, 2, c(10, -7), "+"), n_grid = 2
    ),
    tolerance = 1e-12
  )
})


test_that("no-effect B variance has the requested grid size", {
  no_effect <- list(cs = list(), alpha = list(), fitted_func = list())
  X <- cbind(c(0, 1, 2), c(2, 0, 1))
  observed <- fsusieR:::compute_B_pv(no_effect, X, n_grid = 7)
  expect_identical(dim(observed), c(3L, 7L))
  expect_true(all(observed == 0))
})


test_that("partial objective is finite for valid moments", {
  observed <- fsusieR:::pois_fsusie_partial_objective(
    Y = matrix(c(0, 1, 2, 3), 2),
    scaling = c(1, 2),
    Mu_pm = matrix(0, 2, 2),
    Mu_pv = matrix(0.5, 2, 2),
    B_pm = matrix(0, 2, 2),
    B_pv = matrix(0.25, 2, 2),
    sigma2 = 1
  )
  expect_true(is.finite(observed))
})


test_that("log1p warm start estimates latent sigma2 on the original scale", {
  fitted <- matrix(log1p(4), 2, 2)
  response <- fitted + matrix(c(-1, 1, -1, 1), 2, 2) * sqrt(0.32)
  observed <- fsusieR:::estimate_pois_sigma2_warm_start(
    response_log1p = response,
    fitted_log1p = fitted,
    fitted_variance = matrix(0, 2, 2),
    scaling = c(1, 1),
    fallback_sigma2 = 0.5
  )

  expect_true(observed$used)
  expect_equal(observed$poisson_variance, 0.16, tolerance = 1e-12)
  expect_equal(observed$attenuation_sq, 0.64, tolerance = 1e-12)
  expect_equal(observed$sigma2, 0.25, tolerance = 1e-12)
})


test_that("log1p warm start has an explicit s2 fallback", {
  observed <- fsusieR:::estimate_pois_sigma2_warm_start(
    response_log1p = matrix(0, 2, 2),
    fitted_log1p = matrix(0, 2, 2),
    fitted_variance = matrix(0, 2, 2),
    scaling = c(1, 1),
    fallback_sigma2 = 0.3
  )
  disabled <- fsusieR:::estimate_pois_sigma2_warm_start(
    response_log1p = matrix(1, 2, 2),
    fitted_log1p = matrix(1, 2, 2),
    fitted_variance = matrix(0, 2, 2),
    scaling = c(1, 1),
    fallback_sigma2 = 0.4,
    enabled = FALSE
  )

  expect_false(observed$used)
  expect_equal(observed$sigma2, 0.3)
  expect_false(disabled$used)
  expect_equal(disabled$sigma2, 0.4)
})


test_that("Poisson VGA remains stable for the high-count regression case", {
  observed <- fsusieR::vga_pois_solver(
    init_val = 3,
    x = 15849,
    s = 1,
    beta = 3,
    sigma2 = 3.5736,
    tol = 1e-8
  )

  expect_equal(observed$m, 9.67071235, tolerance = 1e-7)
  expect_equal(observed$v, 6.310178e-5, tolerance = 1e-10)
  expect_true(fsusieR:::vga_pois_solution_is_valid(
    observed, x = 15849, s = 1, beta = 3, sigma2 = 3.5736,
    tol = 1e-7
  ))
})


test_that("Newton and bisection solve mixed-count Poisson VGA rows", {
  counts <- c(0, 1, 100, 15849, 1e6)
  prior_mean <- rep(3, length(counts))

  for (method in c("newton", "bisection")) {
    observed <- fsusieR::vga_pois_solver(
      init_val = prior_mean,
      x = counts,
      s = 1,
      beta = prior_mean,
      sigma2 = 3.5736,
      method = method
    )
    expect_true(fsusieR:::vga_pois_solution_is_valid(
      observed, counts, 1, prior_mean, 3.5736
    ))
    expect_true(all(observed$v > 0 & observed$v <= 3.5736))
    expect_lt(max(observed$m), 20)
  }
})


test_that("Poisson VGA validation rejects the former false solution", {
  false_solution <- list(m = 1453.77587, v = 6.47529602e-5)
  expect_false(fsusieR:::vga_pois_solution_is_valid(
    false_solution, x = 15849, s = 1, beta = 3, sigma2 = 3.5736
  ))
})


test_that("bisection resolves roots much smaller than its absolute tolerance", {
  target <- c(1/25000, 1/1e6)
  observed <- fsusieR::bisection(
    function(value, target) value - target,
    lower = c(0, 0),
    upper = c(4, 4),
    target = target,
    tol = 1e-5
  )
  expect_lt(max(abs(observed/target - 1)), 2e-5)
})


test_that("inverse DWT matrix matches the package coefficient ordering", {
  set.seed(11)
  Y <- matrix(rnorm(24), nrow = 3)
  specifications <- list(
    list(filter.number = 1, family = "DaubExPhase", tolerance = 1e-10),
    list(filter.number = 10, family = "DaubLeAsymm", tolerance = 1e-7)
  )
  for (specification in specifications) {
    transformed <- DWT2(
      Y,
      filter.number = specification$filter.number,
      family = specification$family
    )
    coefficients <- cbind(transformed$D, transformed$C)
    inverse_dwt <- fsusieR:::pois_inverse_dwt_matrix(
      n_grid = ncol(Y),
      filter.number = specification$filter.number,
      family = specification$family
    )
    expect_equal(
      coefficients %*% inverse_dwt,
      Y,
      tolerance = specification$tolerance
    )
  }
  expect_error(fsusieR:::pois_inverse_dwt_matrix(12), "power of two")
})


test_that("structured moments use raw single-effect posterior moments", {
  X <- matrix(c(
    0, 2,
    1, 0,
    2, 1
  ), ncol = 2, byrow = TRUE)
  response <- matrix(c(
    0.5, -1.5,
    1.0, -0.5,
    1.5, -1.0
  ), ncol = 2, byrow = TRUE)
  alpha <- c(0.25, 0.75)
  conditional_mean <- matrix(c(1, 2, -1, 0.5), nrow = 2,
                             byrow = TRUE)
  conditional_variance <- matrix(c(0.2, 0.3, 0.4, 0.1), nrow = 2,
                                 byrow = TRUE)
  mock <- list(
    alpha = list(alpha),
    fitted_wc = list(conditional_mean),
    fitted_wc2 = list(conditional_variance),
    csd_X = apply(X, 2, sd)
  )
  X_scaled <- scale(X)
  expected_effect_mean <- X_scaled %*%
    sweep(conditional_mean, 1, alpha, "*")
  expected_effect_second <- X_scaled^2 %*%
    sweep(conditional_mean^2 + conditional_variance, 1, alpha, "*")

  observed <- fsusieR:::compute_pois_structured_moments(
    susiF.obj = mock,
    X = X,
    response = response,
    inverse_dwt = diag(2)
  )

  expect_equal(
    observed$mean,
    sweep(expected_effect_mean, 2, colMeans(response), "+"),
    tolerance = 1e-12
  )
  expect_equal(
    observed$variance,
    pmax(expected_effect_second - expected_effect_mean^2, 0),
    tolerance = 1e-12
  )
  expect_equal(
    observed$coefficient_mean,
    sweep(sweep(conditional_mean, 1, alpha, "*"),
          1, apply(X, 2, sd), "/"),
    tolerance = 1e-12
  )
  expect_equal(
    fsusieR:::build_B_pm(response, X, observed$coefficient_mean),
    observed$mean,
    tolerance = 1e-12
  )
})


test_that("uninitialized zero-mass effects do not alter structured moments", {
  mock <- list(
    L = 1L,
    alpha = list(c(0, 0)),
    fitted_wc = list(matrix(0, nrow = 2, ncol = 2)),
    fitted_wc2 = list(matrix(1, nrow = 2, ncol = 2)),
    csd_X = c(1, 1)
  )
  X <- matrix(c(0, 1, 1, 0, 2, 1), nrow = 3, byrow = TRUE)
  response <- matrix(1:6, nrow = 3)

  observed <- fsusieR:::compute_pois_structured_moments(
    mock, X, response, diag(2)
  )
  expect_equal(observed$mean,
               matrix(colMeans(response), nrow = 3, ncol = 2, byrow = TRUE))
  expect_equal(observed$variance, matrix(0, nrow = 3, ncol = 2))
  expect_equal(observed$coefficient_mean, matrix(0, nrow = 2, ncol = 2))
})


test_that("Poisson sigma2 M-step includes both posterior variances", {
  Mu_pm <- matrix(c(0, 1, 2, 3), nrow = 2)
  B_pm <- matrix(c(0.5, 0.5, 1.5, 2.5), nrow = 2)
  Mu_pv <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
  B_pv <- matrix(c(0.4, 0.3, 0.2, 0.1), nrow = 2)
  observed <- fsusieR:::pois_sigma2_update(Mu_pm, Mu_pv, B_pm, B_pv)

  expect_equal(observed$residual_component, mean((Mu_pm - B_pm)^2))
  expect_equal(observed$latent_variance_component, mean(Mu_pv))
  expect_equal(observed$structured_variance_component, mean(B_pv))
  expect_equal(
    observed$sigma2,
    mean((Mu_pm - B_pm)^2 + Mu_pv + B_pv),
    tolerance = 1e-12
  )
})


test_that("susiF can hold an input-scale residual variance fixed", {
  set.seed(12)
  X <- cbind(rnorm(30), rnorm(30))
  Y <- tcrossprod(X[, 1], c(0, 0.4, 0.8, 0.4, 0, -0.3, -0.6, -0.3)) +
    matrix(rnorm(30 * 8, sd = 0.5), nrow = 30)

  observed <- suppressWarnings(susiF(
    Y = Y,
    X = X,
    L = 1,
    L_start = 1,
    prior = "mixture_normal",
    maxit = 2,
    standardize_Y = FALSE,
    residual_variance = 0.37,
    post_processing = "none",
    filter_cs = FALSE,
    verbose = FALSE
  ))

  expect_equal(observed$sigma2, 0.37, tolerance = 0)
  expect_true(observed$residual_variance_fixed)
  expect_false(observed$standardize_Y)
  expect_equal(observed$csd_Y_pos, rep(1, ncol(Y)))
  expect_equal(observed$csd_Yf, rep(1, ncol(Y)))
})


test_that("Poisson fSuSiE uses the preceding outer sigma2 in each inner fit", {
  set.seed(21)
  N <- 20
  Tt <- 8
  X <- matrix(rnorm(N * 3), nrow = N)
  eta <- outer(X[, 1], c(0, 0, 0.35, 0.5, 0.35, 0, -0.2, -0.1))
  Y <- matrix(rpois(N * Tt, lambda = exp(eta)), nrow = N)

  observed <- Pois_fSuSiE(
    Y = Y,
    X = X,
    L = 1,
    maxit_outer = 2,
    maxit_inner = 2,
    post_processing = "TI",
    filter.number = 1,
    family = "DaubExPhase",
    warm_start_sigma2 = FALSE,
    s2 = 0.4,
    stable_iterations = 3,
    verbose = FALSE
  )

  expect_equal(observed$n_iter, 2L)
  expect_equal(observed$convergence_trace$inner_sigma2[1],
               observed$sigma2_initial, tolerance = 0)
  expect_equal(observed$convergence_trace$inner_sigma2[2],
               observed$convergence_trace$sigma2[1], tolerance = 0)
  expect_equal(observed$sigma2_used_for_final_susif,
               tail(observed$convergence_trace$inner_sigma2, 1),
               tolerance = 0)
  expect_true(observed$susiF.obj$residual_variance_fixed)
  expect_false(observed$susiF.obj$standardize_Y)
  expect_true(all(is.finite(observed$B_pm)))
  expect_true(all(observed$B_pv >= 0))
  expect_gt(observed$sigma2, 0)
  expect_equal(
    fsusieR:::build_B_pm(observed$Mu_pm, X, observed$est_effect_fm),
    observed$B_pm,
    tolerance = 1e-8
  )
})
