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
