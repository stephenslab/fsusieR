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
