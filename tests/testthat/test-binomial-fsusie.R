test_that("binomial inputs retain success counts and trial totals", {
  Y <- matrix(c(0, 1, 2, 3, 1, 0, 2, 1), nrow = 4)
  X <- cbind(c(0, 1, 2, 1), c(2, 0, 1, 2))

  scalar <- fsusieR:::validate_binomial_fsusie_inputs(Y, 3, X, 2)
  by_sample <- fsusieR:::validate_binomial_fsusie_inputs(
    Y, c(3, 4, 5, 6), X, 2
  )
  expect_equal(scalar$trials, matrix(3, 4, 2))
  expect_equal(by_sample$trials[, 1], c(3, 4, 5, 6))
  expect_error(
    fsusieR:::validate_binomial_fsusie_inputs(Y, 2, X, 2),
    "no larger"
  )
  expect_error(
    fsusieR:::validate_binomial_fsusie_inputs(Y, 3.5, X, 2),
    "positive integer"
  )
})


test_that("Gauss-Hermite logistic-normal moments respect symmetry", {
  variance <- c(0.01, 0.5, 2, 10)
  observed <- fsusieR:::binomial_logistic_normal_terms(
    rep(0, length(variance)), variance, quadrature_order = 30
  )
  expect_equal(observed$probability, rep(0.5, length(variance)),
               tolerance = 1e-12)
  expect_true(all(observed$probability_variance > 0))
  expect_true(all(observed$softplus_log_variance_derivative > 0))
})


test_that("binomial quadrature gradient matches finite differences", {
  m <- c(-1.2, 0.4, 2)
  log_v <- log(c(0.15, 0.8, 1.4))
  x <- c(0, 4, 9)
  trials <- c(3, 8, 10)
  beta <- c(-0.4, 0.1, 0.8)
  sigma2 <- c(0.7, 1.1, 2)
  order <- 30L
  objective <- function(theta) {
    current_m <- theta[1:3]
    current_v <- exp(theta[4:6])
    terms <- fsusieR:::binomial_logistic_normal_terms(
      current_m, current_v, order
    )
    sum(x * current_m - trials * terms$softplus -
          ((current_m - beta)^2 + current_v) / (2 * sigma2) +
          log(current_v) / 2)
  }
  theta <- c(m, log_v)
  epsilon <- 1e-6
  numerical <- vapply(seq_along(theta), function(index) {
    step <- rep(0, length(theta))
    step[index] <- epsilon
    (objective(theta + step) - objective(theta - step)) / (2 * epsilon)
  }, numeric(1))
  terms <- fsusieR:::binomial_logistic_normal_terms(m, exp(log_v), order)
  analytical <- c(
    x - trials * terms$probability - (m - beta) / sigma2,
    0.5 - exp(log_v) / (2 * sigma2) -
      trials * terms$softplus_log_variance_derivative
  )
  expect_equal(analytical, numerical, tolerance = 2e-6)
})


test_that("binomial VGA solves interior and boundary count cases", {
  successes <- c(0, 1, 5, 9, 10, 30)
  totals <- c(1, 2, 10, 10, 10, 30)
  prior_mean <- c(-1, 0, 0.5, 1, 2, -2)
  prior_variance <- c(0.2, 0.5, 1, 2, 4, 3)
  observed <- vga_binomial_solver(
    init_val = prior_mean,
    x = successes,
    trials = totals,
    beta = prior_mean,
    sigma2 = prior_variance,
    quadrature_order = 30
  )
  expect_true(fsusieR:::vga_binomial_solution_is_valid(
    observed, successes, totals, prior_mean, prior_variance,
    quadrature_order = 30
  ))
  expect_true(all(observed$v > 0))
  expect_true(all(observed$v <= prior_variance))
  expect_gt(observed$m[5], observed$m[1])
})


test_that("zero-trial VGA coordinates return their Gaussian prior", {
  observed <- vga_binomial_solver(
    init_val = c(10, 0),
    x = c(0, 1),
    trials = c(0, 2),
    beta = c(-0.5, 0.2),
    sigma2 = c(0.7, 1.1),
    quadrature_order = 20
  )
  expect_equal(observed$m[1], -0.5)
  expect_equal(observed$v[1], 0.7)
})


test_that("binomial fSuSiE keeps a requested splitting variance fixed", {
  set.seed(210)
  n <- 24L
  tt <- 8L
  X <- matrix(rnorm(n * 3), nrow = n)
  effect <- c(0, 0, 0.5, 0.8, 0.5, 0, -0.3, -0.2)
  eta <- -0.8 + outer(X[, 1], effect)
  trials <- matrix(sample(12:25, n * tt, replace = TRUE), nrow = n)
  Y <- matrix(
    rbinom(n * tt, size = as.numeric(trials),
           prob = stats::plogis(as.numeric(eta))),
    nrow = n
  )

  observed <- suppressWarnings(Binomial_fSuSiE(
    Y = Y,
    trials = trials,
    X = X,
    L = 1,
    maxit_outer = 2,
    maxit_inner = 2,
    stable_iterations = 3,
    s2 = 0.4,
    estimate_sigma2 = FALSE,
    post_processing = "TI",
    filter.number = 1,
    family = "DaubExPhase",
    quadrature_order = 20,
    verbose = FALSE
  ))

  expect_identical(observed$family, "binomial")
  expect_equal(observed$sigma2, 0.4, tolerance = 0)
  expect_false(observed$estimate_sigma2)
  expect_true(observed$susiF.obj$residual_variance_fixed)
  expect_false(observed$susiF.obj$standardize_Y)
  expect_equal(dim(observed$posterior_probability_mean), c(n, tt))
  expect_true(all(observed$posterior_probability_mean > 0 &
                    observed$posterior_probability_mean < 1))
  expect_equal(
    observed$posterior_success_mean,
    observed$trials * observed$posterior_probability_mean,
    tolerance = 1e-12
  )
  expect_equal(nrow(observed$sigma2_subcycle_trace), 0L)
  expect_equal(observed$solver_failures$row, integer(0))
})


test_that("binomial fSuSiE estimates sigma2 with replicated trials", {
  set.seed(211)
  n <- 20L
  tt <- 8L
  X <- matrix(rnorm(n * 3), nrow = n)
  effect <- c(0, 0.25, 0.5, 0.25, 0, -0.2, -0.3, 0)
  eta <- -0.5 + outer(X[, 1], effect) + matrix(rnorm(n * tt, sd = 0.25), n)
  trials <- matrix(20, nrow = n, ncol = tt)
  Y <- matrix(
    rbinom(n * tt, 20, stats::plogis(as.numeric(eta))),
    nrow = n
  )

  observed <- suppressWarnings(Binomial_fSuSiE(
    Y = Y,
    trials = trials,
    X = X,
    L = 1,
    maxit_outer = 2,
    maxit_inner = 2,
    stable_iterations = 3,
    s2 = 0.3,
    warm_start_sigma2 = FALSE,
    sigma2_subcycles = 2,
    post_processing = "TI",
    filter.number = 1,
    family = "DaubExPhase",
    quadrature_order = 20,
    verbose = FALSE
  ))

  expect_true(observed$estimate_sigma2)
  expect_gt(observed$sigma2, 0)
  expect_true(all(is.finite(observed$Mu_pm)))
  expect_true(all(observed$Mu_pv > 0))
  expect_identical(observed$convergence_trace$sigma2_subcycles, c(1L, 2L))
  expect_identical(observed$sigma2_subcycle_trace$iteration,
                   c(1L, 2L, 2L))
  expect_equal(observed$solver_failures$row, integer(0))
})
