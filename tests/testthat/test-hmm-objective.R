test_that("state deletion cannot decrease the fitted HMM objective", {
  # Previously, the default 0.05 pruning allowance removed the useful second
  # state at step 5 and lowered the objective by approximately 0.0025.
  for (inference in c("exact", "variational")) {
    fit <- fit_ash_hmm(
      y = 0.26, se = 1, mu = c(0, 0.5), prior_sd = 0,
      estimate_init = TRUE, learn_state_means = FALSE, maxiter = 10,
      inference = inference)
    expect_true(all(diff(fit$history$objective) >= 0))
    expect_equal(tail(fit$history$states, 1), 2)
    expect_equal(fit$training_objective, tail(fit$history$objective, 1))
  }
})

test_that("objective-safe deletion still removes unsupported HMM states", {
  fit <- fit_ash_hmm(
    y = c(0.1, -0.1, 0), se = rep(0.1, 3),
    mu = c(0, 10, -10), prior_sd = 0,
    transition_prior = c(1.5, 1, 1),
    maxiter = 1, prune_start = 1, prune_every = 1,
    prune_max_fraction = 0.9)
  expect_equal(tail(fit$history$states, 1), 1)
  expect_true(all(diff(fit$history$objective) >= 0))
  expect_equal(fit$training_objective, tail(fit$history$objective, 1))
})
