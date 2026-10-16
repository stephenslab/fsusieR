test_that("SMASH retains uncertainty in a constant effect, including reflected grids", {
  set.seed(131)
  n <- 40L
  x <- as.numeric(scale(seq_len(n)))
  # Orthogonal residual columns with sample variance one give the same
  # standardized slope and SE at every position, so all wavelet details vanish.
  q <- qr.Q(qr(cbind(1, x, matrix(rnorm(n * 8L), nrow = n))))
  residual <- q[, 3:10] * sqrt(n - 1)

  for (m in c(8L, 5L)) {
    scales <- seq(0.5, 2, length.out = m)
    # Five positions are padded to [2, 1, 1, 2, 3, 4, 5, 5], then mirrored.
    weights <- if (m == 8L) rep(1 / 8, 8) else c(2, 2, 1, 1, 2) / 8
    expected_var <- sum(weights^2) / (n - 1) * scales^2

    for (beta in c(0.02, 1)) {
      y <- sweep(beta * x + residual[, seq_len(m)], 2, scales, "*")
      fit <- univariate_smash_regression(y, matrix(x, ncol = 1L))

      expect_equal(fit$effect_estimate, beta * scales, tolerance = 1e-10)
      expect_equal(fit$fitted_var, expected_var, tolerance = 1e-10)
      expect_equal(unname(fit$cred_band[1, ] - fit$cred_band[2, ]),
                   2 * qnorm(0.975) * sqrt(expected_var), tolerance = 1e-10)
      if (beta == 0.02) {
        expect_true(all(fit$cred_band[2, ] < 0 & fit$cred_band[1, ] > 0))
      } else {
        expect_true(all(fit$cred_band[2, ] > 0))
      }
    }
  }
})

test_that("noise-only SMASH bands do not collapse onto a random nonzero mean", {
  set.seed(13)
  n <- 80L
  p <- 6L
  m <- 32L
  x <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  y <- matrix(rnorm(n * m, sd = 0.5), nrow = n)

  fit <- susiF(y, x, L = 2, post_processing = "smash", verbose = FALSE)

  expect_length(fit$cred_band, 1L)
  expect_true(all(is.finite(fit$fitted_var[[1]])))
  expect_true(all(fit$fitted_var[[1]] > 0))
  # Inspect the bands directly to check uncertainty independently of region extraction.
  expect_false(any(fit$cred_band[[1]][1, ] < 0 |
                     fit$cred_band[[1]][2, ] > 0))
})
