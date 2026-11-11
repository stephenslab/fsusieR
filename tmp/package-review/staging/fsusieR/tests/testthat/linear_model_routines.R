library(testthat)
library(matrixStats)
lm_reference <- function(Y, X, idx = NULL) {
  if (!is.null(idx)) {
    Y <- Y[idx, , drop = FALSE]
    X <- X[idx, , drop = FALSE]
  }

  p <- ncol(X)
  J <- ncol(Y)

  Bhat <- matrix(NA_real_, p, J)
  Shat <- matrix(NA_real_, p, J)

  for (j in seq_len(J)) {
    for (k in seq_len(p)) {
      fit <- lm(Y[, j] ~ X[, k] - 1)
      s   <- summary(fit)
      Bhat[k, j] <- coef(s)[1, "Estimate"]
      Shat[k, j] <- coef(s)[1, "Std. Error"]
    }
  }

  list(Bhat = Bhat, Shat = Shat)
}
generate_test_data <- function(n = 80, p = 12, J = 20, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * J), n, J)

  # center & scale to match assumptions
  X <- colScale(X)
  Y <- colScale(Y)

  list(X = X, Y = Y)
}
test_that("Bhat matches lm without subsetting", {
  dat <- generate_test_data()
  fast <- cal_Bhat_Shat(dat$Y, dat$X)
  ref  <- lm_reference(dat$Y, dat$X)

  expect_equal(
    fast$Bhat,
    ref$Bhat,
    tolerance = 1e-10
  )
})
test_that("Shat matches lm standard errors", {
  dat <- generate_test_data()
  fast <- cal_Bhat_Shat(dat$Y, dat$X)
  ref  <- lm_reference(dat$Y, dat$X)

  expect_equal(
    fast$Shat,
    ref$Shat,
    tolerance = 1e-10
  )
})
test_that("ind_analysis vector matches lm on subset", {
  dat <- generate_test_data()
  idx <- sample(seq_len(nrow(dat$X)), 50)

  fast <- cal_Bhat_Shat(dat$Y, dat$X, ind_analysis = idx)
  ref  <- lm_reference(dat$Y, dat$X, idx = idx)

  expect_equal(fast$Bhat, ref$Bhat, tolerance = 1e-10)
  expect_equal(fast$Shat, ref$Shat, tolerance = 1e-10)
})
test_that("ind_analysis list works per-response", {
  dat <- generate_test_data(J = 6)
  ind_list <- lapply(seq_len(ncol(dat$Y)), function(j) {
    sample(seq_len(nrow(dat$Y)), 60)
  })

  fast <- cal_Bhat_Shat(dat$Y, dat$X, ind_analysis = ind_list)

  for (j in seq_len(ncol(dat$Y))) {
    ref_j <- lm_reference(dat$Y[, j, drop = FALSE],
                          dat$X,
                          idx = ind_list[[j]])

    expect_equal(
      fast$Bhat[, j],
      ref_j$Bhat[, 1],
      tolerance = 1e-10
    )

    expect_equal(
      fast$Shat[, j],
      ref_j$Shat[, 1],
      tolerance = 1e-10
    )
  }
})

test_that("lowc_wc masks coefficients and sets Shat to 1", {
  dat <- generate_test_data()
  lowc <- c(2, 5)

  fast <- cal_Bhat_Shat(dat$Y, dat$X, lowc_wc = lowc)

  expect_true(all(fast$Bhat[, lowc] == 0))
  expect_true(all(fast$Shat[, lowc] == 1))
})
test_that("large problem remains numerically stable", {
  dat <- generate_test_data(n = 300, p = 80, J = 50)
  fast <- cal_Bhat_Shat(dat$Y, dat$X)

  expect_true(all(is.finite(fast$Bhat)))
  expect_true(all(is.finite(fast$Shat)))
})
