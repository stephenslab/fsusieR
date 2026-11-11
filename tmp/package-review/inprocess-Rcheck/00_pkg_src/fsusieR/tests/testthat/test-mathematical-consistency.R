test_that("cal_Bhat_Shat uses fixed sigma2 without breaking legacy calls", {
  X <- scale(matrix(c(
    0, 1,
    1, 0,
    2, 1,
    3, 3,
    4, 2
  ), ncol = 2, byrow = TRUE))
  Y <- scale(matrix(c(
    0.2, -0.1,
    0.5,  0.4,
   -0.3,  0.2,
    0.7, -0.5,
   -0.1,  0.8
  ), ncol = 2, byrow = TRUE))

  fixed <- cal_Bhat_Shat(Y, X, sigma2 = c(0.5, 2))
  expected_shat <- sqrt(outer(1 / colSums(X^2), c(0.5, 2)))
  expect_equal(fixed$Shat, expected_shat, tolerance = 1e-12)

  # The historical third positional argument was v1. It must remain ignored
  # rather than accidentally being interpreted as sigma2.
  legacy <- cal_Bhat_Shat(Y, X, rep(1, nrow(X)))
  fallback <- cal_Bhat_Shat(Y, X)
  expect_equal(legacy, fallback, tolerance = 1e-12)
})

test_that("posterior coefficients skipped by the approximation have zero SD", {
  prior <- structure(list(NULL), class = "mixture_normal")
  Bhat <- matrix(c(0.1, -0.2, 0.3, -0.4), nrow = 2)
  Shat <- matrix(1, nrow = 2, ncol = 2)

  observed <- post_mat_sd(
    prior,
    Bhat,
    Shat,
    lBF = c(0, 0),
    lowc_wc = 2L,
    indx_lst = list(1:2),
    e = 1
  )
  expect_equal(observed, matrix(0, nrow = 2, ncol = 2))
})

test_that("L equals one partial residual does not subtract its own effect", {
  D <- matrix(c(0.2, -0.1, 0.4, 0.7, -0.3, 0.5), ncol = 2)
  C <- c(0.1, -0.2, 0.3)
  X <- matrix(c(1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  obj <- structure(list(
    L = 1L,
    G_prior = structure(list(NULL), class = "mixture_normal"),
    alpha = list(c(0.6, 0.4)),
    fitted_wc = list(matrix(c(0.5, -0.2, 0.1, 0.3, 0.4, -0.1),
                            nrow = 2))
  ), class = "susiF")

  observed <- cal_partial_resid(
    obj,
    l = 0L,
    X = X,
    D = D,
    C = C,
    indx_lst = list(1:2, 3L)
  )
  expect_equal(unname(observed), unname(cbind(D, C)))
})

test_that("fSuSiE objective is expected log likelihood minus exact KL", {
  X <- matrix(c(1, 0, 0, 2, 1, -1, 2, 1), ncol = 2, byrow = TRUE)
  D <- matrix(c(0.2, -0.3, 0.5, 0.7), ncol = 1)
  C <- c(0.1, 0.4, -0.2, 0.3)
  Y <- cbind(D, C)
  obj <- structure(list(
    L = 1L,
    G_prior = structure(list(NULL), class = "mixture_normal"),
    alpha = list(c(0.6, 0.4)),
    fitted_wc = list(matrix(c(0.5, -0.2, 0.1, 0.3), nrow = 2)),
    fitted_wc2 = list(matrix(c(0.10, 0.05, 0.04, 0.02), nrow = 2)),
    sigma2 = 0.7,
    lBF = list(c(log(2), log(0.5))),
    KL = 0.25
  ), class = "susiF")

  EF <- get_post_F(obj, 1L)
  EF2 <- get_post_F2(obj, 1L)
  expected_kl <-
    sum(Y * (X %*% EF)) / obj$sigma2 -
    sum(colSums(X^2) * rowSums(EF2)) / (2 * obj$sigma2) -
    log(mean(exp(obj$lBF[[1]])))

  observed_kl <- cal_KL_l(
    obj,
    l = 1L,
    X = X,
    D = D,
    C = C,
    indx_lst = list(1L, 2L)
  )
  expect_equal(observed_kl, expected_kl, tolerance = 1e-12)

  expected_signal_ss <- sum(colSums(X^2) * rowSums(EF2))
  expected_post <-
    -length(Y) / 2 * log(2 * pi * obj$sigma2) -
    (sum(Y^2) - 2 * sum(Y * (X %*% EF)) + expected_signal_ss) /
      (2 * obj$sigma2)
  expect_equal(
    loglik_SFR_post(obj, 1L, Y, X),
    expected_post,
    tolerance = 1e-12
  )

  expect_equal(
    get_objective(obj, Y, X, D, C, list(1L, 2L)),
    Eloglik(obj, Y, X) - obj$KL,
    tolerance = 1e-12
  )
})
