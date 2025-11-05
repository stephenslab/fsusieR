test_that("Rcpp C++ implementations match base R and matrixStats", {
  skip_if_not_installed("matrixStats")
  library(matrixStats)

  set.seed(1)
  X <- matrix(rnorm(1e7), ncol = 1000)  # 10,000 x 1000 matrix

  # ---- colSums ----
  ref_sums <- matrixStats::colSums2(X)
  cpp_sums <- fsusieR:::colSumsCpp(X)
  expect_equal(cpp_sums, ref_sums, tolerance = 1e-12,
               info = "colSumsCpp should match matrixStats::colSums2")

  # ---- colVars ----
  ref_vars <- matrixStats::colVars(X)
  cpp_vars <- fsusieR:::colVarsCpp(X)
  expect_equal(cpp_vars, ref_vars, tolerance = 1e-12,
               info = "colVarsCpp should match matrixStats::colVars")

  # ---- covCpp ----
  smallX <- X[, 1:10]
  ref_cov <- cov(smallX)
  cpp_cov <- fsusieR:::covCpp(smallX)
  expect_equal(cpp_cov, ref_cov, tolerance = 1e-12,
               info = "covCpp should match base::cov for 10 columns")

  # ---- check positive semi-definiteness for covCpp ----
  eigen_cpp <- eigen(cpp_cov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigen_cpp > -1e-10),
              info = "covCpp should produce a positive semidefinite matrix")
})
