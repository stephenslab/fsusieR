test_that("fSuSiE restores every SNP-indexed field to named input columns", {
  obj <- structure(list(
    P = 2L,
    L = 1L,
    alpha = list(c(0.8, 0.2)),
    pip = c(0.8, 0.2),
    lBF = list(c(2, 1)),
    fitted_wc = list(matrix(1:6, nrow = 2)),
    fitted_wc2 = list(matrix(7:12, nrow = 2)),
    alpha_hist = list(list(c(0.8, 0.2))),
    csd_X = c(1, 3),
    d = c(9, 9),
    cs = list(1:2),
    cov_lev = 0.95,
    fitted_func = list(rep(0, 3))
  ), class = "susiF")

  restored <- fsusieR:::rename_format_output(
    obj,
    names_colX = paste0("snp", 1:4),
    tidx = c(2L, 4L),
    kept_index = c(1L, 3L),
    original_P = 4L
  )

  expect_equal(restored$P, 4L)
  expect_equal(restored$n_cs, 1L)
  expect_equal(restored$cs_size, 2L)
  expect_equal(restored$fitted_P, 2L)
  expect_equal(restored$variable_index, c(1L, 3L))
  expect_equal(restored$removed_variable_index, c(2L, 4L))
  expect_equal(unname(restored$alpha[[1]]), c(0.8, 0, 0.2, 0))
  expect_equal(unname(restored$pip), c(0.8, 0, 0.2, 0))
  expect_equal(unname(restored$lBF[[1]]), c(2, -Inf, 1, -Inf))
  expect_equal(restored$cs[[1]], structure(c(1L, 3L),
                                           names = c("snp1", "snp3")))
  expect_equal(nrow(restored$fitted_wc[[1]]), 4L)
  expect_equal(unname(restored$fitted_wc[[1]][c(2, 4), ]), matrix(0, 2, 3))
  expect_equal(names(restored$pip), paste0("snp", 1:4))
})

test_that("fSuSiE restores unnamed input columns by numeric position", {
  obj <- structure(list(
    P = 2L,
    L = 1L,
    alpha = list(c(0.7, 0.3)),
    pip = c(0.7, 0.3),
    lBF = list(c(1, 0)),
    fitted_wc = list(matrix(1:4, nrow = 2)),
    fitted_wc2 = list(matrix(5:8, nrow = 2)),
    alpha_hist = list(),
    csd_X = c(1, 1),
    d = c(9, 9),
    cs = list(),
    cov_lev = 0.95,
    fitted_func = list(rep(0, 2))
  ), class = "susiF")

  restored <- fsusieR:::rename_format_output(
    obj,
    names_colX = NULL,
    tidx = c(1L, 3L),
    kept_index = c(2L, 4L),
    original_P = 4L
  )

  expect_null(names(restored$pip))
  expect_equal(restored$cs[[1]], c(2L, 4L))
  expect_equal(unname(restored$pip), c(0, 0.7, 0, 0.3))
  expect_equal(nrow(restored$fitted_wc[[1]]), 4L)
})
