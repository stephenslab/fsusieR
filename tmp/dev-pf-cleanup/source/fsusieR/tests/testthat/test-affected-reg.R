test_that("HMM regions skip empty effects and keep the original effect index", {
  obj <- list(
    outing_grid = c(100, 105, 130, 131, 200, 205, 250, 300),
    cs = list(1L, 2L, 3L),
    lfsr_func = list(rep(0.5, 8),
                     c(0.01, 0.5, 0.01, 0.01, 0.5, 0.5, 0.5, 0.01),
                     rep(0.5, 8))
  )

  expect_warning(reg <- affected_reg(obj), NA)
  rownames(reg) <- NULL
  expect_equal(reg, data.frame(CS = c(2, 2, 2),
                               Start = c(100, 130, 300),
                               End = c(100, 131, 300)))
})

test_that("SMASH regions exclude zero and skip empty effects", {
  empty_band <- rbind(up = rep(1, 8), low = rep(-1, 8))
  obj <- list(
    outing_grid = c(100, 105, 130, 131, 200, 205, 250, 300),
    cs = list(1L, 2L, 3L),
    cred_band = list(empty_band,
                     rbind(up = c(-1, 0, 2, 2, 2, 0, -1, 1),
                           low = c(-2, -1, 1, 1, -1, -1, -2, 0)),
                     empty_band)
  )

  expect_warning(reg <- affected_reg(obj), NA)
  rownames(reg) <- NULL
  expect_equal(reg, data.frame(CS = c(2, 2, 2),
                               Start = c(100, 130, 250),
                               End = c(100, 131, 250)))
})

test_that("no affected regions return an empty three-column data frame", {
  expected <- data.frame(CS = integer(), Start = numeric(), End = numeric())
  for (n_effects in 0:2) {
    obj <- list(outing_grid = 1:4, cs = as.list(seq_len(n_effects)))

    obj$lfsr_func <- rep(list(rep(0.5, 4)), n_effects)
    expect_warning(reg <- affected_reg(obj), NA)
    expect_identical(reg, expected)

    obj$lfsr_func <- NULL
    obj$cred_band <- rep(list(rbind(up = rep(1, 4), low = rep(-1, 4))),
                         n_effects)
    expect_warning(reg <- affected_reg(obj), NA)
    expect_identical(reg, expected)
  }
})

test_that("HMM uses the requested strict LFSR threshold", {
  obj <- list(outing_grid = 1:4, cs = list(1L),
              lfsr_func = list(c(0.05, 0.04, 0.1, 0.5)))

  expect_equal(affected_reg(obj), data.frame(CS = 1L, Start = 2L, End = 2L))
  expect_equal(affected_reg(obj, lfsr_thresh = 0.1),
               data.frame(CS = 1L, Start = 1L, End = 2L))
  expect_warning(reg <- affected_reg(obj, lfsr_thresh = 0), NA)
  expect_equal(nrow(reg), 0L)
})

test_that("SMASH keeps adjacent positive and negative regions separate", {
  obj <- list(outing_grid = 1:4, cs = list(1L),
              cred_band = list(rbind(up = c(2, 2, -1, -1),
                                     low = c(1, 1, -2, -2))))

  expect_warning(reg <- affected_reg(obj), NA)
  rownames(reg) <- NULL
  expect_equal(reg, data.frame(CS = c(1L, 1L), Start = c(1L, 3L),
                               End = c(2L, 4L)))
})
