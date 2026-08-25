test_that("susiF keeps internally requested HMM post-processing quiet", {
  postprocess_verbose <- logical(0)

  local_mocked_bindings(
    update_cal_pip = function(obj, ...) obj,
    name_cs = function(obj, ...) obj,
    update_cal_fit_func = function(obj, ..., verbose) {
      postprocess_verbose <<- c(postprocess_verbose, verbose)
      obj
    },
    update_cal_indf = function(obj, ...) obj,
    rename_format_output = function(obj, ...) obj,
    .package = "fsusieR"
  )

  obj <- structure(list(cs = list(1L)), class = c("susiF", "list"))
  common_args <- list(
    obj = obj,
    Y = matrix(0, nrow = 2L, ncol = 2L),
    X = matrix(c(0, 1), ncol = 1L),
    indx_lst = list(1L, 2L),
    outing_grid = 1:2,
    filter_cs = FALSE,
    tidx = integer(0),
    kept_index = 1L,
    original_P = 1L,
    pos = 1:2,
    verbose = TRUE
  )

  do.call(out_prep.susiF,
          c(common_args, list(post_processing = "HMM")))
  do.call(out_prep.susiF,
          c(common_args, list(post_processing = "TI")))

  expect_identical(postprocess_verbose, c(FALSE, TRUE))
})


test_that("single-effect HMM regression forwards its verbosity flag", {
  hmm_verbose <- logical(0)
  local_mocked_bindings(
    fit_hmm = function(x, sd, verbose, ...) {
      hmm_verbose <<- c(hmm_verbose, verbose)
      list(posterior = list(lfsr = rep(1, length(x)),
                            mean = rep(0, length(x))))
    },
    .package = "fsusieR"
  )

  obj <- structure(
    list(cs = list(1L), alpha = list(c(1, 0))),
    class = c("susiF", "list")
  )
  set.seed(41)
  x <- cbind(seq_len(12), rev(seq_len(12)))
  y <- outer(seq_len(12), c(0.1, 0.2, 0.3, 0.4)) +
    matrix(rnorm(48, sd = 0.5), nrow = 12)

  HMM_regression(obj, y, x, verbose = FALSE)
  invisible(capture.output(HMM_regression(obj, y, x, verbose = TRUE)))

  expect_identical(hmm_verbose, c(FALSE, TRUE))
})
