make_dynamic_susif <- function(L=3L, P=3L) {
  alpha <- lapply(seq_len(L), function(l) {
    value <- seq_len(P) + l
    value / sum(value)
  })
  prior <- structure(list(), class="mixture_normal")
  structure(
    list(
      L=L,
      L_max=max(L, P),
      P=P,
      alpha=alpha,
      lBF=lapply(seq_len(L), function(l) seq_len(P) + l),
      fitted_wc=lapply(seq_len(L), function(l) matrix(l, nrow=P, ncol=4)),
      fitted_wc2=lapply(seq_len(L), function(l) matrix(l + 1, nrow=P, ncol=4)),
      cs=lapply(seq_len(L), function(l) as.integer((l - 1L) %% P + 1L)),
      est_pi=lapply(seq_len(L), function(l) c(0.5, 0.5)),
      est_sd=lapply(seq_len(L), function(l) c(0, 1)),
      cred_band=lapply(seq_len(L), function(l) matrix(0, nrow=2, ncol=4)),
      lfsr_wc=lapply(seq_len(L), function(l) rep(1, 4)),
      KL=seq_len(L),
      pip=rep(0, P),
      G_prior=prior,
      tol_null_prior=0.001,
      cov_lev=0.8,
      n_wac=4L,
      csd_X=rep(1, P),
      d=rep(9, P),
      alpha_hist=list(alpha),
      n_expand=0L,
      greedy=TRUE,
      greedy_backfit_update=FALSE,
      ELBO=numeric()
    ),
    class="susiF"
  )
}


test_that("fSuSiE discard keeps every effect-indexed field aligned", {
  obj <- make_dynamic_susif()
  retained_alpha <- obj$alpha[c(1, 3)]

  result <- discard_cs(obj, cs=2L, out_prep=TRUE)

  expect_equal(result$L, 2L)
  expect_equal(result$alpha, retained_alpha)
  for (field in c(
    "alpha", "lBF", "fitted_wc", "fitted_wc2", "cs",
    "est_pi", "est_sd", "cred_band", "lfsr_wc", "KL"
  )) {
    expect_length(result[[field]], result$L)
  }
  expect_equal(result$KL, c(1L, 3L))
  expect_length(result$pip, result$P)

  protected <- discard_cs(obj, cs=seq_len(obj$L), out_prep=TRUE)
  expect_equal(protected$L, 1L)
  expect_length(protected$cs, 1L)

  expect_error(discard_cs(obj, cs=0L), "invalid effect index")
  expect_error(discard_cs(obj, cs=4L), "invalid effect index")
})


test_that("fSuSiE expansion respects P and keeps dynamic fields aligned", {
  obj <- make_dynamic_susif(L=1L, P=3L)
  obj$L_max <- 5L

  expanded <- expand_susiF_obj(obj, L_extra=10L)

  expect_equal(expanded$L, 3L)
  for (field in c(
    "alpha", "lBF", "fitted_wc", "fitted_wc2", "cs",
    "est_pi", "est_sd", "cred_band", "lfsr_wc", "KL"
  )) {
    expect_length(expanded[[field]], expanded$L)
  }
  expect_true(all(lengths(expanded$alpha) == expanded$P))
  expect_equal(expand_susiF_obj(expanded, 2L)$L, 3L)
  expect_error(expand_susiF_obj(obj, -1L), "non-negative integer")

  finalized <- obj
  finalized$column_index_space <- "original"
  finalized$fitted_P <- finalized$P
  expect_error(expand_susiF_obj(finalized, 1L), "after original-column")
})


test_that("fSuSiE maps dynamic output without recalculating CS membership", {
  obj <- make_dynamic_susif()
  obj$cs <- list(c(1L, 3L), c(1L, 2L), 2L)
  obj <- discard_cs(obj, cs=2L, out_prep=TRUE)
  obj$fitted_func <- lapply(seq_len(obj$L), function(l) rep(l, 4))
  obj$fitted_var <- lapply(seq_len(obj$L), function(l) rep(0, 4))

  restored <- fsusieR:::rename_format_output(
    obj,
    names_colX=paste0("snp", 1:5),
    tidx=c(2L, 4L),
    kept_index=c(1L, 3L, 5L),
    original_P=5L
  )

  expect_equal(restored$L, 2L)
  expect_equal(restored$P, 5L)
  expect_equal(restored$fitted_P, 3L)
  expect_equal(restored$variable_index, c(1L, 3L, 5L))
  expect_equal(restored$removed_variable_index, c(2L, 4L))
  expect_equal(
    restored$cs[[1]],
    structure(c(1L, 5L), names=c("snp1", "snp5"))
  )
  expect_equal(
    restored$cs[[2]],
    structure(3L, names="snp3")
  )
  expect_true(all(vapply(restored$alpha, length, integer(1)) == 5L))
  expect_true(all(vapply(restored$lBF, length, integer(1)) == 5L))
  expect_true(all(vapply(restored$fitted_wc, nrow, integer(1)) == 5L))
  expect_true(all(vapply(restored$fitted_wc2, nrow, integer(1)) == 5L))
  expect_equal(unname(restored$pip[c(2, 4)]), c(0, 0))
  expect_equal(restored$n_cs, restored$L)
  expect_equal(restored$cs_size, c(2L, 1L))
  expect_identical(fsusieR:::rename_format_output(restored), restored)
})


test_that("fSuSiE restores unnamed input columns by numeric position", {
  obj <- make_dynamic_susif(L=1L, P=2L)
  obj$cs <- list(c(1L, 2L))

  restored <- fsusieR:::rename_format_output(
    obj,
    names_colX=NULL,
    tidx=c(1L, 3L),
    kept_index=c(2L, 4L),
    original_P=4L
  )

  expect_null(names(restored$pip))
  expect_equal(restored$cs[[1]], c(2L, 4L))
  expect_equal(unname(restored$pip[c(1, 3)]), c(0, 0))
})


test_that("fSuSiE rejects an all-constant design before fitting", {
  expect_error(
    susiF(
      Y=matrix(rnorm(32), nrow=8),
      X=matrix(1, nrow=8, ncol=2),
      verbose=FALSE
    ),
    "All columns of X are constant"
  )
})

test_that("fSuSiE singleton credible sets pass the purity criterion", {
  obj <- make_dynamic_susif(L=2L, P=2L)
  obj$cs <- list(1L, 2L)
  obj$est_pi <- lapply(seq_len(obj$L), function(l) c(0.5, 0.5))

  expect_length(
    which_dummy_cs(obj, X=matrix(c(0, 1, 0, 1), nrow=2)),
    0L
  )
})
