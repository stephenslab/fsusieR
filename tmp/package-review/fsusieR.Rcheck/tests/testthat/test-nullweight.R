test_that("mixture updates match ASH independently of the SNP cutoff", {
  # Identical SNP profiles carry the same evidence at every cutoff. Their
  # total weight should remain one, giving the same fit as ASH on one profile.
  b <- c(3, 0, -2, 0.2, 1, 4, -4, 0.1)
  Bhat <- matrix(rep(b, each = 40), nrow = 40)
  Shat <- matrix(1, nrow = 40, ncol = length(b))
  indx <- gen_wavelet_indx(3)
  mixsd <- c(0, 1, 3)
  g <- list(fitted_g = ashr::normalmix(c(0.8, 0.1, 0.1), 0, mixsd))
  control <- list(verbose = FALSE, eps = 1e-10, numiter.em = 40)

  for (prior in c("mixture_normal", "mixture_normal_per_scale")) {
    groups <- if (prior == "mixture_normal") list(seq_along(b)) else indx
    G_prior <- structure(rep(list(g), length(groups)), class = prior)

    for (nullweight in c(1, 10)) {
      expected <- unlist(lapply(groups, function(j) {
        ashr::ash(b[j], rep(1, length(j)), g = g$fitted_g,
                  mixcompdist = "normal", nullweight = nullweight,
                  control = control, outputlevel = 0)$fitted_g$pi
      }))

      for (cutoff in c(1, 10, 20, 40)) {
        fit <- EM_pi(G_prior, Bhat, Shat, indx,
                     control_mixsqp = control, lowc_wc = NULL,
                     nullweight = nullweight, max_SNP_EM = cutoff,
                     tol_null_prior = 0)
        expect_equal(as.numeric(unlist(fit$tpi_k)), as.numeric(expected),
                     tolerance = 1e-6,
                     info = paste(prior, nullweight, cutoff))
      }
    }
  }
})

test_that("susiF rejects nullweight values outside the ASH convention", {
  for (value in list(0, 0.1, -1, Inf, NA_real_, numeric(0), c(1, 10), "10")) {
    expect_error(susiF(nullweight = value),
                 "nullweight.*finite number at least 1")
  }
})

test_that("EBmvFR rejects nullweight values outside the ASH convention", {
  for (value in list(0, 0.1, -1, Inf, NA_real_, numeric(0), c(1, 10), "10")) {
    expect_error(EBmvFR(nullweight = value),
                 "nullweight.*finite number at least 1")
  }
})
