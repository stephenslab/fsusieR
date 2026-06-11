#' @title Poisson functional SuSiE
#'
#' @description Fine-mapping of inhomogeneous Poisson process data via the
#'   split-variational empirical-Bayes framework of Denault, Xie and
#'   Stephens. Counts \eqn{y_{i,t}\sim\mathrm{Pois}(s_i\exp(\mu_{i,t}))}
#'   are modelled with a Gaussian latent log-intensity
#'   \eqn{\mu_{i,t}\sim N(\alpha_{0,t}+b_{i,t},\sigma^{2})}, where the
#'   per-position intercept \eqn{\alpha_{0,t}} and the SNP effects
#'   \eqn{b_{i,t}=\sum_j x_{i,j}F_{j,t}} are estimated by an inner fSuSiE
#'   fit to the variational posterior mean of \eqn{M}.
#'
#' @param Y N x T integer matrix of counts.
#' @param X N x J covariate matrix. Required.
#' @param L Upper bound on the number of single effects.
#' @param scaling Length-N vector of per-sample scaling factors \eqn{s_i}.
#'   Defaults to ones.
#' @param reflect Logical; if `TRUE` and `ncol(Y)` is not a power of two,
#'   the data are mirrored to the next power of two before wavelet
#'   processing. Set automatically when needed.
#' @param maxit_outer Maximum number of outer (coordinate-ascent) iterations.
#' @param maxit_inner Maximum number of inner IBSS iterations forwarded to
#'   [susiF()].
#' @param tol Convergence tolerance on the per-iteration change in `B_pm`
#'   (used as a backup criterion to the ELBO check).
#' @param elbo_tol Relative tolerance on the ELBO increment used for the
#'   primary convergence check (default `1e-4`).
#' @param control_mixsqp,nullweight,cov_lev,min_purity,cor_small,post_processing,thresh_lowcount
#'   Forwarded to [susiF()].
#' @param verbose Logical; print outer-iteration progress.
#' @param diagnostic_plot Logical; if `TRUE`, plot fit diagnostics each
#'   outer iteration. Replaces the old `print` argument.
#' @param True_intensity Optional N x T matrix of ground-truth log
#'   intensities, used only when `diagnostic_plot = TRUE`. Not used in
#'   any computation.
#' @param print Deprecated alias for `diagnostic_plot`.
#' @param Z Reserved for future support of unpenalised technical
#'   covariates. Must be `NULL` (errors otherwise).
#' @param update_Mu_each_iter Logical; if `FALSE`, the latent-mean update
#'   is skipped after the first outer iteration. Mainly for debugging.
#' @param s2 Initial value of \eqn{\sigma^{2}}. Default `0.5`.
#'
#' @return A list with components `Mu_pm`, `Mu_pv`, `B_pm`, `B_pv`,
#'   `est_effect_fm`, `sigma2`, `fitted_latent` (lead-SNP reconstructed
#'   latent log-intensity), `fitted` (`= exp(fitted_latent)`), `elbo_trace`,
#'   `converged`, `n_iter`, `susiF.obj`. See Details.
#'
#' @details
#' The algorithm is a two-block coordinate ascent on the split-variational
#' ELBO
#' \deqn{F = \sum_{i,t}\bigl[y_{i,t}\bar\mu_{i,t}-s_i e^{\bar\mu_{i,t}+\nu_{i,t}/2}\bigr]
#'   -\tfrac{NT}{2}\log\sigma^{2}
#'   -\tfrac{1}{2\sigma^{2}}\!\sum_{i,t}\!\bigl[(\bar\mu_{i,t}-\bar b_{i,t})^{2}+\nu_{i,t}+\mathrm{Var}_q[b_{i,t}]\bigr]
#'   +\tfrac{1}{2}\sum_{i,t}\log\nu_{i,t}
#'   -\sum_r\!\mathrm{KL}(q_r\|g_r) + \text{const.}}
#' At each outer iteration we
#' \enumerate{
#'   \item update the per-sample VGA posterior on the latent log-rate
#'     via the Newton solver [vga_pois_solver()];
#'   \item refit [susiF()] on `Mu_pm` to get an updated PIP-weighted SNP
#'     effect matrix and rebuild `B_pm = colMeans(Mu_pm) +
#'     (X - colMeans(X)) %*% est_effect_fm`. This re-parametrises the
#'     per-position intercept \eqn{\alpha_{0,t}} as the column mean of the
#'     susiF response, and uses the PIP-weighted prediction so that all
#'     SNPs in a CS contribute (not just the leading SNP). This drives the
#'     optimisation; the *reported* `fitted` value instead uses the lead-SNP
#'     post-processed reconstruction (`build_fitted_leadSNP`), consistent
#'     with susiF's `ind_fitted_func`;
#'   \item compute the variance contribution `Var_q[b_{i,t}]` (see below)
#'     and update \eqn{\sigma^{2}} from its closed-form ELBO maximiser.
#' }
#' `Var_q[b_{i,t}]` captures the SuSiE selection uncertainty about which
#' SNP within an SER is causal:
#' \deqn{\mathrm{Var}_q[b_{i,t}] \approx
#'   \sum_r (m^{(r)}_t)^{2}\,\bigl(A^{(r)}_i - (M^{(r)}_i)^{2}\bigr),
#'   \quad A^{(r)}_i=\sum_j x_{i,j}^{2}\alpha^{(r)}_j,\;M^{(r)}_i=\sum_j x_{i,j}\alpha^{(r)}_j,}
#' which is zero when the CS is one-hot (high purity) and grows with
#' selection uncertainty. The effect-magnitude posterior variance is
#' absorbed into `sigma2` for simplicity.
#'
#' The per-SER KL terms are taken from `susiF.obj$KL` when available.
#' Note that those KLs are computed with susiF's own internal
#' \eqn{\sigma^{2}_{\text{inner}}}; they enter the ELBO as a soft
#' regulariser whose absolute scale is approximate. The σ²-update and
#' VGA blocks of the ELBO are exact under the model assumptions.
#'
#' Convergence: we declare convergence on the relative ELBO change,
#' `|ELBO_new - ELBO_old| / |ELBO_old| < elbo_tol`. The per-iteration
#' `max(abs(B_pm - B_pm_old))` is still reported but is **only a
#' diagnostic** -- it is not used to gate convergence, because within-CS
#' SNP label-switching keeps it above `tol` indefinitely and previously
#' forced every fit to run all `maxit_outer` iterations.
#'
#' Differences from the manuscript: the per-position intercept
#' \eqn{\alpha_{0,t}} is implicitly handled by susiF's own column
#' centring of its response (rather than via an EBNM prior, as the paper
#' describes in Eq. 35). The technical-covariate matrix `Z` is not yet
#' supported.
#'
#' @seealso [susiF()], [vga_pois_solver()]
#'
#' @export
Pois_fSuSiE <- function(Y,
                        X,
                        L = 3,
                        scaling = NULL,
                        reflect = FALSE,
                        maxit_outer = 10,
                        maxit_inner = 10,
                        tol = 1e-3,
                        elbo_tol = 1e-4,
                        control_mixsqp = list(verbose = FALSE,
                                              eps = 1e-6,
                                              numiter.em = 4),
                        nullweight = 0.10,
                        cov_lev = 0.95,
                        min_purity = 0.5,
                        cor_small = FALSE,
                        post_processing = "smash",
                        thresh_lowcount = 0,
                        verbose = TRUE,
                        diagnostic_plot = FALSE,
                        True_intensity = NULL,
                        print = NULL,
                        Z = NULL,
                        update_Mu_each_iter = TRUE,
                        s2 = 0.5) {

  ## ------------------------------------------------------------------ ##
  ## 1. Input validation                                                ##
  ## ------------------------------------------------------------------ ##
  if (missing(X) || is.null(X)) stop("Please provide X matrix")
  if (missing(Y) || is.null(Y)) stop("Please provide Y matrix")
  if (!is.null(Z))
    stop("Z covariate handling is not yet implemented; please pass Z = NULL.")

  if (!is.null(print)) diagnostic_plot <- isTRUE(print)

  tidx <- which(apply(X, 2, var) == 0)
  if (length(tidx) > 0) {
    warning(sprintf("Removed %d constant columns from X", length(tidx)))
    X <- X[, -tidx, drop = FALSE]
  }

  Jlog <- log2(ncol(Y))
  if ((Jlog %% 1) != 0) reflect <- TRUE
  if (reflect) {
    tl <- lapply(seq_len(nrow(Y)), function(i) reflect_vec(Y[i, ]))
    Y <- do.call(rbind, lapply(seq_along(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx
  } else {
    idx_out <- seq_len(ncol(Y))
  }

  N  <- nrow(Y)
  Tt <- ncol(Y)
  if (is.null(scaling)) {
    scaling <- rep(1, N)
  } else if (length(scaling) != N) {
    stop("`scaling` must have length equal to nrow(Y)")
  }

  ## ------------------------------------------------------------------ ##
  ## 2. Initialisation                                                  ##
  ##    Warm-start B_pm from fSuSiE on log1p(Y/s) and rebuild the       ##
  ##    intercept ourselves via build_B_pm() (see helper comment for    ##
  ##    why susiF.obj$ind_fitted_func cannot be used directly).         ##
  ## ------------------------------------------------------------------ ##
  Y_warm <- log1p(sweep(Y, 1, scaling, "/"))
  susiF.obj <- susiF(Y_warm, X = X, L = L, verbose = FALSE)

  est_effect_fm <- matrix(0, nrow = ncol(X), ncol = Tt)
  if (length(susiF.obj$cs) > 0) {
    est_effect_fm <- reconstruct_effect(susiF.obj)
    B_pm <- build_B_pm(Y_warm, X, est_effect_fm)
  } else {
    B_pm <- matrix(colMeans(Y_warm), nrow = N, ncol = Tt, byrow = TRUE)
  }
  B_pv <- compute_B_pv(susiF.obj, X)            # N x T, SNP-selection var

  Mu_pm  <- B_pm
  Mu_pv  <- matrix(1 / Tt, nrow = N, ncol = Tt)
  sigma2 <- max(s2, 1e-6)

  ## ------------------------------------------------------------------ ##
  ## 3. Outer coordinate ascent                                         ##
  ## ------------------------------------------------------------------ ##
  converged   <- FALSE
  iter        <- 1
  elbo_trace  <- numeric(0)
  elbo_prev   <- -Inf

  while (!converged && iter <= maxit_outer) {

    if (verbose) message("=== Pois_fSuSiE outer iter ", iter, " ===")
    B_pm_old <- B_pm

    ## 3a. Per-row VGA: q(mu_i) ~ N(Mu_pm[i,], diag(Mu_pv[i,])) given
    ##     prior mean B_pm[i,], variance sigma2.
    if (update_Mu_each_iter || iter == 1) {
      for (i in seq_len(N)) {
        sol <- try(vga_pois_solver(init_val = Mu_pm[i, ],
                                   x        = Y[i, ],
                                   s        = scaling[i],
                                   beta     = B_pm[i, ],
                                   sigma2   = sigma2),
                   silent = TRUE)
        if (inherits(sol, "try-error") ||
            !all(is.finite(sol$m)) || !all(is.finite(sol$v)) ||
            any(sol$v <= 0)) {
          ## Solver failed or returned a non-finite / invalid estimate.
          ## Fall back to the log1p transform so downstream susiF/ash never
          ## sees NaN/Inf (which previously crashed ashr::squarem with
          ## "missing value where TRUE/FALSE needed").
          Mu_pm[i, ] <- log1p(Y[i, ] / scaling[i])
          Mu_pv[i, ] <- sigma2
        } else {
          Mu_pm[i, ] <- sol$m
          Mu_pv[i, ] <- sol$v
        }
      }
      ## Final safety net: clamp the latent log-rate to a sane range and
      ## guarantee strictly-positive, finite variances before refitting.
      Mu_pm[!is.finite(Mu_pm)] <- 0
      Mu_pm <- pmin(pmax(Mu_pm, -30), 30)
      Mu_pv[!is.finite(Mu_pv) | Mu_pv <= 0] <- sigma2
    }

    ## 3b. Refit fSuSiE on Mu_pm. susiF column-centres its response;
    ##     we rebuild B_pm via build_B_pm() so the intercept is restored.
    susiF.obj <- susiF(Y               = Mu_pm,
                       X               = X,
                       L               = L,
                       tol             = tol,
                       maxit           = maxit_inner,
                       control_mixsqp  = control_mixsqp,
                       nullweight      = nullweight,
                       cov_lev         = cov_lev,
                       min_purity      = min_purity,
                       cor_small       = cor_small,
                       post_processing = post_processing,
                       thresh_lowcount = thresh_lowcount,
                       cal_obj         = FALSE,
                       verbose         = FALSE)

    if (length(susiF.obj$cs) > 0) {
      est_effect_fm <- reconstruct_effect(susiF.obj)
      B_pm <- build_B_pm(Mu_pm, X, est_effect_fm)
    } else {
      est_effect_fm[] <- 0
      B_pm <- matrix(colMeans(Mu_pm), nrow = N, ncol = Tt, byrow = TRUE)
    }
    B_pv <- compute_B_pv(susiF.obj, X)

    ## 3c. Closed-form sigma^2 maximiser of the joint ELBO:
    ##       sigma^2 = mean( (Mu_pm - B_pm)^2 + Mu_pv + B_pv ).
    sigma2 <- max(mean((Mu_pm - B_pm)^2 + Mu_pv + B_pv), 1e-6)

    ## 3d. Joint ELBO (drops Poisson y!-constants; sums with the susiF
    ##     per-SER KL when available).
    elbo <- pois_fsusie_elbo(Y, scaling, Mu_pm, Mu_pv,
                             B_pm, B_pv, sigma2, susiF.obj)
    elbo_trace <- c(elbo_trace, elbo)

    if (verbose)
      message(sprintf("  sigma2 = %.5f  ELBO = %.4f", sigma2, elbo))

    if (diagnostic_plot)
      pois_fsusie_diagnostic_plot(Y, Mu_pm, susiF.obj, True_intensity, B_pm)

    ## 3e. Convergence: relative ELBO change small (the principled VI
    ##     stopping rule). max|dB_pm| is reported only as a DIAGNOSTIC and
    ##     no longer gates convergence: within a credible set the selected
    ##     SNP can keep switching between near-equivalent (highly correlated)
    ##     SNPs, so the PIP-weighted B_pm wiggles at the 1e-2 level
    ##     indefinitely without changing the ELBO. Gating on max|dB| < tol
    ##     therefore almost never triggered and forced every fit to run the
    ##     full maxit_outer iterations (each a costly susiF refit), which was
    ##     the main reason benchmark jobs hit their wall-clock limit.
    if (iter > 1) {
      d_elbo <- elbo - elbo_prev
      rel    <- abs(d_elbo) / (abs(elbo_prev) + 1e-12)
      d_B    <- max(abs(B_pm - B_pm_old))
      if (verbose)
        message(sprintf("  dELBO = %.4g  rel = %.3g  max|dB| = %.3g (diag)",
                        d_elbo, rel, d_B))
      converged <- (rel < elbo_tol)
    }
    elbo_prev <- elbo
    iter <- iter + 1
  }

  if (!converged && verbose)
    warning("Pois_fSuSiE did not converge in ", maxit_outer, " iterations")

  ## ------------------------------------------------------------------ ##
  ## 4. Finalise output                                                 ##
  ## ------------------------------------------------------------------ ##
  susiF.obj <- fsusieR::update_cal_pip(susiF.obj)

  ## Fitted value reported to the user: lead-SNP reconstruction from the
  ## POST-PROCESSED effect curves, with the Y/Mu offset added back (see
  ## build_fitted_leadSNP). The internal prior-mean B_pm above (PIP-weighted)
  ## is left unchanged and still drives the VGA / sigma2 / ELBO updates.
  fitted_latent    <- build_fitted_leadSNP(Mu_pm, X, susiF.obj)
  fitted_intensity <- exp(fitted_latent)

  list(Mu_pm         = Mu_pm,
       Mu_pv         = Mu_pv,
       B_pm          = B_pm,
       B_pv          = B_pv,
       est_effect_fm = est_effect_fm,
       sigma2        = sigma2,
       fitted_latent = fitted_latent[, idx_out, drop = FALSE],
       fitted        = fitted_intensity[, idx_out, drop = FALSE],
       elbo_trace    = elbo_trace,
       converged     = converged,
       n_iter        = iter - 1,
       susiF.obj     = susiF.obj)
}


## ======================================================================
## Internal helpers
## ======================================================================

# PIP-weighted reconstruction of the J x T SNP-effect matrix from a
# fitted susiF object. Iterates over credible sets so that pruned /
# low-purity SERs do not contribute.
reconstruct_effect <- function(susiF.obj) {
  J  <- length(susiF.obj$alpha[[1]])
  Tt <- length(susiF.obj$fitted_func[[1]])
  out <- matrix(0, nrow = J, ncol = Tt)
  if (length(susiF.obj$cs) == 0) return(out)
  for (l in seq_along(susiF.obj$cs)) {
    out <- out + outer(susiF.obj$alpha[[l]],
                       as.numeric(susiF.obj$fitted_func[[l]]))
  }
  out
}



# Build B_pm = alpha_0 + (X - colMeans(X)) %*% est_effect_fm
# where alpha_0_t = colMeans(response)_t is the per-position intercept.
# This is the correct prior-mean for the latent log-rate: it reparametrises
# alpha_0 + X %*% F as colMeans(response) + X_centered %*% F so that the
# intercept matches the column means of the susiF response.
#
# This PIP-weighted prior-mean drives the optimisation. (The reported
# fitted value uses the lead-SNP post-processed reconstruction instead --
# see build_fitted_leadSNP() / susiF.obj$ind_fitted_func, which is now
# populated by out_prep.susiF.)
build_B_pm <- function(response, X, est_effect_fm) {
  X_centered <- sweep(X, 2, colMeans(X), "-")
  effect_pred <- X_centered %*% est_effect_fm
  intercept_t <- colMeans(response)
  sweep(effect_pred, 2, intercept_t, "+")
}

# Lead-SNP reconstruction of the per-individual latent log-intensity, using
# susiF's POST-PROCESSED effect curves (susiF.obj$fitted_func) rather than the
# raw PIP-weighted estimate. This is the analogue of susiF's ind_fitted_func
# and is what we report as the fitted value (output only; the optimisation
# still uses the PIP-weighted build_B_pm above).
#
# Offset handling mirrors build_B_pm / update_cal_indf.susiF: the post-
# processed curves are slopes fit on column-centred data, so we add back the
# per-position mean of the response (colMeans) and centre X by its column
# means before applying each lead-SNP effect:
#   yhat_i(t) = mean_t(response) + sum_l ( X[i,j_l] - mean(X[,j_l]) ) * f_l(t)
# with j_l = which.max(alpha_l).
build_fitted_leadSNP <- function(response, X, susiF.obj) {
  N  <- nrow(X)
  Tt <- ncol(response)
  out <- matrix(colMeans(response), nrow = N, ncol = Tt, byrow = TRUE)
  if (length(susiF.obj$cs) == 0 || length(susiF.obj$fitted_func) == 0)
    return(out)
  X_centered <- sweep(X, 2, colMeans(X), "-")
  for (l in seq_along(susiF.obj$fitted_func)) {
    j   <- which.max(susiF.obj$alpha[[l]])
    f_l <- as.numeric(susiF.obj$fitted_func[[l]])
    out <- out + outer(X_centered[, j], f_l)
  }
  out
}

# Approximate Var_q[b_{i,t}] under the leading-SNP-curve approximation,
# keeping the SNP-selection variance term and dropping the smaller
# effect-magnitude variance term (absorbed into sigma^2). Returns an
# N x T matrix; zero in the high-purity limit.
compute_B_pv <- function(susiF.obj, X) {
  N <- nrow(X)
  if (length(susiF.obj$cs) == 0 ||
      length(susiF.obj$fitted_func) == 0) {
    return(matrix(0, nrow = N, ncol = 1))
  }
  Tt  <- length(susiF.obj$fitted_func[[1]])
  Bpv <- matrix(0, nrow = N, ncol = Tt)
  X2  <- X * X
  for (l in seq_along(susiF.obj$cs)) {
    alpha_l <- as.numeric(susiF.obj$alpha[[l]])
    m_l     <- as.numeric(susiF.obj$fitted_func[[l]])
    A_l     <- as.numeric(X2 %*% alpha_l)        # N-vector
    M_l     <- as.numeric(X  %*% alpha_l)        # N-vector
    sel_var <- pmax(A_l - M_l^2, 0)               # SNP-selection variance
    Bpv     <- Bpv + outer(sel_var, m_l^2)
  }
  pmax(Bpv, 0)
}

# Joint variational ELBO (up to constants in Y).
# Components:
#   Poisson: sum(y*mu - s*exp(mu + v/2))           (drop log(y!), y log s)
#   Gaussian fit + entropy of q(mu):
#       -(NT/2) log sigma2
#       - (1/(2 sigma2)) sum((mu - b)^2 + nu + Var_b)
#       + 0.5 sum log nu                            (drop (NT/2)(1+log(2 pi)))
#   - sum(susiF.obj$KL)  (if available; computed with susiF's inner sigma2)
pois_fsusie_elbo <- function(Y, scaling, Mu_pm, Mu_pv,
                             B_pm, B_pv, sigma2,
                             susiF.obj = NULL) {
  N  <- nrow(Y)
  Tt <- ncol(Y)
  NT <- N * Tt
  s_mat <- matrix(scaling, nrow = N, ncol = Tt)

  pois_term <- sum(Y * Mu_pm) - sum(s_mat * exp(Mu_pm + Mu_pv / 2))

  resid_sq  <- sum((Mu_pm - B_pm)^2 + Mu_pv + B_pv)
  gauss_fit <- -(NT / 2) * log(sigma2) - resid_sq / (2 * sigma2)
  entropy   <- 0.5 * sum(log(pmax(Mu_pv, 1e-300)))

  kl_F <- 0
  if (!is.null(susiF.obj) && !is.null(susiF.obj$KL)) {
    kl_F <- sum(susiF.obj$KL, na.rm = TRUE)
  }

  pois_term + gauss_fit + entropy - kl_F
}


# Diagnostic plot when diagnostic_plot = TRUE.
pois_fsusie_diagnostic_plot <- function(Y, Mu_pm, susiF.obj,
                                        True_intensity, B_pm) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (is.null(True_intensity)) {
    graphics::par(mfrow = c(2, 1))
    plot(log1p(Y), Mu_pm, xlab = "log1p(Y)", ylab = "Mu_pm",
         main = "Latent log-intensity recovery")
    graphics::abline(a = 0, b = 1)
    if (length(susiF.obj$fitted_func) >= 1)
      plot(susiF.obj$fitted_func[[1]], type = "l",
           main = "fSuSiE fitted_func[[1]]")
  } else {
    graphics::par(mfrow = c(2, 2))
    plot(log1p(Y), Mu_pm, main = "Mu_pm vs log1p(Y)"); graphics::abline(0, 1)
    trmse <- mean((c(B_pm) - c(True_intensity))^2)
    plot(True_intensity, B_pm,
         main = sprintf("B_pm vs truth  (MSE %.4f)", trmse))
    graphics::abline(0, 1)
    if (length(susiF.obj$fitted_func) >= 1)
      plot(susiF.obj$fitted_func[[1]], type = "l",
           main = "fSuSiE fitted_func[[1]]")
    if (length(susiF.obj$fitted_func) >= 2)
      plot(susiF.obj$fitted_func[[2]], type = "l",
           main = "fSuSiE fitted_func[[2]]")
  }
}
