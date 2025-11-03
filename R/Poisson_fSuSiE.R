#' @export
Pois_fSuSiE <- function(Y,
                        Z = NULL,
                        X = NULL,
                        L = 3,
                        scaling = NULL,
                        ebps_method = c("ebpm",
                                        'pois_mean_split',
                                        'ind_pois_mean_split',
                                        'ind_ebps',
                                        'ind_poisson_smoothing',
                                        'nugget'),
                        reflect = FALSE,
                        verbose = TRUE,
                        tol = 1e-3,
                        maxit_outer = 10,
                        maxit_inner = 100,
                        control_mixsqp = list(verbose = FALSE,
                                              eps = 1e-6,
                                              numiter.em = 4),
                        nullweight = 10,
                        cov_lev = 0.95,
                        min_purity = 0.5,
                        cor_small = TRUE,
                        post_processing = "HMM",
                            print=TRUE,
                        update_Mu_each_iter = TRUE) {

  # Validate inputs
  if (is.null(X) && is.null(Z)) {
    stop("Please provide at least one of X or Z matrix")
  }

  # Determine fit approach
  has_X <- !is.null(X)
  has_Z <- !is.null(Z)

  if (!has_X) {
    fit_approach <- "penalized"
    if (verbose) print("No X provided: performing penalized regression only")
  } else if (!has_Z) {
    fit_approach <- "fine_mapping"
    if (verbose) print("No Z provided: performing fine-mapping only")
  } else {
    fit_approach <- "both"
    if (verbose) print("Both X and Z provided: performing joint fine-mapping and penalized regression")
  }

  # Match arguments
  ebps_method <- match.arg(ebps_method)

  # Handle reflection for non-dyadic length
  J <- log2(ncol(Y))
  if ((J %% 1) != 0) reflect <- TRUE

  if (reflect) {
    tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i, ]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx
  } else {
    idx_out <- 1:ncol(Y)
  }

  # Remove constant columns
  if (has_X) {
    tidx <- which(apply(X, 2, var) == 0)
    if (length(tidx) > 0) {
      warning(paste("Removed", length(tidx), "constant columns from X"))
      X <- X[, -tidx, drop = FALSE]
    }

    est_effect_fm <- matrix(0, nrow =ncol(X), ncol = ncol(Y))      # X effects
  }else{

    est_effect_fm <-0     # X effects
  }

  if (has_Z) {
    tidx <- which(apply(Z, 2, var) == 0)
    if (length(tidx) > 0) {
      warning(paste("Removed", length(tidx), "constant columns from Z"))
      Z <- Z[, -tidx, drop = FALSE]
    }
  }

  # Setup scaling
  if (is.null(scaling)) {
    scaling <- rep(1, nrow(Y))
  } else {
    if (length(scaling) != nrow(Y)) {
      stop("scaling must have length equal to nrow(Y)")
    }
  }

  # ============================================================================
  # STEP 1 (INITIALIZATION): Compute initial Mu_pm via empirical Bayes
  # ============================================================================
  if (verbose) cat("Initializing Mu_pm via", ebps_method, "...\n")

  Mu_pm_init <- log(Y + 1)

  if (ebps_method == "ebpm") {
    tt <- ebpm_normal(c(Y), s = rep(scaling, ncol(Y)))
    Mu_pm <- matrix(tt$posterior$mean_log, byrow = FALSE, ncol = ncol(Y))

  } else if (ebps_method == "pois_mean_split") {
    tt <- pois_mean_split(c(Y), s = rep(scaling, ncol(Y)),
                          mu_pm_init = c(Mu_pm_init))
    Mu_pm <- matrix(tt$posterior$mean_log, byrow = FALSE, ncol = ncol(Y))

  } else if (ebps_method == "ind_pois_mean_split") {
    Mu_pm <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
    for (i in 1:nrow(Y)) {
      Mu_pm[i, ] <- pois_mean_split(Y[i, ], s = scaling[i],
                                    mu_pm_init = log(Y[i, ] + 1))$posterior$mean_log
    }

  } else if (ebps_method == "ind_ebps") {
    Mu_pm <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
    for (i in 1:nrow(Y)) {
      Mu_pm[i, ] <- ebps(Y[i, ], s = scaling[i])$posterior$mean_log
    }

  } else if (ebps_method == "ind_poisson_smoothing") {
    Mu_pm <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
    for (i in 1:nrow(Y)) {
      Mu_pm[i, ] <- log(pois_smooth_split(Y[i, ], s = scaling[i],
                                          Eb_init = log(Y[i, ] + 1))$Emean)
    }

  } else if (ebps_method == "nugget") {
    Mu_pm <- fit_latent_nugget(Y)$Y
  }

  # Initialize components
  alpha_0 <- rowMeans(Mu_pm)  # Intercept per sample
  Theta_pm <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))  # Z effects

  sigma2 <-  .01  # Initial variance

  print(dim(est_effect_fm))
  # ============================================================================
  # COORDINATE ASCENT ALGORITHM (Algorithm 5 from supplement)
  # ============================================================================
  converged <- FALSE
  iter <- 1
  elbo_hist <- numeric(maxit_outer)

  while (!converged && iter <= maxit_outer) {

    if (verbose) cat("\n=== Outer iteration", iter, "===\n")

    # Store old values for convergence
    alpha_0_old <- alpha_0
    Theta_pm_old <- Theta_pm
    B_pm_old <- est_effect_fm

    # ------------------------------------------------------------------------
    # STEP 1: Update Mu_pm | current estimates (optional, expensive)
    # ------------------------------------------------------------------------
    if (update_Mu_each_iter && iter > 1) {
      if (verbose) cat("Step 1: Updating Mu_pm...\n")

      # Re-solve Poisson-Normal mean problem with current mean structure
      b_it <- matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) + Theta_pm + B_pm

      # This would require re-running split-VA step
      # For efficiency, paper suggests skipping this for fine-mapping
      # Mu_pm <- update_Mu_pm_split_VA(Y, b_it, sigma2, scaling)
      tt <-    pois_mean_GP(x=c(Y),
                            prior_mean = c(B_pm),
                            s =  rep( scaling, ncol(Y)),
                            prior_var = sigma2 )
      Mu_pm <- matrix( tt$posterior$posteriorMean_latent,byrow = FALSE, ncol=ncol(Y))
      Mu_pv <- matrix( tt$posterior$posteriorVar_latent ,byrow = FALSE, ncol=ncol(Y))
    }

    # Transform to wavelet space (if needed by downstream functions)
    # A = Mu_pm %*% W (wavelet transform)

    # ------------------------------------------------------------------------
    # STEP 2a: Update Theta (Z effects) via EBmvFR
    # ------------------------------------------------------------------------
    if (has_Z) {
      if (verbose) cat("Step 2a: Updating Z effects (Theta)...\n")

      # Residualize: A_{-B} = A - A_0 - XB
      resid_for_Z <- Mu_pm - matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) - B_pm

      EBmvFR.obj <- EBmvFR(
        Y = resid_for_Z,
        X = Z,
        tol = tol,
        control_mixsqp = control_mixsqp,
        nullweight = nullweight,
        cal_obj = FALSE,
        verbose = FALSE,
        maxit = maxit_inner,
        max_step_EM = 1
      )

      Theta_pm <- Z %*% EBmvFR.obj$fitted_func

      if (verbose) {
        cat("  Z effects: max =", round(max(abs(Theta_pm)), 4), "\n")
      }
    } else {
      EBmvFR.obj <- NULL
    }

    # ------------------------------------------------------------------------
    # STEP 2b: Update B (X effects) via fSuSiE
    # ------------------------------------------------------------------------
    if (has_X) {
      if (verbose) cat("Step 2b: Updating X effects (B) via fSuSiE...\n")

      # Residualize: A_{-Theta} = A - A_0 - Z*Theta
      resid_for_X <- Mu_pm - matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) - Theta_pm

      susiF.obj <- susiF(
        Y = resid_for_X,
        X = X,
        L = L,
        tol = tol,
        control_mixsqp = control_mixsqp,
        nullweight = nullweight,
        cal_obj = FALSE,
        verbose = FALSE,
        cov_lev = cov_lev,
        min_purity = min_purity,
        maxit = maxit_inner,
        cor_small = cor_small,
        post_processing = post_processing
      )

      # Reconstruct fitted values
      if (length(susiF.obj$cs) > 0) {
        est_effect_fm=Reduce("+", lapply(1:length(susiF.obj$cs), function(l) {
          t(susiF.obj$fitted_func[[l]] %*% t(susiF.obj$alpha[[l]]))
        }))
        print(dim(est_effect_fm))

        print(dim(B_pm_old))
        B_pm <- X %*% est_effect_fm
      } else {
        B_pm <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
      }

      if (verbose) {
        cat("  Found", length(susiF.obj$cs), "credible sets\n")
        cat("  X effects: max =", round(max(abs(B_pm)), 4), "\n")
      }
    } else {
      susiF.obj <- NULL
    }

    # ------------------------------------------------------------------------
    # STEP 2c: Update intercept alpha_0 via EBNM
    # ------------------------------------------------------------------------
    if (verbose) cat("Step 2c: Updating intercept (alpha_0)...\n")

    # Residualize: A_{-(Theta,B)} = A - Z*Theta - X*B
    resid_for_alpha0 <- Mu_pm - Theta_pm - B_pm
    alpha_0 <- rowMeans(resid_for_alpha0)  # Simple mean across positions

    # Could use EBNM for more sophisticated shrinkage
    # alpha_0 <- ebnm(rowMeans(resid_for_alpha0), s = sqrt(sigma2 / ncol(Y)))$posterior$mean

    # ------------------------------------------------------------------------
    # STEP 3: Update variance sigma2
    # ------------------------------------------------------------------------
    if (verbose) cat("Step 3: Updating sigma2...\n")

    residuals <- Mu_pm - matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) - Theta_pm - B_pm
    sigma2= var(c(residuals))
    if (verbose) cat("  sigma2 =", round(sigma2, 6), "\n")

    # ------------------------------------------------------------------------
    # Check convergence
    # ------------------------------------------------------------------------
    diff_alpha <- max(abs(alpha_0 - alpha_0_old))
    diff_Theta <- max(abs(Theta_pm - Theta_pm_old))
    diff_B <- max(abs(est_effect_fm - B_pm_old))
    max_diff <- max(diff_alpha, diff_Theta, diff_B)

    if (verbose) {
      cat("Convergence: max change =", round(max_diff, 8), "\n")
      cat("  |Δalpha_0| =", round(diff_alpha, 8), "\n")
      cat("  |ΔTheta| =", round(diff_Theta, 8), "\n")
      cat("  |ΔB| =", round(diff_B, 8), "\n")
    }
    if (print){
par (mfrow=c(3,1))
      plot( log1p(Y ),  (Mu_pm  ))
      plot(susiF.obj$fitted_func[[1]])
      plot(susiF.obj$fitted_func[[2]])

      abline(a=0,b=1)
      par (mfrow=c(1,1))
    }

    converged <- max_diff < tol
    iter <- iter + 1
  }

  if (!converged && verbose) {
    warning("Algorithm did not converge in ", maxit_outer, " iterations")
  }

  # ============================================================================
  # Finalize and return results
  # ============================================================================
  if (has_X && !is.null(susiF.obj)) {
    susiF.obj <- fsusieR::update_cal_pip(susiF.obj)
  }

  # Compute fitted intensities
  fitted_log_intensity <- matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) + Theta_pm + B_pm
  fitted_intensity <- exp(fitted_log_intensity)

  # Prepare output
  out <- list(
    Mu_pm = Mu_pm,
    alpha_0 = alpha_0,
    Theta_pm = Theta_pm,
    B_pm = B_pm,
    sigma2 = sigma2,
    fitted = fitted_intensity[, idx_out],
    converged = converged,
    n_iter = iter - 1
  )

  if (has_X) out$susiF.obj <- susiF.obj
  if (has_Z) out$EBmvFR.obj <- EBmvFR.obj

  return(out)
}
