#' @export
Pois_fSuSiE2 <- function(Y,
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
                        post_processing = "smash",
                        print=TRUE,
                        update_Mu_each_iter = TRUE,
                        True_intensity=NULL,
                        s2=1) {

  # Validate inputs
  if (is.null(X) && is.null(Z)) {
    stop("Please provide at least one of X or Z matrix")
  }

  has_X <- !is.null(X)
  has_Z <- !is.null(Z)

  ebps_method <- match.arg(ebps_method)

  # Handle reflection for non-dyadic length
  n_bins <- ncol(Y)
  N <- nrow(Y)
  J <- log2(n_bins)
  if ((J %% 1) != 0) reflect <- TRUE

  if (reflect) {
    tl <- lapply(1:N, function(i) reflect_vec(Y[i, ]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx
  } else {
    idx_out <- 1:n_bins
  }

  # Setup scaling
  if (is.null(scaling)) {
    scaling <- rep(1, N)
  }

  # ============================================================================
  # STEP 1 (INITIALIZATION)
  # ============================================================================
  if (verbose) cat("Initializing Mu_pm via", ebps_method, "...\n")

  # Use your chosen initialization to get Mu_pm and Mu_pv (posterior mean/var of log-lambda)
  # This part remains mostly as you had it, but ensure we keep variance
  Mu_pv <- matrix(1/n_bins, nrow = N, ncol = ncol(Y))

  if (ebps_method == "ebpm") {
    tt <- ebpm_normal(c(Y), s = rep(scaling, ncol(Y)))
    Mu_pm <- matrix(tt$posterior$mean_log, byrow = FALSE, ncol = ncol(Y))
    Mu_pv <- matrix(tt$posterior$var_log, byrow = FALSE, ncol = ncol(Y))
  } else {
    # Defaulting to a safe VGA init for other methods
    Mu_pm <- log(Y + 0.1)
  }

  # Initialize components
  alpha_0 <- rowMeans(Mu_pm)
  Theta_pm <- matrix(0, nrow = N, ncol = ncol(Y))
  B_pm <- matrix(0, nrow = N, ncol = ncol(Y))
  est_effect_fm <- if(has_X) matrix(0, nrow = ncol(X), ncol = ncol(Y)) else 0
  sigma2 <- 0.1

  # ============================================================================
  # COORDINATE ASCENT
  # ============================================================================
  converged <- FALSE
  iter <- 1

  while (!converged && iter <= maxit_outer) {

    if (verbose) cat("\n=== Outer iteration", iter, "===\n")

    alpha_0_old <- alpha_0
    Theta_pm_old <- Theta_pm
    B_pm_old <- B_pm # Fixed: now comparing fitted values to fitted values

    # ------------------------------------------------------------------------
    # STEP 1: Update Mu_pm (VGA Step)
    # ------------------------------------------------------------------------
    if (update_Mu_each_iter && iter > 1) {
      if (verbose) cat("Step 1: Updating latent Mu_pm (VGA)...\n")
      # current mean structure: b_it = alpha_0 + Theta + B
      current_mean <- matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) + Theta_pm + B_pm

      for (i in 1:N) {
        # Fast Newton/Bisection solver from your vga_pois_solver
        res_vga <- vga_pois_solver(init_val = Mu_pm[i,],
                                   x = Y[i,],
                                   s = scaling[i],
                                   beta = current_mean[i,],
                                   sigma2 = sigma2)
        Mu_pm[i,] <- res_vga$m
        Mu_pv[i,] <- res_vga$v
      }
    }

    # ------------------------------------------------------------------------
    # STEP 2a: Update Theta (Z effects)
    # ------------------------------------------------------------------------
    if (has_Z) {
      resid_for_Z <- Mu_pm - matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) - B_pm
      EBmvFR.obj <- EBmvFR(Y = resid_for_Z, X = Z, maxit = maxit_inner, verbose = FALSE)
      Theta_pm <- Z %*% EBmvFR.obj$fitted_func
    }

    # ------------------------------------------------------------------------
    # STEP 2b: Update B (X effects) via fSuSiE
    # ------------------------------------------------------------------------
    if (has_X) {
      resid_for_X <- Mu_pm - matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) - Theta_pm

      susiF.obj <- susiF(
        Y = resid_for_X,
        X = X,
        L = L,
        verbose = FALSE,
        post_processing = post_processing,
      )

      if (length(susiF.obj$cs) > 0) {
        # Correct reconstruction of effects
        est_effect_fm <- Reduce("+", lapply(1:length(susiF.obj$alpha), function(l) {
          susiF.obj$alpha[[l]] %*% t(susiF.obj$fitted_func[[l]])
        }))

        B_pm <- X %*%  (est_effect_fm)
      } else {
        B_pm <- matrix(0, nrow = N, ncol = ncol(Y))
      }
    }

    # ------------------------------------------------------------------------
    # STEP 2c: Update intercept alpha_0 (The FIX)
    # ------------------------------------------------------------------------
    resid_for_alpha0 <- Mu_pm - Theta_pm - B_pm
    alpha_0 <- rowMeans(resid_for_alpha0) # Removed the 0*

    # ------------------------------------------------------------------------
    # STEP 3: Update variance sigma2 (The FIX)
    # ------------------------------------------------------------------------
    # sigma2 is the nugget variance in the Split-VA framework
    sigma2 <- mean((Mu_pm - B_pm - Theta_pm - matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)))^2 + Mu_pv)
    if (verbose) cat("  sigma2 =", round(sigma2, 6), "\n")

    # ------------------------------------------------------------------------
    # Convergence Check
    # ------------------------------------------------------------------------
    max_diff <- max(
      max(abs(alpha_0 - alpha_0_old)),
      max(abs(Theta_pm - Theta_pm_old)),
      max(abs(B_pm - B_pm_old))
    )

    if (verbose) cat("Convergence: max change =", round(max_diff, 8), "\n")

    converged <- max_diff < tol
    iter <- iter + 1
  }

  # Finalize
  if (has_X && !is.null(susiF.obj)) {
    susiF.obj <- fsusieR::update_cal_pip(susiF.obj)
  }

  fitted_log_intensity <- matrix(rep(alpha_0, ncol(Y)), ncol = ncol(Y)) + Theta_pm + B_pm

  return(list(
    Mu_pm = Mu_pm,
    alpha_0 = alpha_0,
    Theta_pm = Theta_pm,
    B_pm = B_pm,
    sigma2 = sigma2,
    fitted = exp(fitted_log_intensity)[, idx_out],
    susiF.obj = if(has_X) susiF.obj else NULL,
    converged = converged
  ))
}
