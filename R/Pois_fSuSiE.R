#' @title Poisson functional SuSiE
#'
#' @description Fine-mapping of inhomogeneous Poisson process data with a
#'   Gaussian variational approximation for the latent log-intensity and an
#'   inner [susiF()] fit for the structured genetic mean.
#'
#' @param Y N by T numeric matrix of non-negative integer counts.
#' @param X N by J numeric covariate matrix.
#' @param L Upper bound on the number of single effects.
#' @param scaling Length-N vector of finite, strictly positive sample-specific
#'   Poisson scaling factors. Defaults to one.
#' @param reflect Logical; mirror the observations before wavelet processing.
#'   Reflection is enabled automatically when `ncol(Y)` is not a power of two.
#' @param maxit_outer Maximum number of outer iterations.
#' @param maxit_inner Maximum number of inner [susiF()] iterations.
#' @param tol Positive convergence tolerance forwarded to [susiF()].
#' @param elbo_tol Deprecated alias for `outer_tol`. Kept for compatibility.
#' @param control_mixsqp,nullweight,cov_lev,min_purity,cor_small,post_processing,thresh_lowcount
#'   Forwarded to [susiF()]. `post_processing = "none"` is not supported by
#'   this approximation because its `fitted_func` is already alpha-collapsed.
#' @param filter.number,family Wavelet filter passed to [susiF()].
#' @param verbose Logical; report outer-iteration progress.
#' @param diagnostic_plot Logical; draw diagnostic plots each outer iteration.
#' @param True_intensity Optional N by T matrix of true latent log-intensities,
#'   used only by `diagnostic_plot`.
#' @param print Deprecated alias for `diagnostic_plot`.
#' @param Z Reserved for future dense-covariate support. Must be `NULL`.
#' @param update_Mu_each_iter Logical; update the latent Gaussian posterior at
#'   every outer iteration. Setting this to `FALSE` is intended for debugging.
#' @param s2 Positive fallback initial value of the noise variance `sigma2`.
#' @param warm_start_sigma2 Logical; if `TRUE` (the default), use the existing
#'   `maxit_inner`-iteration `susiF(log1p(Y/scaling))` warm fit to initialize
#'   `sigma2`. The estimate is computed on the original log1p scale, subtracts
#'   first-order Poisson sampling variance, and corrects log1p attenuation.
#' @param outer_tol Positive convergence tolerance for the maximum relative
#'   change in `Mu_pm`, `B_pm`, and `sigma2`.
#' @param stable_iterations Number of consecutive iterations below `outer_tol`
#'   required for convergence.
#' @param solver_failure Either `"error"` (the default) or `"log1p"`. The latter
#'   records a failed VGA row and uses a log1p fallback instead of stopping.
#'
#' @return A list containing posterior moments, fitted values, `sigma2`, an
#'   iteration-level `convergence_trace`, solver failures, and the inner
#'   `susiF.obj`. `est_effect_fm` is the raw variational posterior mean effect
#'   matrix in the original X units; `postprocessed_effect_fm` is retained for
#'   reporting. `sigma2_used_for_final_susif` records the variance used to fit
#'   the final inner posterior, before its subsequent M-step. The
#'   `partial_objective_trace` contains only the Poisson latent and Gaussian
#'   coupling terms. The legacy `elbo_trace` component is retained as an alias,
#'   but is not a complete joint ELBO.
#'
#' @details
#' The noise-variance update is
#' \deqn{\sigma^2 = \mathrm{mean}\{(\bar\mu-\bar B)^2 + V_\mu + V_B\}.}
#' Here `B` denotes the complete structured Gaussian mean, including the
#' plug-in intercept. `V_B` is computed from the raw fSuSiE variational
#' factors, including both SNP-selection uncertainty and conditional
#' wavelet-coefficient uncertainty.
#'
#' During each Gaussian block update, [susiF()] centers but does not scale the
#' latent response and holds its residual variance fixed at the current outer
#' `sigma2`. Post-processing is used only for reported effect curves; it does
#' not feed the variational moments or the variance update. The
#' `partial_objective_trace` still excludes the fSuSiE prior/KL contribution,
#' so convergence is based on parameter changes rather than this diagnostic.
#'
#' @seealso [susiF()], [vga_pois_solver()]
#'
#' @export
Pois_fSuSiE <- function(Y,
                        X,
                        L = 3,
                        scaling = NULL,
                        reflect = FALSE,
                        maxit_outer = 20,
                        maxit_inner = 10,
                        tol = 1e-3,
                        elbo_tol = NULL,
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
                        s2 = 0.5,
                        warm_start_sigma2 = TRUE,
                        outer_tol = 1e-3,
                        stable_iterations = 2L,
                        solver_failure = c("error", "log1p"),
                        filter.number = 10,
                        family = "DaubLeAsymm") {

  if (missing(X) || is.null(X)) stop("Please provide X matrix", call. = FALSE)
  if (missing(Y) || is.null(Y)) stop("Please provide Y matrix", call. = FALSE)
  if (!is.null(Z)) {
    stop("Z covariate handling is not implemented; please pass Z = NULL.",
         call. = FALSE)
  }
  if (!is.null(print)) diagnostic_plot <- isTRUE(print)
  if (!is.null(elbo_tol)) {
    warning("`elbo_tol` is deprecated; using it as `outer_tol`.",
            call. = FALSE)
    outer_tol <- elbo_tol
  }
  post_processing <- match.arg(post_processing, c("smash", "TI", "HMM"))
  solver_failure <- match.arg(solver_failure)

  checked <- validate_pois_fsusie_inputs(Y, X, L, scaling)
  Y <- checked$Y
  X <- checked$X
  L <- checked$L
  scaling <- checked$scaling

  positive_scalar <- function(x) {
    length(x) == 1L && is.finite(x) && x > 0
  }
  positive_integer <- function(x) {
    positive_scalar(x) && x == as.integer(x)
  }
  if (!is.logical(warm_start_sigma2) || length(warm_start_sigma2) != 1L ||
      is.na(warm_start_sigma2)) {
    stop("`warm_start_sigma2` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!positive_integer(maxit_outer) || !positive_integer(maxit_inner) ||
      !positive_integer(stable_iterations) || !positive_scalar(tol) ||
      !positive_scalar(outer_tol) || !positive_scalar(s2)) {
    stop("Iteration counts, tolerances, and `s2` must be strictly positive.",
         call. = FALSE)
  }
  maxit_outer <- as.integer(maxit_outer)
  maxit_inner <- as.integer(maxit_inner)
  stable_iterations <- as.integer(stable_iterations)

  original_grid_size <- ncol(Y)
  true_intensity_internal <- NULL
  if (!is.null(True_intensity)) {
    true_intensity_internal <- as.matrix(True_intensity)
    if (!is.numeric(true_intensity_internal) ||
        !identical(dim(true_intensity_internal), dim(Y)) ||
        any(!is.finite(true_intensity_internal))) {
      stop("`True_intensity` must be a finite numeric matrix with dim(Y).",
           call. = FALSE)
    }
  }

  must_reflect <- (log2(original_grid_size) %% 1) != 0
  if (isTRUE(reflect) || must_reflect) {
    reflected <- lapply(seq_len(nrow(Y)), function(i) reflect_vec(Y[i, ]))
    Y <- do.call(rbind, lapply(reflected, `[[`, "x"))
    output_index <- reflected[[1L]]$idx
    if (!is.null(true_intensity_internal)) {
      reflected_truth <- lapply(seq_len(nrow(true_intensity_internal)),
                                function(i) {
                                  reflect_vec(true_intensity_internal[i, ])$x
                                })
      true_intensity_internal <- do.call(rbind, reflected_truth)
    }
  } else {
    output_index <- seq_len(original_grid_size)
  }

  N <- nrow(Y)
  Tt <- ncol(Y)
  inverse_dwt <- pois_inverse_dwt_matrix(
    n_grid = Tt,
    filter.number = filter.number,
    family = family
  )
  Y_warm <- log1p(sweep(Y, 1L, scaling, "/"))
  susiF.obj <- fit_pois_susif(
    response = Y_warm,
    X = X,
    L = L,
    tol = tol,
    maxit_inner = maxit_inner,
    control_mixsqp = control_mixsqp,
    nullweight = nullweight,
    cov_lev = cov_lev,
    min_purity = min_purity,
    cor_small = cor_small,
    post_processing = post_processing,
    thresh_lowcount = thresh_lowcount,
    filter.number = filter.number,
    family = family,
    standardize_Y = TRUE,
    residual_variance = NULL
  )

  postprocessed_effect_fm <- reconstruct_effect(susiF.obj, ncol(X), Tt)
  est_effect_fm <- postprocessed_effect_fm
  B_pm <- build_B_pm(Y_warm, X, postprocessed_effect_fm)
  B_pv <- compute_B_pv(susiF.obj, X, Tt)
  Mu_pm <- B_pm
  Mu_pv <- matrix(1 / Tt, nrow = N, ncol = Tt)
  warm_start <- estimate_pois_sigma2_warm_start(
    response_log1p = Y_warm,
    fitted_log1p = B_pm,
    fitted_variance = B_pv,
    scaling = scaling,
    fallback_sigma2 = s2,
    enabled = warm_start_sigma2
  )
  sigma2 <- warm_start$sigma2
  if (verbose) {
    message(sprintf("Initial sigma2 = %.5g (%s)", sigma2,
                    warm_start$method))
  }

  partial_objective_trace <- numeric(0)
  convergence_trace <- vector("list", maxit_outer)
  solver_failures <- data.frame(iteration = integer(0),
                                row = integer(0),
                                message = character(0))
  stable_count <- 0L
  converged <- FALSE

  for (iter in seq_len(maxit_outer)) {
    if (verbose) message("=== Pois_fSuSiE outer iter ", iter, " ===")
    old_Mu_pm <- Mu_pm
    old_B_pm <- B_pm
    old_sigma2 <- sigma2

    if (isTRUE(update_Mu_each_iter) || iter == 1L) {
      for (i in seq_len(N)) {
        failure_message <- NULL
        sol <- tryCatch(
          vga_pois_solver(init_val = Mu_pm[i, ],
                          x = Y[i, ],
                          s = scaling[i],
                          beta = B_pm[i, ],
                          sigma2 = sigma2),
          error = function(e) {
            failure_message <<- conditionMessage(e)
            NULL
          }
        )
        invalid <- !vga_pois_solution_is_valid(
          sol,
          x = Y[i, ],
          s = scaling[i],
          beta = B_pm[i, ],
          sigma2 = sigma2,
          tol = 1e-4
        )

        if (invalid) {
          if (is.null(failure_message)) {
            failure_message <- paste(
              "invalid VGA output: moments do not satisfy the Poisson",
              "variational score equations"
            )
          }
          solver_failures <- rbind(
            solver_failures,
            data.frame(iteration = iter, row = i, message = failure_message)
          )
          if (solver_failure == "error") {
            stop(sprintf("VGA failed at iteration %d, row %d: %s",
                         iter, i, failure_message), call. = FALSE)
          }
          Mu_pm[i, ] <- log1p(Y[i, ] / scaling[i])
          Mu_pv[i, ] <- sigma2
        } else {
          Mu_pm[i, ] <- sol$m
          Mu_pv[i, ] <- sol$v
        }
      }
    }

    sigma2_for_inner <- sigma2
    susiF.obj <- fit_pois_susif(
      response = Mu_pm,
      X = X,
      L = L,
      tol = tol,
      maxit_inner = maxit_inner,
      control_mixsqp = control_mixsqp,
      nullweight = nullweight,
      cov_lev = cov_lev,
      min_purity = min_purity,
      cor_small = cor_small,
      post_processing = post_processing,
      thresh_lowcount = thresh_lowcount,
      filter.number = filter.number,
      family = family,
      standardize_Y = FALSE,
      residual_variance = sigma2_for_inner
    )
    if (!isTRUE(all.equal(susiF.obj$sigma2, sigma2_for_inner,
                          tolerance = 0))) {
      stop("The inner susiF fit changed its fixed residual variance.",
           call. = FALSE)
    }
    postprocessed_effect_fm <- reconstruct_effect(susiF.obj, ncol(X), Tt)
    structured_moments <- compute_pois_structured_moments(
      susiF.obj = susiF.obj,
      X = X,
      response = Mu_pm,
      inverse_dwt = inverse_dwt
    )
    B_pm <- structured_moments$mean
    B_pv <- structured_moments$variance
    est_effect_fm <- structured_moments$coefficient_mean

    sigma2_update <- pois_sigma2_update(Mu_pm, Mu_pv, B_pm, B_pv)
    sigma2 <- sigma2_update$sigma2
    residual_component <- sigma2_update$residual_component
    latent_variance_component <- sigma2_update$latent_variance_component
    structured_variance_component <-
      sigma2_update$structured_variance_component

    partial_objective <- pois_fsusie_partial_objective(
      Y, scaling, Mu_pm, Mu_pv, B_pm, B_pv, sigma2
    )
    partial_objective_trace <- c(partial_objective_trace, partial_objective)

    delta_Mu <- relative_rms_change(Mu_pm, old_Mu_pm)
    delta_B <- relative_rms_change(B_pm, old_B_pm)
    delta_sigma2 <- abs(sigma2 - old_sigma2) / max(old_sigma2, 1e-8)
    maximum_change <- max(delta_Mu, delta_B, delta_sigma2)
    convergence_trace[[iter]] <- data.frame(
      iteration = iter,
      sigma2 = sigma2,
      inner_sigma2 = sigma2_for_inner,
      noise_sd = sqrt(sigma2),
      residual_component = residual_component,
      latent_variance_component = latent_variance_component,
      structured_variance_component = structured_variance_component,
      relative_rms_Mu = delta_Mu,
      relative_rms_B = delta_B,
      relative_sigma2 = delta_sigma2,
      maximum_change = maximum_change,
      partial_objective = partial_objective
    )

    stable_count <- if (maximum_change < outer_tol) stable_count + 1L else 0L
    if (verbose) {
      message(sprintf(
        "  sigma2=%.5g  max relative change=%.3g  partial objective=%.5g",
        sigma2, maximum_change, partial_objective
      ))
    }
    if (diagnostic_plot) {
      pois_fsusie_diagnostic_plot(
        Y, Mu_pm, susiF.obj, true_intensity_internal, B_pm
      )
    }
    if (stable_count >= stable_iterations) {
      converged <- TRUE
      break
    }
  }

  n_iter <- length(partial_objective_trace)
  convergence_trace <- do.call(rbind, convergence_trace[seq_len(n_iter)])
  if (!converged && verbose) {
    warning("Pois_fSuSiE did not converge in ", maxit_outer, " iterations",
            call. = FALSE)
  }
  susiF.obj <- update_cal_pip(susiF.obj)

  posterior_count_mean <- sweep(
    exp(pmin(Mu_pm + Mu_pv / 2, 700)), 1L, scaling, "*"
  )
  crop <- function(value) value[, output_index, drop = FALSE]
  fitted_latent <- crop(B_pm)

  list(
    Mu_pm = crop(Mu_pm),
    Mu_pv = crop(Mu_pv),
    B_pm = crop(B_pm),
    B_pv = crop(B_pv),
    est_effect_fm = est_effect_fm[, output_index, drop = FALSE],
    postprocessed_effect_fm =
      postprocessed_effect_fm[, output_index, drop = FALSE],
    sigma2 = sigma2,
    sigma2_used_for_final_susif = susiF.obj$sigma2,
    sigma2_initial = warm_start$sigma2,
    sigma2_warm_start = warm_start,
    fitted_latent = fitted_latent,
    fitted = exp(pmin(fitted_latent, 700)),
    posterior_count_mean = crop(posterior_count_mean),
    partial_objective_trace = partial_objective_trace,
    elbo_trace = partial_objective_trace,
    elbo_trace_is_partial = TRUE,
    convergence_trace = convergence_trace,
    converged = converged,
    n_iter = n_iter,
    solver_failures = solver_failures,
    susiF.obj = susiF.obj,
    approximation = list(
      exact_joint_cavi = FALSE,
      reason = paste(
        "The Gaussian fSuSiE block and sigma2 M-step now match the stated",
        "mean-field updates; the intercept remains a plug-in column mean and",
        "each inner block is optimized for at most maxit_inner iterations."
      ),
      B_moments = paste(
        "B_pm and B_pv use the raw alpha, conditional wavelet means, and",
        "conditional wavelet variances; post-processing is output-only."
      ),
      objective = paste(
        "partial_objective_trace excludes the fSuSiE prior/KL contribution and",
        "is a diagnostic, not the complete joint ELBO."
      )
    ),
    reflected_internal_grid = Tt != original_grid_size,
    output_grid_index = output_index
  )
}


## Internal helpers -------------------------------------------------------

validate_pois_fsusie_inputs <- function(Y, X, L, scaling) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  if (!is.numeric(Y) || !is.numeric(X)) {
    stop("`Y` and `X` must be numeric matrices.", call. = FALSE)
  }
  if (nrow(Y) < 2L || ncol(Y) < 2L) {
    stop("`Y` must have at least two rows and two columns.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y) || ncol(X) < 1L) {
    stop("`X` must have nrow(X) == nrow(Y) and at least one column.",
         call. = FALSE)
  }
  if (any(!is.finite(Y)) || any(!is.finite(X))) {
    stop("`Y` and `X` must contain only finite values.", call. = FALSE)
  }
  if (any(Y < 0) || any(abs(Y - round(Y)) > 1e-8)) {
    stop("`Y` must contain non-negative integer Poisson counts.",
         call. = FALSE)
  }
  X_var <- apply(X, 2L, stats::var)
  if (any(!is.finite(X_var)) || any(X_var <= .Machine$double.eps)) {
    stop("`X` contains a constant column; remove it before fitting.",
         call. = FALSE)
  }
  if (length(L) != 1L || !is.finite(L) || L != as.integer(L) ||
      L < 1L || L > ncol(X)) {
    stop("`L` must be an integer between 1 and ncol(X).", call. = FALSE)
  }
  if (is.null(scaling)) scaling <- rep(1, nrow(Y))
  if (!is.numeric(scaling) || length(scaling) != nrow(Y) ||
      any(!is.finite(scaling)) || any(scaling <= 0)) {
    stop("`scaling` must contain one finite, strictly positive value per row.",
         call. = FALSE)
  }
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  list(Y = Y, X = X, L = as.integer(L), scaling = as.numeric(scaling))
}


## Initialize latent log-rate variance from the already-computed log1p fit.
## The internal susiF sigma2 is deliberately not used: it is estimated after
## position/wavelet standardization and is not on the latent log-rate scale.
estimate_pois_sigma2_warm_start <- function(response_log1p,
                                            fitted_log1p,
                                            fitted_variance,
                                            scaling,
                                            fallback_sigma2,
                                            enabled = TRUE) {
  fallback_sigma2 <- max(as.numeric(fallback_sigma2), 1e-8)
  fallback <- list(
    sigma2 = fallback_sigma2,
    method = "user-supplied s2 fallback",
    used = FALSE,
    transformed_residual_mse = NA_real_,
    structured_variance = NA_real_,
    poisson_variance = NA_real_,
    attenuation_sq = NA_real_,
    untruncated_sigma2 = NA_real_
  )
  if (!isTRUE(enabled)) {
    fallback$method <- "user-supplied s2"
    return(fallback)
  }

  valid <- is.matrix(response_log1p) && is.matrix(fitted_log1p) &&
    is.matrix(fitted_variance) &&
    identical(dim(response_log1p), dim(fitted_log1p)) &&
    identical(dim(response_log1p), dim(fitted_variance)) &&
    length(scaling) == nrow(response_log1p) &&
    all(is.finite(response_log1p)) && all(is.finite(fitted_log1p)) &&
    all(is.finite(fitted_variance)) && all(fitted_variance >= 0) &&
    all(is.finite(scaling)) && all(scaling > 0)
  if (!valid) return(fallback)

  fitted_rate <- pmax(expm1(pmin(fitted_log1p, 700)), 0)
  scaling_matrix <- matrix(scaling, nrow(response_log1p),
                           ncol(response_log1p))
  transformed_residual_mse <- mean((response_log1p - fitted_log1p)^2)
  structured_variance <- mean(fitted_variance)
  poisson_variance <- mean(
    fitted_rate / (scaling_matrix * (1 + fitted_rate)^2)
  )
  attenuation_sq <- mean((fitted_rate / (1 + fitted_rate))^2)
  if (!is.finite(attenuation_sq) ||
      attenuation_sq <= sqrt(.Machine$double.eps)) return(fallback)

  untruncated_sigma2 <-
    (transformed_residual_mse + structured_variance - poisson_variance) /
    attenuation_sq
  if (!is.finite(untruncated_sigma2)) return(fallback)

  list(
    sigma2 = max(untruncated_sigma2, 1e-8),
    method = "log1p-fSuSiE delta correction",
    used = TRUE,
    transformed_residual_mse = transformed_residual_mse,
    structured_variance = structured_variance,
    poisson_variance = poisson_variance,
    attenuation_sq = attenuation_sq,
    untruncated_sigma2 = untruncated_sigma2
  )
}


fit_pois_susif <- function(response,
                           X,
                           L,
                           tol,
                           maxit_inner,
                           control_mixsqp,
                           nullweight,
                           cov_lev,
                           min_purity,
                           cor_small,
                           post_processing,
                           thresh_lowcount,
                           filter.number = 10,
                           family = "DaubLeAsymm",
                           standardize_Y = TRUE,
                           residual_variance = NULL) {
  fit <- susiF(
    Y = response,
    X = X,
    L = L,
    L_start = min(3L, L),
    tol = tol,
    maxit = maxit_inner,
    control_mixsqp = control_mixsqp,
    nullweight = nullweight,
    cov_lev = cov_lev,
    min_purity = min_purity,
    cor_small = cor_small,
    post_processing = post_processing,
    thresh_lowcount = thresh_lowcount,
    filter.number = filter.number,
    family = family,
    standardize_Y = standardize_Y,
    residual_variance = residual_variance,
    filter_cs = FALSE,
    cal_obj = FALSE,
    verbose = FALSE
  )
  if (is.null(fit)) {
    stop("susiF returned NULL (almost all wavelet coefficients were filtered).",
         call. = FALSE)
  }
  fit
}


## Matrix mapping wavelet coefficients in the package's c(D, C) ordering
## back to the observation grid. A row vector of coefficients c is inverted by
## c %*% inverse_dwt.
pois_inverse_dwt_matrix <- function(n_grid,
                                    filter.number = 10,
                                    family = "DaubLeAsymm") {
  if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 2L ||
      n_grid != as.integer(n_grid) || log2(n_grid) %% 1 != 0) {
    stop("`n_grid` must be an integer power of two greater than one.",
         call. = FALSE)
  }
  n_grid <- as.integer(n_grid)
  ## GenW stores coefficients as c(C, D), whereas fSuSiE stores c(D, C).
  ## Reorder the rows of the inverse transform to match fSuSiE exactly.
  inverse_cd <- t(wavethresh::GenW(
    n = n_grid,
    filter.number = filter.number,
    family = family
  ))
  inverse_dwt <- inverse_cd[c(seq.int(2L, n_grid), 1L), , drop = FALSE]
  if (!identical(dim(inverse_dwt), c(n_grid, n_grid)) ||
      any(!is.finite(inverse_dwt))) {
    stop("Could not construct the inverse wavelet-transform matrix.",
         call. = FALSE)
  }
  inverse_dwt
}


pois_raw_effects <- function(susiF.obj) {
  effect_counts <- c(
    length(susiF.obj$alpha),
    length(susiF.obj$fitted_wc),
    length(susiF.obj$fitted_wc2)
  )
  if (length(unique(effect_counts)) != 1L ||
      (!is.null(susiF.obj$L) && susiF.obj$L != effect_counts[1L])) {
    stop("The raw susiF effect fields are not aligned.", call. = FALSE)
  }
  n_effect <- effect_counts[1L]
  if (n_effect == 0L) integer(0) else seq_len(n_effect)
}


## First two pointwise moments of the complete structured mean under the raw
## fSuSiE variational distribution. For each single effect, gamma selects one
## SNP for the entire curve, so the calculation retains the covariance across
## wavelet coefficients induced by that shared SNP selection.
compute_pois_structured_moments <- function(susiF.obj,
                                            X,
                                            response,
                                            inverse_dwt) {
  X <- as.matrix(X)
  response <- as.matrix(response)
  inverse_dwt <- as.matrix(inverse_dwt)
  N <- nrow(X)
  P <- ncol(X)
  Tt <- ncol(response)
  if (nrow(response) != N ||
      !identical(dim(inverse_dwt), c(Tt, Tt))) {
    stop("Incompatible dimensions in the structured-moment calculation.",
         call. = FALSE)
  }

  x_scale <- as.numeric(susiF.obj$csd_X)
  if (length(x_scale) != P || any(!is.finite(x_scale)) ||
      any(x_scale <= 0)) {
    x_scale <- apply(X, 2L, stats::sd)
  }
  if (any(!is.finite(x_scale)) || any(x_scale <= 0)) {
    stop("Cannot reconstruct moments for constant columns of X.",
         call. = FALSE)
  }
  X_centered <- sweep(X, 2L, colMeans(X), "-")
  X_scaled <- sweep(X_centered, 2L, x_scale, "/")

  structured_mean <- matrix(0, nrow = N, ncol = Tt)
  structured_variance <- matrix(0, nrow = N, ncol = Tt)
  coefficient_mean <- matrix(0, nrow = P, ncol = Tt)
  inverse_dwt_sq <- inverse_dwt^2

  for (l in pois_raw_effects(susiF.obj)) {
    alpha_l <- as.numeric(susiF.obj$alpha[[l]])
    mean_wc <- as.matrix(susiF.obj$fitted_wc[[l]])
    variance_wc <- as.matrix(susiF.obj$fitted_wc2[[l]])
    if (length(alpha_l) != P ||
        !identical(dim(mean_wc), c(P, Tt)) ||
        !identical(dim(variance_wc), c(P, Tt)) ||
        any(!is.finite(alpha_l)) || any(alpha_l < 0) ||
        any(!is.finite(mean_wc)) || any(!is.finite(variance_wc))) {
      stop("Invalid raw posterior moments in the susiF object.",
           call. = FALSE)
    }
    alpha_mass <- sum(alpha_l)
    if (alpha_mass <= 1e-12) {
      if (any(abs(mean_wc) > 1e-12)) {
        stop("An uninitialized susiF effect has non-zero posterior means.",
             call. = FALSE)
      }
      next
    }
    if (abs(alpha_mass - 1) > 1e-6) {
      stop("A fitted susiF effect has alpha values that do not sum to one.",
           call. = FALSE)
    }

    conditional_mean <- mean_wc %*% inverse_dwt
    conditional_variance <- pmax(variance_wc, 0) %*% inverse_dwt_sq
    weighted_mean <- sweep(conditional_mean, 1L, alpha_l, "*")
    weighted_second <- sweep(
      conditional_mean^2 + conditional_variance,
      1L,
      alpha_l,
      "*"
    )
    effect_mean <- X_scaled %*% weighted_mean
    effect_second <- (X_scaled^2) %*% weighted_second

    structured_mean <- structured_mean + effect_mean
    structured_variance <- structured_variance +
      pmax(effect_second - effect_mean^2, 0)
    coefficient_mean <- coefficient_mean +
      sweep(weighted_mean, 1L, x_scale, "/")
  }

  structured_mean <- sweep(
    structured_mean,
    2L,
    colMeans(response),
    "+"
  )
  list(
    mean = structured_mean,
    variance = structured_variance,
    coefficient_mean = coefficient_mean
  )
}


pois_sigma2_update <- function(Mu_pm, Mu_pv, B_pm, B_pv,
                               minimum_sigma2 = 1e-8) {
  Mu_pm <- as.matrix(Mu_pm)
  Mu_pv <- as.matrix(Mu_pv)
  B_pm <- as.matrix(B_pm)
  B_pv <- as.matrix(B_pv)
  if (length(minimum_sigma2) != 1L || !is.finite(minimum_sigma2) ||
      minimum_sigma2 <= 0) {
    stop("`minimum_sigma2` must be a finite positive scalar.", call. = FALSE)
  }
  if (!identical(dim(Mu_pm), dim(Mu_pv)) ||
      !identical(dim(Mu_pm), dim(B_pm)) ||
      !identical(dim(Mu_pm), dim(B_pv)) ||
      any(!is.finite(Mu_pm)) || any(!is.finite(Mu_pv)) ||
      any(!is.finite(B_pm)) || any(!is.finite(B_pv)) ||
      any(Mu_pv < 0) || any(B_pv < 0)) {
    stop("Invalid posterior moments in the sigma2 update.", call. = FALSE)
  }
  residual_component <- mean((Mu_pm - B_pm)^2)
  latent_variance_component <- mean(Mu_pv)
  structured_variance_component <- mean(B_pv)
  untruncated_sigma2 <- residual_component + latent_variance_component +
    structured_variance_component
  if (!is.finite(untruncated_sigma2) || untruncated_sigma2 <= 0) {
    stop("The sigma2 update produced a non-positive or non-finite value.",
         call. = FALSE)
  }
  list(
    sigma2 = max(untruncated_sigma2, minimum_sigma2),
    untruncated_sigma2 = untruncated_sigma2,
    residual_component = residual_component,
    latent_variance_component = latent_variance_component,
    structured_variance_component = structured_variance_component
  )
}


active_pois_effects <- function(susiF.obj) {
  if (is.null(susiF.obj$cs) || length(susiF.obj$cs) == 0L ||
      is.null(susiF.obj$alpha) || is.null(susiF.obj$fitted_func)) {
    return(integer(0))
  }
  n_effect <- min(length(susiF.obj$alpha), length(susiF.obj$fitted_func))
  if (n_effect == 0L) integer(0) else seq_len(n_effect)
}


## PIP-weighted reconstruction under the post-processed lead-curve
## approximation. This is not an exact raw posterior-moment reconstruction.
reconstruct_effect <- function(susiF.obj, n_snp, n_grid) {
  out <- matrix(0, nrow = n_snp, ncol = n_grid)
  for (l in active_pois_effects(susiF.obj)) {
    alpha_l <- as.numeric(susiF.obj$alpha[[l]])
    mean_l <- as.numeric(susiF.obj$fitted_func[[l]])
    if (length(alpha_l) != n_snp || length(mean_l) != n_grid) {
      stop("Incompatible alpha/fitted_func dimensions in the susiF object.",
           call. = FALSE)
    }
    out <- out + tcrossprod(alpha_l, mean_l)
  }
  out
}


build_B_pm <- function(response, X, est_effect_fm) {
  X_centered <- sweep(X, 2L, colMeans(X), "-")
  effect_prediction <- X_centered %*% est_effect_fm
  sweep(effect_prediction, 2L, colMeans(response), "+")
}


## Approximate Var[B_it] under the same alpha-distributed lead-curve model.
## Uses centered X and the available post-processed curve variance:
## Var = E[X_j^2] (m_t^2 + v_t) - E[X_j]^2 m_t^2.
compute_B_pv <- function(susiF.obj, X, n_grid) {
  X <- as.matrix(X)
  N <- nrow(X)
  P <- ncol(X)
  B_pv <- matrix(0, nrow = N, ncol = n_grid)
  X_centered <- sweep(X, 2L, colMeans(X), "-")

  for (l in active_pois_effects(susiF.obj)) {
    alpha_l <- as.numeric(susiF.obj$alpha[[l]])
    mean_l <- as.numeric(susiF.obj$fitted_func[[l]])
    if (length(alpha_l) != P || length(mean_l) != n_grid ||
        any(!is.finite(alpha_l)) || any(alpha_l < 0)) {
      stop("Invalid alpha/fitted_func values in the susiF object.",
           call. = FALSE)
    }

    curve_variance <- rep(0, n_grid)
    if (!is.null(susiF.obj$fitted_var) &&
        length(susiF.obj$fitted_var) >= l &&
        !is.null(susiF.obj$fitted_var[[l]])) {
      candidate <- as.numeric(susiF.obj$fitted_var[[l]])
      if (length(candidate) == n_grid && all(is.finite(candidate))) {
        curve_variance <- pmax(candidate, 0)
      }
    }

    second_X <- as.numeric((X_centered^2) %*% alpha_l)
    mean_X <- as.numeric(X_centered %*% alpha_l)
    effect_second_moment <- mean_l^2 + curve_variance
    B_pv <- B_pv + tcrossprod(second_X, effect_second_moment) -
      tcrossprod(mean_X^2, mean_l^2)
  }
  pmax(B_pv, 0)
}


pois_fsusie_partial_objective <- function(Y,
                                          scaling,
                                          Mu_pm,
                                          Mu_pv,
                                          B_pm,
                                          B_pv,
                                          sigma2) {
  if (!is.finite(sigma2) || sigma2 <= 0 || any(Mu_pv <= 0)) return(-Inf)
  N <- nrow(Y)
  Tt <- ncol(Y)
  NT <- N * Tt
  scaling_matrix <- matrix(scaling, nrow = N, ncol = Tt)
  log_rate <- pmin(Mu_pm + Mu_pv / 2, 700)
  poisson <- sum(Y * Mu_pm - scaling_matrix * exp(log_rate))
  residual_second_moment <- sum((Mu_pm - B_pm)^2 + Mu_pv + B_pv)
  gaussian <- -(NT / 2) * log(2 * pi * sigma2) -
    residual_second_moment / (2 * sigma2)
  entropy_Mu <- 0.5 * sum(log(2 * pi * exp(1) * Mu_pv))
  poisson + gaussian + entropy_Mu
}


relative_rms_change <- function(new, old) {
  sqrt(mean((new - old)^2)) /
    max(sqrt(mean(old^2)), sqrt(.Machine$double.eps))
}


## Retained internal helper for callers that need a lead-SNP reconstruction.
build_fitted_leadSNP <- function(response, X, susiF.obj) {
  N <- nrow(X)
  Tt <- ncol(response)
  out <- matrix(colMeans(response), nrow = N, ncol = Tt, byrow = TRUE)
  if (length(active_pois_effects(susiF.obj)) == 0L) return(out)
  X_centered <- sweep(X, 2L, colMeans(X), "-")
  for (l in active_pois_effects(susiF.obj)) {
    lead <- which.max(susiF.obj$alpha[[l]])
    out <- out + tcrossprod(X_centered[, lead],
                            as.numeric(susiF.obj$fitted_func[[l]]))
  }
  out
}


pois_fsusie_diagnostic_plot <- function(Y, Mu_pm, susiF.obj,
                                         True_intensity, B_pm) {
  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)
  if (is.null(True_intensity)) {
    graphics::par(mfrow = c(2, 1))
    graphics::plot(log1p(Y), Mu_pm,
                   xlab = "log1p(Y)", ylab = "Mu_pm",
                   main = "Latent log-intensity recovery")
    graphics::abline(a = 0, b = 1)
    if (length(susiF.obj$fitted_func) >= 1L) {
      graphics::plot(susiF.obj$fitted_func[[1L]], type = "l",
                     main = "fSuSiE fitted_func[[1]]")
    }
  } else {
    graphics::par(mfrow = c(2, 2))
    graphics::plot(log1p(Y), Mu_pm, main = "Mu_pm vs log1p(Y)")
    graphics::abline(0, 1)
    mse <- mean((B_pm - True_intensity)^2)
    graphics::plot(True_intensity, B_pm,
                   main = sprintf("B_pm vs truth (MSE %.4f)", mse))
    graphics::abline(0, 1)
    if (length(susiF.obj$fitted_func) >= 1L) {
      graphics::plot(susiF.obj$fitted_func[[1L]], type = "l",
                     main = "fSuSiE fitted_func[[1]]")
    }
    if (length(susiF.obj$fitted_func) >= 2L) {
      graphics::plot(susiF.obj$fitted_func[[2L]], type = "l",
                     main = "fSuSiE fitted_func[[2]]")
    }
  }
}
