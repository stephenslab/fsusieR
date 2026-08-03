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
#' @param verbose Logical; report outer-iteration progress.
#' @param diagnostic_plot Logical; draw diagnostic plots each outer iteration.
#' @param True_intensity Optional N by T matrix of true latent log-intensities,
#'   used only by `diagnostic_plot`.
#' @param print Deprecated alias for `diagnostic_plot`.
#' @param Z Reserved for future dense-covariate support. Must be `NULL`.
#' @param update_Mu_each_iter Logical; update the latent Gaussian posterior at
#'   every outer iteration. Setting this to `FALSE` is intended for debugging.
#' @param s2 Positive initial value of the noise variance `sigma2`.
#' @param outer_tol Positive convergence tolerance for the maximum relative
#'   change in `Mu_pm`, `B_pm`, and `sigma2`.
#' @param stable_iterations Number of consecutive iterations below `outer_tol`
#'   required for convergence.
#' @param solver_failure Either `"error"` (the default) or `"log1p"`. The latter
#'   records a failed VGA row and uses a log1p fallback instead of stopping.
#'
#' @return A list containing posterior moments, fitted values, `sigma2`, an
#'   iteration-level `convergence_trace`, solver failures, and the inner
#'   `susiF.obj`. `partial_objective_trace` contains only the Poisson latent and
#'   Gaussian coupling terms. The legacy `elbo_trace` component is retained as
#'   an alias, but is not a complete joint ELBO.
#'
#' @details
#' The noise-variance update is
#' \deqn{\sigma^2 = \mathrm{mean}\{(\bar\mu-\bar B)^2 + V_\mu + V_B\}.}
#' The approximation to `V_B` uses centered `X` and includes both SNP-selection
#' uncertainty and the available post-processed effect-curve variance.
#'
#' This implementation is deliberately labelled a hybrid approximation. The
#' public [susiF()] update standardizes its response and re-estimates an inner
#' Gaussian residual variance, so its KL term is not on the scale of the outer
#' Poisson model. Consequently `partial_objective_trace` is a diagnostic and
#' convergence is based on parameter changes rather than a claimed joint ELBO.
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
                        outer_tol = 1e-3,
                        stable_iterations = 2L,
                        solver_failure = c("error", "log1p")) {

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
    thresh_lowcount = thresh_lowcount
  )

  est_effect_fm <- reconstruct_effect(susiF.obj, ncol(X), Tt)
  B_pm <- build_B_pm(Y_warm, X, est_effect_fm)
  B_pv <- compute_B_pv(susiF.obj, X, Tt)
  Mu_pm <- B_pm
  Mu_pv <- matrix(1 / Tt, nrow = N, ncol = Tt)
  sigma2 <- max(s2, 1e-8)

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
        invalid <- is.null(sol) || is.null(sol$m) || is.null(sol$v) ||
          !all(is.finite(sol$m)) || !all(is.finite(sol$v)) || any(sol$v <= 0)

        if (invalid) {
          if (is.null(failure_message)) failure_message <- "invalid VGA output"
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
      thresh_lowcount = thresh_lowcount
    )
    est_effect_fm <- reconstruct_effect(susiF.obj, ncol(X), Tt)
    B_pm <- build_B_pm(Mu_pm, X, est_effect_fm)
    B_pv <- compute_B_pv(susiF.obj, X, Tt)

    residual_component <- mean((Mu_pm - B_pm)^2)
    latent_variance_component <- mean(Mu_pv)
    structured_variance_component <- mean(B_pv)
    sigma2 <- residual_component + latent_variance_component +
      structured_variance_component
    if (!is.finite(sigma2) || sigma2 <= 0) {
      stop("The sigma2 update produced a non-positive or non-finite value.",
           call. = FALSE)
    }
    sigma2 <- max(sigma2, 1e-8)

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
    sigma2 = sigma2,
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
        "susiF re-estimates a residual variance on its internally standardized",
        "response instead of using the outer Poisson-FSuSiE sigma2."
      ),
      B_moments = paste(
        "Post-processed lead-SNP curves are alpha-distributed; B_pv includes",
        "selection and available curve-magnitude variance but remains approximate."
      ),
      objective = paste(
        "partial_objective_trace excludes the incompatible inner KL and is",
        "a diagnostic, not a complete joint ELBO."
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
                           thresh_lowcount) {
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
    cal_obj = FALSE,
    verbose = FALSE
  )
  if (is.null(fit)) {
    stop("susiF returned NULL (almost all wavelet coefficients were filtered).",
         call. = FALSE)
  }
  fit
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
