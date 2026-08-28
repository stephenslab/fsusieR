#' @title Binomial functional SuSiE
#'
#' @description
#' Fine-mapping of aggregated binary functional data with a binomial-logistic
#' likelihood, a Gaussian variational approximation for the latent log-odds,
#' and an inner [susiF()] fit for the structured genetic mean.
#'
#' @param Y N by T numeric matrix of integer success counts.
#' @param trials A positive integer scalar, a length-N vector, or an N by T
#'   matrix giving the corresponding trial totals. The success counts and
#'   trial totals, rather than their ratio, must be supplied.
#' @param X N by J numeric covariate matrix.
#' @param L Upper bound on the number of single effects.
#' @param reflect Logical; mirror observations before wavelet processing.
#'   Reflection is enabled automatically when `ncol(Y)` is not a power of two.
#' @param maxit_outer Maximum number of split-VA outer iterations.
#' @param maxit_inner Maximum number of inner [susiF()] iterations.
#' @param tol Positive convergence tolerance forwarded to [susiF()].
#' @param outer_tol Positive tolerance for the maximum relative change in the
#'   latent mean, structured mean, and splitting variance.
#' @param stable_iterations Number of consecutive stable outer iterations
#'   required for convergence.
#' @param sigma2_subcycles Non-negative number of conditional latent-moment and
#'   splitting-variance subcycles before each Gaussian fSuSiE update. The first
#'   outer iteration uses at most one subcycle.
#' @param s2 Positive initial splitting variance on the latent log-odds scale.
#' @param estimate_sigma2 Logical; estimate the splitting variance when `TRUE`.
#'   Set this to `FALSE` to keep it fixed at `s2`, particularly when trial
#'   totals are mostly one and the variance is weakly identified.
#' @param warm_start_sigma2 Logical; initialize the splitting variance from an
#'   empirical-logit fSuSiE fit when possible. Ignored if
#'   `estimate_sigma2 = FALSE`.
#' @param quadrature_order Positive integer number of Gauss--Hermite nodes used
#'   for logistic-normal expectations.
#' @param solver_failure Either `"error"` or `"empirical_logit"`. The latter
#'   records a failed local update and substitutes a finite empirical logit.
#' @param control_mixsqp,nullweight,cov_lev,min_purity,cor_small,post_processing,thresh_lowcount
#'   Forwarded to [susiF()]. `post_processing = "none"` is not supported.
#' @param filter.number,family Wavelet filter passed to [susiF()].
#' @param verbose Logical; report outer-iteration progress.
#' @param diagnostic_plot Logical; draw diagnostic plots each outer iteration.
#' @param True_probability Optional N by T matrix of true probabilities, used
#'   only for diagnostic plots.
#' @param Z Reserved for future dense-covariate support. Must be `NULL`.
#' @param update_Mu_each_iter Logical; update latent Gaussian moments in every
#'   outer iteration. Setting this to `FALSE` is intended for debugging.
#'
#' @return A list with latent and structured posterior moments, fitted
#'   probabilities and success counts, the splitting variance, convergence
#'   diagnostics, and the inner `susiF.obj`.
#'
#' @details
#' The fitted model is
#' \deqn{Y_{it} \mid \mu_{it} \sim
#'       \mathrm{Binomial}(n_{it},\mathrm{logit}^{-1}(\mu_{it})),}
#' \deqn{\mu_{it} \mid B_{it} \sim N(B_{it},\sigma^2),}
#' where `B` has the sum-of-single-functions prior implemented by [susiF()].
#' The local variational distribution is Gaussian. Logistic-normal
#' expectations and their derivatives are evaluated with Gauss--Hermite
#' quadrature. The Gaussian block centers but does not scale the latent
#' response and holds its residual variance fixed at the current `sigma2`.
#'
#' Estimating `sigma2` from Bernoulli data can be weakly identified. At a zero
#' latent prior mean, the marginal success probability is one half for every
#' symmetric logistic-normal variance. Replicated trials add the
#' extra-binomial component needed to learn latent heterogeneity. The function
#' warns when all trial totals equal one.
#'
#' @seealso [susiF()], [Pois_fSuSiE()], [vga_binomial_solver()]
#' @export
Binomial_fSuSiE <- function(Y,
                            trials,
                            X,
                            L = 3,
                            reflect = FALSE,
                            maxit_outer = 20,
                            maxit_inner = 10,
                            tol = 1e-3,
                            outer_tol = 1e-3,
                            stable_iterations = 2L,
                            sigma2_subcycles = 5L,
                            s2 = 0.5,
                            estimate_sigma2 = TRUE,
                            warm_start_sigma2 = TRUE,
                            quadrature_order = 20L,
                            solver_failure = c("error", "empirical_logit"),
                            control_mixsqp = list(verbose = FALSE,
                                                  eps = 1e-6,
                                                  numiter.em = 4),
                            nullweight = 0.10,
                            cov_lev = 0.95,
                            min_purity = 0.5,
                            cor_small = FALSE,
                            post_processing = "smash",
                            thresh_lowcount = 0,
                            filter.number = 10,
                            family = "DaubLeAsymm",
                            verbose = TRUE,
                            diagnostic_plot = FALSE,
                            True_probability = NULL,
                            Z = NULL,
                            update_Mu_each_iter = TRUE) {
  if (missing(X) || is.null(X)) stop("Please provide X matrix", call. = FALSE)
  if (missing(Y) || is.null(Y)) stop("Please provide Y matrix", call. = FALSE)
  if (missing(trials) || is.null(trials)) {
    stop("Please provide binomial trial totals in `trials`.", call. = FALSE)
  }
  if (!is.null(Z)) {
    stop("Z covariate handling is not implemented; please pass Z = NULL.",
         call. = FALSE)
  }
  post_processing <- match.arg(post_processing, c("smash", "TI", "HMM"))
  solver_failure <- match.arg(solver_failure)

  checked <- validate_binomial_fsusie_inputs(Y, trials, X, L)
  Y <- checked$Y
  trials <- checked$trials
  X <- checked$X
  L <- checked$L

  positive_scalar <- function(value) {
    length(value) == 1L && is.finite(value) && value > 0
  }
  positive_integer <- function(value) {
    positive_scalar(value) && value == as.integer(value)
  }
  nonnegative_integer <- function(value) {
    length(value) == 1L && is.finite(value) && value >= 0 &&
      value == as.integer(value)
  }
  scalar_logical <- function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
  }
  if (!positive_integer(maxit_outer) || !positive_integer(maxit_inner) ||
      !positive_integer(stable_iterations) ||
      !positive_integer(quadrature_order) || quadrature_order < 3L ||
      !nonnegative_integer(sigma2_subcycles) || !positive_scalar(tol) ||
      !positive_scalar(outer_tol) || !positive_scalar(s2)) {
    stop(paste(
      "Iteration counts, tolerances, and `s2` must be strictly positive;",
      "`quadrature_order` must be at least three and `sigma2_subcycles`",
      "must be a non-negative integer."
    ), call. = FALSE)
  }
  if (!scalar_logical(estimate_sigma2) ||
      !scalar_logical(warm_start_sigma2) ||
      !scalar_logical(update_Mu_each_iter)) {
    stop(paste(
      "`estimate_sigma2`, `warm_start_sigma2`, and",
      "`update_Mu_each_iter` must be TRUE or FALSE."
    ), call. = FALSE)
  }
  maxit_outer <- as.integer(maxit_outer)
  maxit_inner <- as.integer(maxit_inner)
  stable_iterations <- as.integer(stable_iterations)
  quadrature_order <- as.integer(quadrature_order)
  sigma2_subcycles <- as.integer(sigma2_subcycles)

  if (isTRUE(estimate_sigma2) && all(trials == 1)) {
    warning(paste(
      "All trial totals are one. The binomial splitting variance can be",
      "weakly identified; consider `estimate_sigma2 = FALSE` or supply",
      "aggregated success and trial counts."
    ), call. = FALSE)
  }

  original_grid_size <- ncol(Y)
  true_probability_internal <- NULL
  if (!is.null(True_probability)) {
    true_probability_internal <- as.matrix(True_probability)
    if (!is.numeric(true_probability_internal) ||
        !identical(dim(true_probability_internal), dim(Y)) ||
        any(!is.finite(true_probability_internal)) ||
        any(true_probability_internal < 0 | true_probability_internal > 1)) {
      stop("`True_probability` must be a finite matrix in [0, 1] with dim(Y).",
           call. = FALSE)
    }
  }

  must_reflect <- (log2(original_grid_size) %% 1) != 0
  if (isTRUE(reflect) || must_reflect) {
    reflected_Y <- lapply(seq_len(nrow(Y)), function(i) reflect_vec(Y[i, ]))
    reflected_trials <- lapply(
      seq_len(nrow(trials)),
      function(i) reflect_vec(trials[i, ])
    )
    Y <- do.call(rbind, lapply(reflected_Y, `[[`, "x"))
    trials <- do.call(rbind, lapply(reflected_trials, `[[`, "x"))
    output_index <- reflected_Y[[1L]]$idx
    if (!is.null(true_probability_internal)) {
      true_probability_internal <- do.call(
        rbind,
        lapply(seq_len(nrow(true_probability_internal)), function(i) {
          reflect_vec(true_probability_internal[i, ])$x
        })
      )
    }
  } else {
    output_index <- seq_len(original_grid_size)
  }

  n_sample <- nrow(Y)
  n_grid <- ncol(Y)
  inverse_dwt <- pois_inverse_dwt_matrix(
    n_grid = n_grid,
    filter.number = filter.number,
    family = family
  )
  empirical_probability <- (Y + 0.5) / (trials + 1)
  empirical_logit <- stats::qlogis(empirical_probability)

  susiF.obj <- fit_binomial_susif(
    response = empirical_logit,
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
  postprocessed_effect_fm <- reconstruct_effect(
    susiF.obj, ncol(X), n_grid
  )
  est_effect_fm <- postprocessed_effect_fm
  B_pm <- build_B_pm(empirical_logit, X, postprocessed_effect_fm)
  B_pv <- compute_B_pv(susiF.obj, X, n_grid)
  Mu_pm <- B_pm
  Mu_pv <- matrix(1 / n_grid, nrow = n_sample, ncol = n_grid)

  warm_start <- estimate_binomial_sigma2_warm_start(
    response_logit = empirical_logit,
    fitted_logit = B_pm,
    fitted_variance = B_pv,
    trials = trials,
    fitted_probability = stats::plogis(B_pm),
    fallback_sigma2 = s2,
    enabled = isTRUE(estimate_sigma2) && isTRUE(warm_start_sigma2)
  )
  sigma2 <- if (isTRUE(estimate_sigma2)) warm_start$sigma2 else s2
  if (!isTRUE(estimate_sigma2)) {
    warm_start$sigma2 <- s2
    warm_start$method <- "fixed user-supplied s2"
    warm_start$used <- FALSE
  }
  if (verbose) {
    message(sprintf("Initial sigma2 = %.5g (%s)", sigma2,
                    warm_start$method))
  }

  partial_objective_trace <- numeric(0)
  convergence_trace <- vector("list", maxit_outer)
  sigma2_subcycle_trace <- vector("list", maxit_outer)
  solver_failures <- data.frame(
    iteration = integer(0),
    subcycle = integer(0),
    row = integer(0),
    message = character(0)
  )
  stable_count <- 0L
  converged <- FALSE

  for (iter in seq_len(maxit_outer)) {
    if (verbose) message("=== Binomial_fSuSiE outer iter ", iter, " ===")
    old_Mu_pm <- Mu_pm
    old_B_pm <- B_pm
    old_sigma2 <- sigma2

    requested_subcycles <- if (!isTRUE(estimate_sigma2) ||
                               sigma2_subcycles == 0L) {
      0L
    } else if (iter == 1L) {
      1L
    } else {
      sigma2_subcycles
    }
    should_update_mu <- isTRUE(update_Mu_each_iter) || iter == 1L
    n_mu_updates <- if (should_update_mu) max(1L, requested_subcycles) else 0L
    subcycle_sigma2 <- numeric(0)
    subcycle_objective <- numeric(0)

    if (n_mu_updates > 0L) {
      for (subcycle_index in seq_len(n_mu_updates)) {
        subcycle_label <- if (requested_subcycles > 0L) {
          subcycle_index
        } else {
          0L
        }
        for (i in seq_len(n_sample)) {
          failure_message <- NULL
          sol <- tryCatch(
            vga_binomial_solver(
              init_val = Mu_pm[i, ],
              x = Y[i, ],
              trials = trials[i, ],
              beta = B_pm[i, ],
              sigma2 = sigma2,
              quadrature_order = quadrature_order
            ),
            error = function(error) {
              failure_message <<- conditionMessage(error)
              NULL
            }
          )
          invalid <- !vga_binomial_solution_is_valid(
            sol,
            x = Y[i, ],
            trials = trials[i, ],
            beta = B_pm[i, ],
            sigma2 = sigma2,
            quadrature_order = quadrature_order,
            tol = 1e-4
          )
          if (invalid) {
            if (is.null(failure_message)) {
              failure_message <- paste(
                "invalid VGA output: moments do not satisfy the binomial",
                "quadrature score equations"
              )
            }
            solver_failures <- rbind(
              solver_failures,
              data.frame(
                iteration = iter,
                subcycle = subcycle_label,
                row = i,
                message = failure_message
              )
            )
            if (solver_failure == "error") {
              stop(sprintf(
                "Binomial VGA failed at iteration %d, subcycle %d, row %d: %s",
                iter, subcycle_label, i, failure_message
              ), call. = FALSE)
            }
            Mu_pm[i, ] <- empirical_logit[i, ]
            Mu_pv[i, ] <- sigma2
          } else {
            Mu_pm[i, ] <- sol$m
            Mu_pv[i, ] <- sol$v
          }
        }

        if (requested_subcycles > 0L) {
          latent_sigma2_update <- pois_sigma2_update(
            Mu_pm, Mu_pv, B_pm, B_pv
          )
          sigma2 <- latent_sigma2_update$sigma2
          subcycle_sigma2[subcycle_index] <- sigma2
          subcycle_objective[subcycle_index] <-
            binomial_fsusie_partial_objective(
              Y, trials, Mu_pm, Mu_pv, B_pm, B_pv, sigma2,
              quadrature_order
            )
        }
      }
    }

    sigma2_subcycle_trace[[iter]] <- data.frame(
      iteration = rep.int(iter, length(subcycle_sigma2)),
      subcycle = seq_along(subcycle_sigma2),
      sigma2 = subcycle_sigma2,
      partial_objective = subcycle_objective
    )

    sigma2_for_inner <- sigma2
    susiF.obj <- fit_binomial_susif(
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
    postprocessed_effect_fm <- reconstruct_effect(
      susiF.obj, ncol(X), n_grid
    )
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
    if (isTRUE(estimate_sigma2)) sigma2 <- sigma2_update$sigma2
    residual_component <- sigma2_update$residual_component
    latent_variance_component <- sigma2_update$latent_variance_component
    structured_variance_component <-
      sigma2_update$structured_variance_component

    partial_objective <- binomial_fsusie_partial_objective(
      Y, trials, Mu_pm, Mu_pv, B_pm, B_pv, sigma2, quadrature_order
    )
    partial_objective_trace <- c(
      partial_objective_trace,
      partial_objective
    )

    delta_Mu <- relative_rms_change(Mu_pm, old_Mu_pm)
    delta_B <- relative_rms_change(B_pm, old_B_pm)
    delta_sigma2 <- if (isTRUE(estimate_sigma2)) {
      abs(sigma2 - old_sigma2) / max(old_sigma2, 1e-8)
    } else {
      0
    }
    maximum_change <- max(delta_Mu, delta_B, delta_sigma2)
    convergence_trace[[iter]] <- data.frame(
      iteration = iter,
      sigma2 = sigma2,
      sigma2_start = old_sigma2,
      sigma2_subcycles = length(subcycle_sigma2),
      inner_sigma2 = sigma2_for_inner,
      sigma2_estimated = isTRUE(estimate_sigma2),
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
        paste0("  sigma2=%.5g  max relative change=%.3g  ",
               "subcycles=%d  partial objective=%.5g"),
        sigma2, maximum_change, length(subcycle_sigma2), partial_objective
      ))
    }
    if (diagnostic_plot) {
      binomial_fsusie_diagnostic_plot(
        Y, trials, Mu_pm, Mu_pv, susiF.obj, true_probability_internal, B_pm,
        quadrature_order
      )
    }
    if (stable_count >= stable_iterations) {
      converged <- TRUE
      break
    }
  }

  n_iter <- length(partial_objective_trace)
  convergence_trace <- do.call(rbind, convergence_trace[seq_len(n_iter)])
  sigma2_subcycle_trace <- do.call(
    rbind,
    sigma2_subcycle_trace[seq_len(n_iter)]
  )
  if (!converged && verbose) {
    warning("Binomial_fSuSiE did not converge in ", maxit_outer,
            " iterations", call. = FALSE)
  }
  susiF.obj <- update_cal_pip(susiF.obj)

  posterior_terms <- binomial_logistic_normal_terms(
    as.numeric(Mu_pm), as.numeric(Mu_pv), quadrature_order
  )
  posterior_probability_mean <- matrix(
    posterior_terms$probability,
    nrow = n_sample,
    ncol = n_grid
  )
  crop <- function(value) value[, output_index, drop = FALSE]
  fitted_latent <- crop(B_pm)
  fitted_probability <- stats::plogis(fitted_latent)
  output_trials <- crop(trials)

  list(
    family = "binomial",
    successes = crop(Y),
    trials = output_trials,
    Mu_pm = crop(Mu_pm),
    Mu_pv = crop(Mu_pv),
    B_pm = crop(B_pm),
    B_pv = crop(B_pv),
    est_effect_fm = est_effect_fm[, output_index, drop = FALSE],
    postprocessed_effect_fm =
      postprocessed_effect_fm[, output_index, drop = FALSE],
    sigma2 = sigma2,
    estimate_sigma2 = isTRUE(estimate_sigma2),
    sigma2_used_for_final_susif = susiF.obj$sigma2,
    sigma2_initial = warm_start$sigma2,
    sigma2_warm_start = warm_start,
    fitted_latent = fitted_latent,
    fitted = fitted_probability,
    fitted_probability = fitted_probability,
    posterior_probability_mean = crop(posterior_probability_mean),
    posterior_success_mean =
      output_trials * crop(posterior_probability_mean),
    quadrature_order = quadrature_order,
    partial_objective_trace = partial_objective_trace,
    elbo_trace = partial_objective_trace,
    elbo_trace_is_partial = TRUE,
    convergence_trace = convergence_trace,
    sigma2_subcycle_trace = sigma2_subcycle_trace,
    converged = converged,
    n_iter = n_iter,
    solver_failures = solver_failures,
    susiF.obj = susiF.obj,
    approximation = list(
      exact_joint_cavi = FALSE,
      likelihood = paste(
        "Logistic-normal expectations use",
        quadrature_order,
        "node Gauss-Hermite quadrature."
      ),
      B_moments = paste(
        "B_pm and B_pv use raw alpha, conditional wavelet means, and",
        "conditional wavelet variances; post-processing is output-only."
      ),
      objective = paste(
        "partial_objective_trace excludes the fSuSiE prior/KL contribution",
        "and is a diagnostic, not the complete joint ELBO."
      )
    ),
    reflected_internal_grid = n_grid != original_grid_size,
    output_grid_index = output_index
  )
}


#' @title Solve the binomial logistic-normal VGA problem
#'
#' @description
#' Optimizes Gaussian variational posterior means and variances for independent
#' binomial observations with Gaussian prior means and variances. Stable
#' Gauss--Hermite quadrature is used for logistic-normal expectations.
#'
#' @param init_val Initial posterior-mean vector.
#' @param x Integer success-count vector.
#' @param trials Integer trial-total vector.
#' @param beta Gaussian prior mean, scalar or vector.
#' @param sigma2 Positive Gaussian prior variance, scalar or vector.
#' @param quadrature_order Number of Gauss--Hermite nodes.
#' @param maxiter Maximum L-BFGS-B iterations.
#' @param tol Projected-gradient tolerance.
#'
#' @return A list containing posterior mean `m`, posterior variance `v`, the
#'   optimizer convergence code, and the optimized local variational objective.
#'
#' @examples
#' successes <- c(0, 3, 8, 10)
#' totals <- rep(10, length(successes))
#' vga_binomial_solver(
#'   init_val = stats::qlogis((successes + 0.5) / (totals + 1)),
#'   x = successes,
#'   trials = totals,
#'   beta = 0,
#'   sigma2 = 0.5
#' )
#' @export
vga_binomial_solver <- function(init_val,
                                x,
                                trials,
                                beta,
                                sigma2,
                                quadrature_order = 20L,
                                maxiter = 300L,
                                tol = 1e-7) {
  n_value <- length(x)
  recycle_to_n <- function(value, name) {
    if (length(value) == 1L) value <- rep(value, n_value)
    if (length(value) != n_value) {
      stop(sprintf("`%s` must have length one or length(x).", name),
           call. = FALSE)
    }
    as.numeric(value)
  }
  init_val <- recycle_to_n(init_val, "init_val")
  x <- recycle_to_n(x, "x")
  trials <- recycle_to_n(trials, "trials")
  beta <- recycle_to_n(beta, "beta")
  sigma2 <- recycle_to_n(sigma2, "sigma2")
  if (n_value < 1L || any(!is.finite(init_val)) || any(!is.finite(x)) ||
      any(!is.finite(trials)) || any(!is.finite(beta)) ||
      any(!is.finite(sigma2)) || any(x < 0) || any(trials < 0) ||
      any(x > trials) || any(abs(x - round(x)) > 1e-8) ||
      any(abs(trials - round(trials)) > 1e-8) || any(sigma2 <= 0) ||
      length(quadrature_order) != 1L || !is.finite(quadrature_order) ||
      quadrature_order < 3L || quadrature_order != as.integer(quadrature_order) ||
      length(maxiter) != 1L || !is.finite(maxiter) || maxiter < 1L ||
      maxiter != as.integer(maxiter) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("Invalid binomial VGA inputs.", call. = FALSE)
  }
  quadrature_order <- as.integer(quadrature_order)

  m <- beta
  v <- sigma2
  active <- trials > 0
  if (!any(active)) {
    return(list(m = m, v = v, convergence = 0L, objective = 0))
  }
  xa <- x[active]
  na <- trials[active]
  ba <- beta[active]
  s2a <- sigma2[active]
  lower_m <- ba + s2a * (xa - na)
  upper_m <- ba + s2a * xa
  initial_m <- pmin(pmax(init_val[active], lower_m), upper_m)
  initial_v <- 1 / (1 / s2a + na / 4)
  minimum_v <- pmax(s2a * 1e-12, .Machine$double.xmin * 1e4)
  initial_log_v <- log(pmin(pmax(initial_v, minimum_v), s2a))
  n_active <- length(xa)

  evaluate <- function(theta, gradient = FALSE) {
    current_m <- theta[seq_len(n_active)]
    current_log_v <- theta[n_active + seq_len(n_active)]
    current_v <- exp(current_log_v)
    terms <- binomial_logistic_normal_terms(
      current_m, current_v, quadrature_order
    )
    value <- sum(
      xa * current_m - na * terms$softplus -
        ((current_m - ba)^2 + current_v) / (2 * s2a) +
        current_log_v / 2
    )
    if (!gradient) return(-value)
    gradient_m <- xa - na * terms$probability -
      (current_m - ba) / s2a
    gradient_log_v <- 0.5 - current_v / (2 * s2a) -
      na * terms$softplus_log_variance_derivative
    -c(gradient_m, gradient_log_v)
  }

  fit <- stats::optim(
    par = c(initial_m, initial_log_v),
    fn = evaluate,
    gr = function(theta) evaluate(theta, gradient = TRUE),
    method = "L-BFGS-B",
    lower = c(lower_m, log(minimum_v)),
    upper = c(upper_m, log(s2a)),
    control = list(maxit = as.integer(maxiter), pgtol = tol, factr = 1e3)
  )
  candidate_m <- fit$par[seq_len(n_active)]
  candidate_v <- exp(fit$par[n_active + seq_len(n_active)])
  m[active] <- candidate_m
  v[active] <- candidate_v
  result <- list(
    m = m,
    v = v,
    convergence = fit$convergence,
    message = fit$message,
    objective = -fit$value
  )
  if (!vga_binomial_solution_is_valid(
    result, x, trials, beta, sigma2, quadrature_order,
    tol = max(20 * tol, 2e-5)
  )) {
    ## A single block optimization can occasionally stop because the sum of
    ## many already-converged objectives changes too little while one
    ## coordinate still has a noticeable score. Retry scalar problems so that
    ## each coordinate has its own convergence criterion.
    scalar_objective <- numeric(n_active)
    scalar_convergence <- integer(n_active)
    for (j in seq_len(n_active)) {
      evaluate_scalar <- function(theta, gradient = FALSE) {
        current_m <- theta[1L]
        current_log_v <- theta[2L]
        current_v <- exp(current_log_v)
        terms <- binomial_logistic_normal_terms(
          current_m, current_v, quadrature_order
        )
        value <- xa[j] * current_m - na[j] * terms$softplus -
          ((current_m - ba[j])^2 + current_v) / (2 * s2a[j]) +
          current_log_v / 2
        if (!gradient) return(-value)
        gradient_m <- xa[j] - na[j] * terms$probability -
          (current_m - ba[j]) / s2a[j]
        gradient_log_v <- 0.5 - current_v / (2 * s2a[j]) -
          na[j] * terms$softplus_log_variance_derivative
        -c(gradient_m, gradient_log_v)
      }
      scalar_fit <- stats::optim(
        par = c(candidate_m[j], log(candidate_v[j])),
        fn = evaluate_scalar,
        gr = function(theta) evaluate_scalar(theta, gradient = TRUE),
        method = "L-BFGS-B",
        lower = c(lower_m[j], log(minimum_v[j])),
        upper = c(upper_m[j], log(s2a[j])),
        control = list(
          maxit = max(500L, 2L * as.integer(maxiter)),
          pgtol = tol,
          factr = 10
        )
      )
      candidate_m[j] <- scalar_fit$par[1L]
      candidate_v[j] <- exp(scalar_fit$par[2L])
      scalar_objective[j] <- -scalar_fit$value
      scalar_convergence[j] <- scalar_fit$convergence
    }
    m[active] <- candidate_m
    v[active] <- candidate_v
    result <- list(
      m = m,
      v = v,
      convergence = max(scalar_convergence),
      message = "scalar L-BFGS-B fallback",
      objective = sum(scalar_objective)
    )
    if (!vga_binomial_solution_is_valid(
      result, x, trials, beta, sigma2, quadrature_order,
      tol = max(20 * tol, 2e-5)
    )) {
      stop(paste(
        "Block and scalar L-BFGS-B failed to solve the binomial VGA",
        "score equations."
      ), call. = FALSE)
    }
  }
  result
}


binomial_gh_cache <- new.env(parent = emptyenv())


binomial_gauss_hermite_rule <- function(order) {
  if (length(order) != 1L || !is.finite(order) || order < 3L ||
      order != as.integer(order)) {
    stop("`order` must be an integer of at least three.", call. = FALSE)
  }
  order <- as.integer(order)
  key <- as.character(order)
  if (exists(key, envir = binomial_gh_cache, inherits = FALSE)) {
    return(get(key, envir = binomial_gh_cache, inherits = FALSE))
  }
  jacobi <- matrix(0, nrow = order, ncol = order)
  off_diagonal <- sqrt(seq_len(order - 1L) / 2)
  jacobi[cbind(seq_len(order - 1L), 2:order)] <- off_diagonal
  jacobi[cbind(2:order, seq_len(order - 1L))] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  index <- base::order(decomposition$values)
  rule <- list(
    nodes = decomposition$values[index],
    weights = decomposition$vectors[1L, index]^2
  )
  assign(key, rule, envir = binomial_gh_cache)
  rule
}


binomial_softplus <- function(value) {
  result <- numeric(length(value))
  positive <- value > 0
  result[positive] <- value[positive] + log1p(exp(-value[positive]))
  result[!positive] <- log1p(exp(value[!positive]))
  result
}


binomial_logistic_normal_terms <- function(m, v, quadrature_order = 20L) {
  m <- as.numeric(m)
  v <- as.numeric(v)
  if (length(m) != length(v) || any(!is.finite(m)) ||
      any(!is.finite(v)) || any(v <= 0)) {
    stop("Invalid logistic-normal moments.", call. = FALSE)
  }
  rule <- binomial_gauss_hermite_rule(quadrature_order)
  root_variance <- sqrt(2) * sqrt(v)
  z <- matrix(m, nrow = length(m), ncol = length(rule$nodes)) +
    tcrossprod(root_variance, rule$nodes)
  probability <- stats::plogis(z)
  weights <- matrix(
    rule$weights,
    nrow = length(m),
    ncol = length(rule$weights),
    byrow = TRUE
  )
  node_matrix <- matrix(
    rule$nodes,
    nrow = length(m),
    ncol = length(rule$nodes),
    byrow = TRUE
  )
  list(
    softplus = rowSums(weights * matrix(
      binomial_softplus(as.numeric(z)),
      nrow = nrow(z),
      ncol = ncol(z)
    )),
    probability = rowSums(weights * probability),
    probability_variance = rowSums(weights * probability * (1 - probability)),
    softplus_log_variance_derivative =
      rowSums(weights * probability * node_matrix) * sqrt(v / 2)
  )
}


vga_binomial_solution_is_valid <- function(result,
                                            x,
                                            trials,
                                            beta,
                                            sigma2,
                                            quadrature_order = 20L,
                                            tol = 1e-4) {
  n_value <- length(x)
  recycle_to_n <- function(value) {
    if (length(value) == 1L) rep(value, n_value) else value
  }
  trials <- recycle_to_n(trials)
  beta <- recycle_to_n(beta)
  sigma2 <- recycle_to_n(sigma2)
  if (is.null(result) || is.null(result$m) || is.null(result$v) ||
      length(result$m) != n_value || length(result$v) != n_value ||
      length(trials) != n_value || length(beta) != n_value ||
      length(sigma2) != n_value || any(!is.finite(result$m)) ||
      any(!is.finite(result$v)) || any(!is.finite(x)) ||
      any(!is.finite(trials)) || any(!is.finite(beta)) ||
      any(!is.finite(sigma2)) || any(result$v <= 0) || any(sigma2 <= 0) ||
      any(result$v > sigma2 * (1 + tol))) {
    return(FALSE)
  }
  inactive <- trials == 0
  if (any(inactive)) {
    inactive_valid <- abs(result$m[inactive] - beta[inactive]) <=
      tol * (1 + abs(beta[inactive])) &
      abs(result$v[inactive] - sigma2[inactive]) <=
      tol * (1 + sigma2[inactive])
    if (!all(inactive_valid)) return(FALSE)
  }
  active <- !inactive
  if (!any(active)) return(TRUE)
  terms <- binomial_logistic_normal_terms(
    result$m[active], result$v[active], quadrature_order
  )
  mean_score <- x[active] - trials[active] * terms$probability -
    (result$m[active] - beta[active]) / sigma2[active]
  variance_score <- 0.5 - result$v[active] / (2 * sigma2[active]) -
    trials[active] * terms$softplus_log_variance_derivative
  mean_scale <- 1 + x[active] + trials[active] * terms$probability +
    abs((result$m[active] - beta[active]) / sigma2[active])
  variance_scale <- 1 + result$v[active] / (2 * sigma2[active]) +
    abs(trials[active] * terms$softplus_log_variance_derivative)
  all(is.finite(mean_score)) && all(is.finite(variance_score)) &&
    max(abs(mean_score) / mean_scale) <= tol &&
    max(abs(variance_score) / variance_scale) <= tol
}


validate_binomial_fsusie_inputs <- function(Y, trials, X, L) {
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
  if (length(trials) == 1L) {
    trials <- matrix(trials, nrow = nrow(Y), ncol = ncol(Y))
  } else if (is.null(dim(trials)) && length(trials) == nrow(Y)) {
    trials <- matrix(trials, nrow = nrow(Y), ncol = ncol(Y))
  } else {
    trials <- as.matrix(trials)
  }
  if (!is.numeric(trials) || !identical(dim(trials), dim(Y))) {
    stop(paste(
      "`trials` must be a scalar, a length-nrow(Y) vector, or a numeric",
      "matrix with dim(Y)."
    ), call. = FALSE)
  }
  if (any(!is.finite(Y)) || any(!is.finite(trials)) ||
      any(!is.finite(X))) {
    stop("`Y`, `trials`, and `X` must contain only finite values.",
         call. = FALSE)
  }
  if (any(Y < 0) || any(abs(Y - round(Y)) > 1e-8)) {
    stop("`Y` must contain non-negative integer success counts.",
         call. = FALSE)
  }
  if (any(trials <= 0) || any(abs(trials - round(trials)) > 1e-8)) {
    stop("`trials` must contain strictly positive integer trial totals.",
         call. = FALSE)
  }
  if (any(Y > trials)) {
    stop("Every success count in `Y` must be no larger than `trials`.",
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
  storage.mode(Y) <- "double"
  storage.mode(trials) <- "double"
  storage.mode(X) <- "double"
  list(Y = Y, trials = trials, X = X, L = as.integer(L))
}


estimate_binomial_sigma2_warm_start <- function(response_logit,
                                                fitted_logit,
                                                fitted_variance,
                                                trials,
                                                fitted_probability,
                                                fallback_sigma2,
                                                enabled = TRUE) {
  fallback_sigma2 <- max(as.numeric(fallback_sigma2), 1e-8)
  fallback <- list(
    sigma2 = fallback_sigma2,
    method = "user-supplied s2 fallback",
    used = FALSE,
    transformed_residual_mse = NA_real_,
    structured_variance = NA_real_,
    binomial_sampling_variance = NA_real_,
    untruncated_sigma2 = NA_real_
  )
  if (!isTRUE(enabled)) {
    fallback$method <- "user-supplied s2"
    return(fallback)
  }
  valid <- is.matrix(response_logit) && is.matrix(fitted_logit) &&
    is.matrix(fitted_variance) && is.matrix(trials) &&
    is.matrix(fitted_probability) &&
    identical(dim(response_logit), dim(fitted_logit)) &&
    identical(dim(response_logit), dim(fitted_variance)) &&
    identical(dim(response_logit), dim(trials)) &&
    identical(dim(response_logit), dim(fitted_probability)) &&
    all(is.finite(response_logit)) && all(is.finite(fitted_logit)) &&
    all(is.finite(fitted_variance)) && all(fitted_variance >= 0) &&
    all(is.finite(trials)) && all(trials > 0) &&
    all(is.finite(fitted_probability)) &&
    all(fitted_probability > 0 & fitted_probability < 1)
  if (!valid) return(fallback)

  probability_floor <- 1 / (2 * (trials + 1))
  stable_probability <- pmin(
    pmax(fitted_probability, probability_floor),
    1 - probability_floor
  )
  transformed_residual_mse <- mean((response_logit - fitted_logit)^2)
  structured_variance <- mean(fitted_variance)
  binomial_sampling_variance <- mean(
    1 / (trials * stable_probability * (1 - stable_probability))
  )
  untruncated_sigma2 <- transformed_residual_mse + structured_variance -
    binomial_sampling_variance
  if (!is.finite(untruncated_sigma2) || untruncated_sigma2 <= 1e-8) {
    return(fallback)
  }
  list(
    sigma2 = untruncated_sigma2,
    method = "empirical-logit fSuSiE delta correction",
    used = TRUE,
    transformed_residual_mse = transformed_residual_mse,
    structured_variance = structured_variance,
    binomial_sampling_variance = binomial_sampling_variance,
    untruncated_sigma2 = untruncated_sigma2
  )
}


fit_binomial_susif <- function(response,
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
  fit_pois_susif(
    response = response,
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
    standardize_Y = standardize_Y,
    residual_variance = residual_variance
  )
}


binomial_fsusie_partial_objective <- function(Y,
                                               trials,
                                               Mu_pm,
                                               Mu_pv,
                                               B_pm,
                                               B_pv,
                                               sigma2,
                                               quadrature_order = 20L) {
  if (!is.finite(sigma2) || sigma2 <= 0 || any(Mu_pv <= 0)) return(-Inf)
  terms <- binomial_logistic_normal_terms(
    as.numeric(Mu_pm), as.numeric(Mu_pv), quadrature_order
  )
  binomial <- sum(
    lchoose(as.numeric(trials), as.numeric(Y)) +
      as.numeric(Y) * as.numeric(Mu_pm) -
      as.numeric(trials) * terms$softplus
  )
  n_value <- length(Y)
  residual_second_moment <- sum((Mu_pm - B_pm)^2 + Mu_pv + B_pv)
  gaussian <- -(n_value / 2) * log(2 * pi * sigma2) -
    residual_second_moment / (2 * sigma2)
  entropy_Mu <- 0.5 * sum(log(2 * pi * exp(1) * Mu_pv))
  binomial + gaussian + entropy_Mu
}


binomial_fsusie_diagnostic_plot <- function(Y,
                                             trials,
                                             Mu_pm,
                                             Mu_pv,
                                             susiF.obj,
                                             True_probability,
                                             B_pm,
                                             quadrature_order) {
  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)
  posterior_terms <- binomial_logistic_normal_terms(
    as.numeric(Mu_pm),
    as.numeric(Mu_pv),
    quadrature_order
  )
  posterior_probability <- matrix(
    posterior_terms$probability,
    nrow = nrow(Y),
    ncol = ncol(Y)
  )
  if (is.null(True_probability)) {
    graphics::par(mfrow = c(2, 1))
    graphics::plot(Y / trials, posterior_probability,
                   xlab = "Observed success proportion",
                   ylab = "Posterior probability",
                   main = "Binomial latent-probability fit")
    graphics::abline(a = 0, b = 1)
    if (length(susiF.obj$fitted_func) >= 1L) {
      graphics::plot(susiF.obj$fitted_func[[1L]], type = "l",
                     main = "fSuSiE fitted_func[[1]]")
    }
  } else {
    graphics::par(mfrow = c(2, 2))
    graphics::plot(Y / trials, posterior_probability,
                   main = "Posterior vs observed probability")
    graphics::abline(0, 1)
    mse <- mean((stats::plogis(B_pm) - True_probability)^2)
    graphics::plot(True_probability, stats::plogis(B_pm),
                   main = sprintf("Fitted vs truth (MSE %.4f)", mse))
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
