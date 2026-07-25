# ash_hmm.R
#
# A dependency-free implementation of an adaptive-shrinkage hidden Markov
# model for heteroskedastic normal observations.
#
# Model:
#   y_t | beta_t ~ N(beta_t, se_t^2)
#   Q_t is a finite-state Markov chain with transition matrix A
#   beta_t | Q_t=m ~ sum_l rho[m,l] G_ml
#
# In signed mode, G_ml is N(mu[m], prior_sd[l]^2). In the default one-sided
# mode, continuous G_ml are those normals truncated to [0, Inf); the zero-scale
# component remains a point mass at mu[m].
#
# prior_sd[1] must be zero, so the first component is a point mass at mu[m].
# Candidate means and prior standard deviations may be generated from y and se.
# By default the null state is exactly delta_0; set null_state = "adaptive"
# to estimate an ash mixture centered at zero for that state as well.
# After a short fixed-grid warm-up, supported non-null means can be learned by
# the same EM sufficient statistics used for rho. Voronoi constraints and
# likelihood-checked state pruning prevent dense candidates from collapsing.


#' Fit an adaptive-shrinkage hidden Markov model
#'
#' @description
#' Fits a finite-state hidden Markov model to heteroskedastic normal effect
#' estimates. Each hidden state has an adaptive-shrinkage prior formed from a
#' point mass at its state center and a grid of normal components centered at
#' the same value. Transition probabilities, state-specific mixture weights,
#' and (optionally) supported state centers are learned by generalized
#' EM/Baum--Welch updates.
#'
#' @details
#' If `mu` is `NULL`, a dense grid of candidate state centers is constructed
#' from `y` and screened before the first forward-backward pass. With
#' `nonnegative_state_means = TRUE`, state 1 is fixed at zero, all non-null
#' state centers are constrained to be strictly positive, and every continuous
#' ash component is a normal distribution truncated to `[0, Inf)`. Consequently
#' the latent effects and posterior means are nonnegative as well. Set
#' `nonnegative_state_means = FALSE` to use ordinary Gaussian ash components on
#' the real line and construct a signed symmetric center grid.
#'
#' Low-occupancy or nearly duplicated states can be proposed for pruning.
#' A proposed deletion is accepted only after recomputing the complete HMM
#' marginal likelihood. After fitting, an optional information-criterion gate
#' can replace the fitted model by the exact all-zero HMM, and a separate
#' change-penalized decoder reports contiguous steps.
#'
#' @param y Numeric vector of noisy effect estimates.
#' @param se Numeric vector of known standard errors. It must have the same
#'   length as `y`, contain only finite values, and be strictly positive.
#' @param mu Optional numeric vector of initial state centers. State 1 is the
#'   zero/hub state. When `NULL`, centers are constructed and screened
#'   automatically. Supply `mu` when using state-sized custom initial values.
#' @param prior_sd Optional nondecreasing numeric vector of ash prior standard
#'   deviations. Its first value must be zero, representing a point mass. When
#'   `NULL`, a Stephens-style geometric scale grid is constructed automatically.
#' @param half_grid Number of nonnegative candidates in the automatic mean
#'   grid, including zero. Signed mode appends their nonzero negatives.
#' @param grid_shape Positive power-grid shape parameter. Values above one put
#'   relatively more candidates near the outer part of the range.
#' @param grid_expansion Positive multiplier applied to the empirical maximum
#'   absolute observation when constructing the automatic mean grid.
#' @param grid_max_abs Optional finite positive endpoint used instead of
#'   `max(abs(y))` when constructing the automatic grid.
#' @param nonnegative_state_means Logical; if `TRUE` (default), use one state
#'   centered at zero, constrain every non-null state center to be strictly
#'   positive, and truncate all continuous ash components below at zero. Thus
#'   both prior and posterior effects have support on `[0, Inf)`. If `FALSE`,
#'   allow signed state centers and effects on the real line.
#' @param positive_state_means Deprecated compatibility alias for
#'   `nonnegative_state_means`. Leave as `NULL` in new code. If both arguments
#'   are supplied, they must agree.
#' @param effect_support Either `"auto"` (default), `"nonnegative"`, or
#'   `"real"`. `"auto"` follows `nonnegative_state_means`. An explicit value
#'   must agree with that argument. The nonnegative model uses exact
#'   truncated-normal marginal emissions and posterior moments; it does not
#'   clip an unconstrained posterior.
#' @param positive_mean_floor Strict lower bound for learned non-null centers in
#'   nonnegative mode. The default is a scale-aware numerical value.
#' @param prefilter Logical controlling independent state-grid screening before
#'   HMM fitting. The default is `TRUE` for an automatic grid and `FALSE` for a
#'   supplied grid.
#' @param min_state_count,min_state_fraction A prefilter candidate is retained
#'   when its independent soft count is at least the larger of these absolute
#'   and fractional thresholds. The zero state is always retained.
#' @param screening_prior_sd Nonnegative extra prior standard deviation used
#'   only in the independent prefilter score.
#' @param screening_block_size Positive integer block size for memory-efficient
#'   prefilter calculations.
#' @param sequence_id Vector identifying independent sequences. Adjacent equal
#'   values belong to one sequence; transitions are not counted across changes
#'   in `sequence_id`.
#' @param forward_backward_engine Forward-backward implementation:
#'   `"auto"` (default) uses the fast scaled matrix recursion for dense
#'   transition masks and the sparse log-domain recursion for sparse masks;
#'   `"scaled"` requests the fast recursion with automatic log-domain fallback;
#'   `"log"` always uses the reference log-domain recursion.
#' @param inference State-inference method. `"exact"` uses ordinary
#'   forward--backward in every EM iteration. `"variational"` uses a
#'   checkerboard coordinate-ascent approximation
#'   `q(Q[1:T]) = prod_t q_t(Q[t])` during training and, by default, performs
#'   one exact forward--backward pass for the returned posterior and marginal
#'   likelihood. Variational inference requires a full transition topology.
#' @param variational_maxiter Maximum checkerboard coordinate sweeps inside
#'   each variational E-step. The small default exploits warm starts between
#'   outer iterations; increase it when a tighter training ELBO matters more
#'   than speed.
#' @param variational_tolerance Positive convergence tolerance for the
#'   variational state ELBO.
#' @param variational_probability_floor Tiny lower bound on variational state
#'   probabilities. It prevents exactly zero transition counts.
#' @param variational_final_exact Logical; after variational training, run one
#'   exact forward--backward pass under the fitted parameters. This is strongly
#'   recommended and is always done when BIC null selection is requested.
#' @param cache_component_emissions Logical; cache the state-by-scale component
#'   log densities and refresh only states whose centers move. This can reduce
#'   density evaluations substantially at a modest memory cost.
#' @param topology Transition topology used when `transition_mask` is `NULL`:
#'   `"full"` allows all transitions; `"hub"` disallows direct transitions
#'   between distinct non-null states.
#' @param transition_mask Optional logical square matrix specifying allowed
#'   transitions between the supplied/retained states.
#' @param init_transition Optional initial row-stochastic transition matrix.
#' @param init_prob Optional initial state-probability vector.
#' @param init_rho Optional initial state-by-scale matrix of ash mixture weights.
#' @param stay_probability Initial self-transition probability used when
#'   constructing a transition matrix automatically. A scalar is recycled by
#'   state.
#' @param null_state Either `"pointmass"` (default), which fixes the null prior
#'   to the Dirac mass at zero, or `"adaptive"`, which learns an ash mixture
#'   centered at zero.
#' @param fixed_pointmass_states Optional integer indices of states whose ash
#'   mixture is fixed entirely on its point component. `NULL` derives the value
#'   from `null_state`; an explicit value is an advanced override.
#' @param shared_mixture Logical; if `TRUE`, all free states share one learned
#'   vector of ash mixture weights. The default learns each state separately.
#' @param estimate_init Logical; estimate initial state probabilities rather
#'   than keeping their initialized values fixed.
#' @param transition_prior Dirichlet parameters for allowed transition
#'   probabilities. Supply a scalar or a state-by-state matrix; values must be
#'   at least one.
#' @param mixture_prior Dirichlet parameters for ash mixture weights. Supply a
#'   scalar, a vector with one entry per scale, or a state-by-scale matrix;
#'   values must be at least one.
#' @param init_prior Dirichlet parameters for the initial distribution, used
#'   only when `estimate_init = TRUE`.
#' @param learn_state_means Logical; learn eligible non-null state centers after
#'   the fixed-grid warm-up.
#' @param fixed_mean_states Integer state indices whose centers are never moved.
#'   State 1 should normally remain fixed.
#' @param mean_update_start First EM iteration at which state-center updates and
#'   center-eligibility decisions are allowed.
#' @param mean_min_effective_count Minimum smoothed state occupancy required for
#'   a center update.
#' @param mean_min_pointmass_weight Minimum fitted weight on the state's point
#'   component required for center learning.
#' @param mean_min_self_transition Minimum self-transition probability required
#'   for center learning. This favors persistent plateau-like states.
#' @param mean_damping Number in `(0, 1]` multiplying each eligible center move.
#' @param mean_bounds Either `"voronoi"` (default), which keeps centers in
#'   nonoverlapping anchor cells, or `"none"`.
#' @param prune_states Logical; enable likelihood-checked dynamic state pruning.
#' @param prune_start First EM iteration eligible for dynamic pruning.
#' @param prune_every Positive number of EM iterations between pruning checks.
#' @param prune_min_state_count,prune_min_state_fraction A state is a
#'   low-occupancy pruning candidate below the larger absolute/fractional cutoff.
#' @param prune_max_fraction Maximum fraction of current states considered in
#'   one deletion batch. It must lie strictly between zero and one.
#' @param prune_max_loglik_loss Maximum allowed decrease in the deletion-check
#'   criterion. This is the full HMM marginal log likelihood under exact
#'   training and the mean-field state ELBO under variational training. The
#'   historical argument name is retained for compatibility.
#' @param merge_distance Distance below which adjacent fitted centers are
#'   considered near duplicates. `NULL` uses `0.05 * median(se)`.
#' @param null_model_selection Either `"bic"` for the strict all-zero safety
#'   comparison or `"none"` to retain the fitted HMM unconditionally.
#' @param null_selection_gamma Nonnegative multiplier for the additional
#'   candidate-grid term in the strict-null information criterion.
#' @param null_bic_margin Nonnegative margin favoring retention of the fitted
#'   HMM over the strict null.
#' @param parameter_tolerance Positive threshold used when counting numerically
#'   active parameters for the strict-null information criterion.
#' @param step_penalty Nonnegative penalty for each state change in the separate
#'   penalized decoder. `NULL` uses `step_penalty_scale * log(length(y))`.
#' @param step_penalty_scale Nonnegative multiplier for the default step penalty.
#'   The default `1.5` accounts for both an added segment level and an unknown
#'   change location, reducing transient one-observation steps.
#' @param maxiter Maximum number of generalized EM iterations.
#' @param tolerance Positive relative convergence tolerance for the penalized
#'   fixed-dimensional objective.
#' @param verbose Logical; print grid selection and per-iteration diagnostics.
#'
#' @return An object of class `ash_hmm_fit`, which is a list containing:
#'   \describe{
#'   \item{call}{The matched function call.}
#'   \item{state_probability}{An observation-by-state matrix of smoothed
#'     probabilities `Pr(Q[t] = m | y)`.}
#'   \item{posterior}{A list containing the marginal posterior `mean`, `sd`,
#'     `probability_ge_zero`, `probability_le_zero`, `probability_zero`, and
#'     local false sign rate `lfsr` for every observation.}
#'   \item{viterbi_state}{The ordinary joint MAP state path under the fitted
#'     transition matrix.}
#'   \item{penalized_state}{The state path from the explicit change-penalized
#'     decoder.}
#'   \item{step_selection}{A list with a segment table, per-sequence counts,
#'     total `step_count`, total `change_count`, `occupied_state_count`, decoded
#'     `state`, and the applied `penalty`.}
#'   \item{boundary_probability}{Posterior probabilities of a state change
#'     between adjacent observations; entries spanning independent sequences
#'     are `NA`.}
#'   \item{fitted}{Fitted centers and scale grid, mixture weights, transition
#'     matrix and mask, initial probabilities, state identifiers and
#'     occupancies, null-state type, effect support, constraints, and
#'     mean-learning settings.}
#'   \item{grid}{Automatic-grid settings, original and retained candidates,
#'     screening statistics, selected scale grid, and pruning history.}
#'   \item{log_likelihood,log_null,log_evidence_ratio}{The fitted HMM marginal
#'     log likelihood, exact all-zero log likelihood, and their difference.}
#'   \item{variational_bound,training_objective,inference}{The final
#'     variational training bound (or `NA` for exact training), the criterion
#'     used during training, and diagnostics recording the inference method,
#'     inner convergence, and whether a final exact smoothing pass was used.}
#'   \item{model_selection}{Details of the optional strict-null comparison,
#'     including criteria, effective dimension, and whether the result was
#'     collapsed to the exact null.}
#'   \item{history}{Per-iteration likelihood, objective, number of states,
#'     moved centers, and pruned-state count.}
#'   \item{mean_history}{Long-form history of state centers and occupancies.}
#'   \item{pruning_history}{One row per removed state, recording its persistent
#'     identifier, original grid index, centers, occupancy, and reason.}
#'   \item{converged}{Logical convergence indicator.}
#'   \item{iterations}{Number of completed generalized EM/model-reduction
#'     iterations recorded after initialization.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' truth <- c(rep(0, 50), rep(2, 50), rep(0, 50))
#' se <- rep(0.5, length(truth))
#' y <- rnorm(length(truth), truth, se)
#'
#' fit <- fit_ash_hmm(y, se, nonnegative_state_means = TRUE)
#' fit$posterior$mean
#' fit$step_selection$segments
#'
#' # Use a signed grid when negative state centers are scientifically possible.
#' signed_fit <- fit_ash_hmm(y - 1, se, nonnegative_state_means = FALSE)
#' }
#'
#' @export
fit_ash_hmm <- function(y, se, mu = NULL, prior_sd = NULL,
                        half_grid = 20L,
                        grid_shape = 3,
                        grid_expansion = 3,
                        grid_max_abs = NULL,
                        nonnegative_state_means = FALSE,
                        positive_state_means = NULL,
                        effect_support = c("auto", "nonnegative", "real"),
                        positive_mean_floor = NULL,
                        prefilter = NULL,
                        min_state_count = 0.5,
                        min_state_fraction = 1e-5,
                        screening_prior_sd = 0,
                        screening_block_size = 5000L,
                        sequence_id = rep(1L, length(y)),
                        forward_backward_engine = c("auto", "scaled", "log"),
                        inference = c("exact", "variational"),
                        variational_maxiter = 1L,
                        variational_tolerance = 1e-5,
                        variational_probability_floor = 1e-12,
                        variational_final_exact = TRUE,
                        cache_component_emissions = TRUE,
                        topology = c( "full" ,"hub"   ),
                        transition_mask = NULL,
                        init_transition = NULL,
                        init_prob = NULL,
                        init_rho = NULL,
                        stay_probability = 0.95,
                        null_state = c("pointmass", "adaptive"),
                        fixed_pointmass_states = NULL,
                        shared_mixture = FALSE,
                        estimate_init = FALSE,
                        transition_prior = 1,
                        mixture_prior = 1,
                        init_prior = 1,
                        learn_state_means = TRUE,
                        fixed_mean_states = 1L,
                        mean_update_start = 3L,
                        mean_min_effective_count = 2,
                        mean_min_pointmass_weight = 0,
                        mean_min_self_transition = 0.93,
                        mean_damping = 1,
                        mean_bounds = c("voronoi", "none"),
                        prune_states = TRUE,
                        prune_start = 5L,
                        prune_every = 5L,
                        prune_min_state_count = 1,
                        prune_min_state_fraction = 1e-4,
                        prune_max_fraction = 0.25,
                        prune_max_loglik_loss = 0.05,
                        merge_distance = NULL,
                        null_model_selection = c("none", "bic"  ),
                        null_selection_gamma = 1,
                        null_bic_margin = 0,
                        parameter_tolerance = 1e-6,
                        step_penalty = NULL,
                        step_penalty_scale = 1.5,
                        maxiter = 20L,
                        tolerance = 1e-5,
                        verbose = FALSE) {
  call <- match.call()
  nonnegative_was_missing <- missing(nonnegative_state_means)
  forward_backward_engine <- match.arg(forward_backward_engine)
  inference <- match.arg(inference)
  topology <- match.arg(topology)
  mean_bounds <- match.arg(mean_bounds)
  null_model_selection <- match.arg(null_model_selection)
  null_state <- match.arg(null_state)
  effect_support <- match.arg(effect_support)
  .ash_hmm_check_data(y, se, sequence_id)
  if (!is.null(positive_state_means)) {
    if (!is.logical(positive_state_means) ||
        length(positive_state_means) != 1L || is.na(positive_state_means)) {
      stop("'positive_state_means' must be NULL, TRUE, or FALSE.",
           call. = FALSE)
    }
    if (!nonnegative_was_missing &&
        !identical(nonnegative_state_means, positive_state_means)) {
      stop(paste("'nonnegative_state_means' and its compatibility alias",
                 "'positive_state_means' must agree when both are supplied."),
           call. = FALSE)
    }
    nonnegative_state_means <- positive_state_means
  }
  if (!is.logical(nonnegative_state_means) ||
      length(nonnegative_state_means) != 1L || is.na(nonnegative_state_means)) {
    stop("'nonnegative_state_means' must be TRUE or FALSE.", call. = FALSE)
  }
  expected_support <- if (nonnegative_state_means) "nonnegative" else "real"
  if (effect_support == "auto") {
    effect_support <- expected_support
  } else if (effect_support != expected_support) {
    stop(paste0("'effect_support = \"", effect_support,
                "\"' conflicts with 'nonnegative_state_means = ",
                nonnegative_state_means, "'."), call. = FALSE)
  }

  automatic_mu <- is.null(mu)
  if (is.null(prefilter)) prefilter <- automatic_mu
  if (!is.logical(prefilter) || length(prefilter) != 1L || is.na(prefilter)) {
    stop("'prefilter' must be TRUE or FALSE.", call. = FALSE)
  }
  if (automatic_mu && (!is.null(transition_mask) || !is.null(init_transition) ||
                       !is.null(init_prob) || !is.null(init_rho))) {
    stop(paste("Supply 'mu' explicitly when using a custom transition mask,",
               "initial transition matrix, initial probabilities, or ash weights."),
         call. = FALSE)
  }

  if (automatic_mu) {
    if (prefilter) {
      grid_selection <- ash_hmm_select_grid(
        y, se,
        half_grid = half_grid,
        shape = grid_shape,
        expansion = grid_expansion,
        max_abs = grid_max_abs,
        min_state_count = min_state_count,
        min_state_fraction = min_state_fraction,
        screening_prior_sd = screening_prior_sd,
        block_size = screening_block_size,
        positive_only = nonnegative_state_means)
      mu <- grid_selection$selected_mu
    } else {
      mu <- ash_mu_grid(y, half_grid = half_grid, shape = grid_shape,
                        expansion = grid_expansion, max_abs = grid_max_abs,
                        positive_only = nonnegative_state_means)
      grid_selection <- list(
        automatic = TRUE,
        full_mu = mu,
        selected_mu = mu,
        selected_index = seq_along(mu),
        removed_index = integer(),
        effective_count = rep(NA_real_, length(mu)),
        mean_weight = rep(NA_real_, length(mu)),
        cutoff = NA_real_,
        settings = list(prefilter = FALSE))
    }
  } else {
    if (!is.numeric(mu) || !length(mu) || anyNA(mu) || any(!is.finite(mu))) {
      stop("'mu' must be a finite numeric vector.", call. = FALSE)
    }
    if (prefilter) {
      if (!is.null(transition_mask) || !is.null(init_transition) ||
          !is.null(init_prob) || !is.null(init_rho)) {
        stop("Prefilter a supplied grid before providing state-sized initial values.",
             call. = FALSE)
      }
      grid_selection <- .ash_hmm_prefilter_mu(
        y, se, mu,
        min_state_count = min_state_count,
        min_state_fraction = min_state_fraction,
        screening_prior_sd = screening_prior_sd,
        block_size = screening_block_size)
      grid_selection$automatic <- FALSE
      mu <- grid_selection$selected_mu
    } else {
      grid_selection <- list(
        automatic = FALSE,
        full_mu = mu,
        selected_mu = mu,
        selected_index = seq_along(mu),
        removed_index = integer(),
        effective_count = rep(NA_real_, length(mu)),
        mean_weight = rep(NA_real_, length(mu)),
        cutoff = NA_real_,
        settings = list(prefilter = FALSE))
    }
  }
  grid_selection$settings$nonnegative_state_means <-
    nonnegative_state_means
  grid_selection$settings$effect_support <- effect_support

  if (nonnegative_state_means) {
    if (mu[1L] != 0 || (length(mu) > 1L && any(mu[-1L] <= 0))) {
      stop(paste("With 'nonnegative_state_means = TRUE', mu[1] must equal zero",
                 "and every non-null state center must be strictly positive."),
           call. = FALSE)
    }
  } else if (mu[1L] != 0) {
    warning("The hub state is state 1, but mu[1] is not zero.")
  }
  automatic_prior_sd <- is.null(prior_sd)
  if (automatic_prior_sd) prior_sd <- ash_sigma_grid(y, se)
  if (!is.numeric(prior_sd) || !length(prior_sd) || anyNA(prior_sd) ||
      any(!is.finite(prior_sd)) || any(prior_sd < 0) || prior_sd[1L] != 0) {
    stop("'prior_sd' must be finite, nonnegative, and start with zero.", call. = FALSE)
  }
  if (is.unsorted(prior_sd, strictly = FALSE)) {
    stop("'prior_sd' must be sorted in nondecreasing order.", call. = FALSE)
  }

  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  maxiter <- as.integer(maxiter)
  variational_maxiter <- as.integer(variational_maxiter)
  if (maxiter < 1L || tolerance <= 0 ||
      length(variational_maxiter) != 1L ||
      is.na(variational_maxiter) || variational_maxiter < 1L ||
      !is.finite(variational_tolerance) ||
      variational_tolerance <= 0 ||
      !is.finite(variational_probability_floor) ||
      variational_probability_floor < 0 ||
      variational_probability_floor >= 1 / m) {
    stop(paste(
      "Require maxiter and variational_maxiter >= 1; positive tolerances;",
      "and 0 <= variational_probability_floor < 1 / number_of_states."),
      call. = FALSE)
  }
  if (!is.logical(variational_final_exact) ||
      length(variational_final_exact) != 1L ||
      is.na(variational_final_exact) ||
      !is.logical(cache_component_emissions) ||
      length(cache_component_emissions) != 1L ||
      is.na(cache_component_emissions)) {
    stop(paste(
      "'variational_final_exact' and 'cache_component_emissions' must be",
      "TRUE or FALSE."), call. = FALSE)
  }
  if (!is.logical(learn_state_means) || length(learn_state_means) != 1L ||
      is.na(learn_state_means) || !is.logical(prune_states) ||
      length(prune_states) != 1L || is.na(prune_states)) {
    stop("'learn_state_means' and 'prune_states' must be TRUE or FALSE.",
         call. = FALSE)
  }
  mean_update_start <- as.integer(mean_update_start)
  prune_start <- as.integer(prune_start)
  prune_every <- as.integer(prune_every)
  if (mean_update_start < 1L || prune_start < 1L || prune_every < 1L) {
    stop("Mean/pruning iteration controls must be positive integers.",
         call. = FALSE)
  }
  scalar_nonnegative <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= 0
  }
  if (is.null(positive_mean_floor)) {
    positive_mean_floor <- sqrt(.Machine$double.eps) *
      max(1, stats::median(se))
  }
  if (!scalar_nonnegative(mean_min_effective_count) ||
      !scalar_nonnegative(mean_min_pointmass_weight) ||
      mean_min_pointmass_weight > 1 ||
      !scalar_nonnegative(mean_min_self_transition) ||
      mean_min_self_transition > 1 ||
      !is.numeric(mean_damping) || length(mean_damping) != 1L ||
      !is.finite(mean_damping) || mean_damping <= 0 || mean_damping > 1 ||
      !scalar_nonnegative(prune_min_state_count) ||
      !scalar_nonnegative(prune_min_state_fraction) ||
      prune_min_state_fraction > 1 ||
      !is.numeric(prune_max_fraction) || length(prune_max_fraction) != 1L ||
      !is.finite(prune_max_fraction) || prune_max_fraction <= 0 ||
      prune_max_fraction >= 1 ||
      !scalar_nonnegative(prune_max_loglik_loss) ||
      !scalar_nonnegative(positive_mean_floor) ||
      !scalar_nonnegative(null_selection_gamma) ||
      !scalar_nonnegative(null_bic_margin) ||
      !scalar_nonnegative(parameter_tolerance) || parameter_tolerance == 0 ||
      !scalar_nonnegative(step_penalty_scale)) {
    stop("Invalid state-mean learning or pruning control.", call. = FALSE)
  }
  if (nonnegative_state_means && positive_mean_floor <= 0) {
    stop("'positive_mean_floor' must be strictly positive.", call. = FALSE)
  }
  if (is.null(step_penalty)) {
    step_penalty <- step_penalty_scale * log(max(2, n))
  }
  if (!scalar_nonnegative(step_penalty)) {
    stop("'step_penalty' must be NULL or a finite nonnegative scalar.",
         call. = FALSE)
  }
  if (is.null(merge_distance)) merge_distance <- 0.05 * stats::median(se)
  if (!scalar_nonnegative(merge_distance)) {
    stop("'merge_distance' must be NULL or a finite nonnegative scalar.",
         call. = FALSE)
  }
  if (verbose && automatic_mu) {
    message(sprintf("automatic mean grid: retained %d of %d states before HMM fitting",
                    length(mu), length(grid_selection$full_mu)))
  }
  if (is.null(fixed_pointmass_states)) {
    fixed_pointmass_states <- if (null_state == "pointmass") 1L else integer()
  }
  fixed_pointmass_states <- sort(unique(as.integer(fixed_pointmass_states)))
  if (any(fixed_pointmass_states < 1L | fixed_pointmass_states > m)) {
    stop("'fixed_pointmass_states' contains an invalid state index.", call. = FALSE)
  }
  fixed_mean_states <- sort(unique(as.integer(fixed_mean_states)))
  if (any(fixed_mean_states < 1L | fixed_mean_states > m)) {
    stop("'fixed_mean_states' contains an invalid state index.", call. = FALSE)
  }
  if (learn_state_means && mean_bounds == "voronoi" && anyDuplicated(mu)) {
    stop("Initial state means must be distinct with Voronoi mean bounds.",
         call. = FALSE)
  }
  mean_anchor <- mu
  state_id <- seq_len(m)
  grid_index <- grid_selection$selected_index
  pruned_grid_index <- integer()
  pruning_history <- data.frame(
    iteration = integer(), state_id = integer(), grid_index = integer(),
    anchor_mu = numeric(), fitted_mu = numeric(), occupancy = numeric(),
    reason = character(), stringsAsFactors = FALSE)

  if (is.null(transition_mask)) {
    transition_mask <- ash_hmm_transition_mask(m, topology)
  } else {
    if (!is.matrix(transition_mask) || any(dim(transition_mask) != c(m, m))) {
      stop("'transition_mask' must be an M by M matrix.", call. = FALSE)
    }
    transition_mask <- transition_mask != 0
    if (any(rowSums(transition_mask) == 0)) {
      stop("Every state needs at least one allowed outgoing transition.", call. = FALSE)
    }
  }
  if (inference == "variational" && !all(transition_mask)) {
    stop(paste(
      "Variational state inference requires topology='full' and a full",
      "transition mask. Use inference='exact' for structural zeros."),
      call. = FALSE)
  }

  if (is.null(init_transition)) {
    A <- .ash_hmm_initial_transition(transition_mask, stay_probability)
  } else {
    A <- .ash_hmm_validate_transition(init_transition, transition_mask)
  }

  if (is.null(init_prob)) init_prob <- .ash_hmm_stationary(A)
  if (!is.numeric(init_prob) || length(init_prob) != m || anyNA(init_prob) ||
      any(!is.finite(init_prob)) || any(init_prob < 0) || sum(init_prob) <= 0) {
    stop("'init_prob' must be a nonnegative length-M vector with positive mass.", call. = FALSE)
  }
  init_prob <- .ash_hmm_normalize(init_prob)

  if (is.null(init_rho)) {
    rho <- matrix(1 / l, m, l)
  } else {
    if (!is.matrix(init_rho) || any(dim(init_rho) != c(m, l)) ||
        anyNA(init_rho) || any(!is.finite(init_rho)) || any(init_rho < 0) ||
        any(rowSums(init_rho) <= 0)) {
      stop("'init_rho' must be a nonnegative M by L matrix with positive row sums.",
           call. = FALSE)
    }
    rho <- init_rho / rowSums(init_rho)
  }
  if (length(fixed_pointmass_states)) {
    rho[fixed_pointmass_states, ] <- 0
    rho[fixed_pointmass_states, 1L] <- 1
  }

  transition_prior <- .ash_hmm_expand_prior(transition_prior, m, m,
                                            "transition_prior")
  mixture_prior <- .ash_hmm_expand_prior(mixture_prior, m, l, "mixture_prior")
  init_prior <- rep(init_prior, length.out = m)
  if (anyNA(init_prior) || any(!is.finite(init_prior)) || any(init_prior < 1)) {
    stop("Every 'init_prior' entry must be finite and at least one.", call. = FALSE)
  }

  mixture_active <- matrix(TRUE, m, l)
  if (length(fixed_pointmass_states)) {
    mixture_active[fixed_pointmass_states, ] <- FALSE
  }

  objective_value <- function(log_likelihood, A, rho, init_prob,
                              transition_prior_current = transition_prior,
                              transition_mask_current = transition_mask,
                              mixture_prior_current = mixture_prior,
                              mixture_active_current = mixture_active,
                              init_prior_current = init_prior) {
    ans <- log_likelihood +
      .ash_hmm_dirichlet_penalty(
        A, transition_prior_current, transition_mask_current) +
      .ash_hmm_dirichlet_penalty(
        rho, mixture_prior_current, mixture_active_current)
    if (estimate_init) {
      use <- init_prior_current > 1
      if (any(use)) {
        if (any(init_prob[use] <= 0)) return(-Inf)
        ans <- ans + sum((init_prior_current[use] - 1) * log(init_prob[use]))
      }
    }
    ans
  }

  infer_states <- function(
    log_emission_current, A_current, init_prob_current,
    initial_probability = NULL,
    mask_current = transition_mask) {
    if (inference == "exact") {
      return(.ash_hmm_forward_backward(
        log_emission_current, A_current, init_prob_current,
        mask_current, sequence_id,
        engine = forward_backward_engine))
    }
    .ash_hmm_variational_states(
      log_emission_current, A_current, init_prob_current,
      mask_current, sequence_id,
      initial_probability = initial_probability,
      maxiter = variational_maxiter,
      tolerance = variational_tolerance,
      probability_floor = variational_probability_floor)
  }

  component_log_emission <- NULL
  if (cache_component_emissions) {
    component_log_emission <-
      .ash_hmm_component_log_emission_array(
        y, se, mu, prior_sd, effect_support)
    log_emission <- .ash_hmm_mix_component_cache(
      component_log_emission, rho)
  } else {
    log_emission <- .ash_hmm_log_emission(
      y, se, mu, prior_sd, rho, effect_support)
  }
  fb <- infer_states(log_emission, A, init_prob)
  objective <- objective_value(fb$log_likelihood, A, rho, init_prob)
  history <- data.frame(iteration = 0L,
                        log_likelihood = fb$log_likelihood,
                        objective = objective,
                        criterion = if (inference == "exact") {
                          "log_likelihood"
                        } else "variational_bound",
                        states = m,
                        means_updated = 0L,
                        states_pruned = 0L)
  mean_history <- list(data.frame(
    iteration = 0L, state_id = state_id, grid_index = grid_index,
    anchor_mu = mean_anchor, mu = mu,
    occupancy = colSums(fb$state_probability)))
  mean_state_enabled <- rep(FALSE, m)
  converged <- FALSE

  for (iter in seq_len(maxiter)) {
    free_states <- setdiff(seq_len(m), fixed_pointmass_states)
    statistics <- if (cache_component_emissions) {
      .ash_hmm_component_statistics_cached(
        y, se, mu, prior_sd, rho, log_emission,
        fb$state_probability, component_log_emission,
        effect_support)
    } else {
      .ash_hmm_component_statistics(
        y, se, mu, prior_sd, rho, log_emission,
        fb$state_probability, effect_support)
    }
    component_counts <- statistics$counts

    rho_new <- rho
    if (length(free_states)) {
      if (shared_mixture) {
        numerator <- colSums(component_counts[free_states, , drop = FALSE]) +
          colSums(mixture_prior[free_states, , drop = FALSE] - 1)
        common <- .ash_hmm_normalize(
          numerator, fallback = colMeans(rho[free_states, , drop = FALSE]))
        rho_new[free_states, ] <- matrix(
          rep(common, each = length(free_states)), length(free_states), l)
      } else {
        for (state in free_states) {
          numerator <- component_counts[state, ] + mixture_prior[state, ] - 1
          rho_new[state, ] <- .ash_hmm_normalize(
            numerator, fallback = rho[state, ])
        }
      }
    }
    if (length(fixed_pointmass_states)) {
      rho_new[fixed_pointmass_states, ] <- 0
      rho_new[fixed_pointmass_states, 1L] <- 1
    }
    A_new <- matrix(0, m, m)
    for (state in seq_len(m)) {
      allowed <- which(transition_mask[state, ])
      numerator <- fb$transition_counts[state, allowed] +
        transition_prior[state, allowed] - 1
      A_new[state, allowed] <- .ash_hmm_normalize(
        numerator, fallback = A[state, allowed])
    }

    # Long, persistent states are the ones for which estimating a precise
    # plateau center is useful. The self-transition gate avoids moving centers
    # that merely tile a smooth bump or are supported by scattered positions.
    current_mean_gate <- rho_new[, 1L] >= mean_min_pointmass_weight &
      diag(A_new) >= mean_min_self_transition
    if (iter == mean_update_start) {
      mean_state_enabled <- current_mean_gate
    } else if (iter > mean_update_start) {
      # Eligibility is established after the warm-up and can switch off, but
      # not on. This stops a weak state from becoming movable merely because
      # later pruning artificially raises its self-transition probability.
      mean_state_enabled <- mean_state_enabled & current_mean_gate
    }
    persistence_damping <- pmin(1, pmax(
      0, (diag(A_new) - mean_min_self_transition) /
        max(1e-8, 1 - mean_min_self_transition)))
    mean_lower_bound <- if (nonnegative_state_means) {
      c(0, rep(positive_mean_floor, max(0L, m - 1L)))
    } else rep(-Inf, m)
    if (effect_support == "nonnegative") {
      mean_update <- .ash_hmm_update_truncated_means(
        y = y, se = se, mu = mu, anchor = mean_anchor,
        prior_sd = prior_sd, rho = rho, log_emission = log_emission,
        gamma = fb$state_probability, statistics = statistics,
        learn = learn_state_means && iter >= mean_update_start,
        fixed = fixed_mean_states,
        minimum_count = mean_min_effective_count,
        eligible = mean_state_enabled,
        damping = mean_damping * persistence_damping,
        minimum = mean_lower_bound,
        bounds = mean_bounds)
    } else {
      mean_update <- .ash_hmm_update_means(
        mu = mu, anchor = mean_anchor, statistics = statistics,
        learn = learn_state_means && iter >= mean_update_start,
        fixed = fixed_mean_states,
        minimum_count = mean_min_effective_count,
        eligible = mean_state_enabled,
        damping = mean_damping * persistence_damping,
        minimum = mean_lower_bound,
        bounds = mean_bounds)
    }
    mu_new <- mean_update$mu

    init_prob_new <- init_prob
    if (estimate_init) {
      numerator <- fb$initial_counts + init_prior - 1
      init_prob_new <- .ash_hmm_normalize(numerator, fallback = init_prob)
    }

    component_log_emission_new <- component_log_emission
    if (cache_component_emissions) {
      component_log_emission_new <-
        .ash_hmm_refresh_component_cache(
          component_log_emission, y, se, mu, mu_new,
          prior_sd, effect_support)
      log_emission_new <- .ash_hmm_mix_component_cache(
        component_log_emission_new, rho_new)
    } else {
      log_emission_new <- .ash_hmm_log_emission(
        y, se, mu_new, prior_sd, rho_new, effect_support)
    }
    fb_new <- infer_states(
      log_emission_new, A_new, init_prob_new,
      initial_probability = fb$state_probability)
    objective_new <- objective_value(
      fb_new$log_likelihood, A_new, rho_new, init_prob_new)

    if (objective_new < objective - 1e-7 * (1 + abs(objective))) {
      warning("The fixed-dimension EM objective decreased beyond numerical tolerance.")
    }
    delta_em <- objective_new - objective
    states_pruned_iter <- 0L

    # Dynamic deletion is deliberately outside the fixed-dimensional EM step.
    # Candidate states must be empty or nearly duplicate. Exact training checks
    # the full HMM marginal likelihood; variational training checks its current
    # mean-field state bound. Groups that fail the selected training-criterion
    # tolerance are halved.
    pruning_due <- prune_states && iter >= prune_start &&
      ((iter - prune_start) %% prune_every == 0L) && m > 1L
    if (pruning_due) {
      occupancy_new <- colSums(fb_new$state_probability)
      protected <- sort(unique(c(1L, fixed_pointmass_states,
                                 fixed_mean_states)))
      cutoff <- max(prune_min_state_count, n * prune_min_state_fraction)
      low_candidates <- setdiff(which(occupancy_new < cutoff), protected)
      close_candidates <- .ash_hmm_close_state_candidates(
        mu_new, occupancy_new, merge_distance, protected)
      candidates <- unique(c(low_candidates, close_candidates))
      if (length(candidates)) {
        candidates <- candidates[order(occupancy_new[candidates])]
        maximum_drop <- max(1L, floor(prune_max_fraction * m))
        maximum_drop <- min(maximum_drop, length(candidates),
                            m - length(unique(protected)))
        number_to_try <- maximum_drop
        accepted <- NULL

        while (number_to_try >= 1L && is.null(accepted)) {
          drop <- candidates[seq_len(number_to_try)]
          keep <- setdiff(seq_len(m), drop)
          mask_try <- transition_mask[keep, keep, drop = FALSE]
          A_try <- A_new[keep, keep, drop = FALSE]
          valid_rows <- rowSums(mask_try) > 0 & rowSums(A_try) > 0
          if (all(valid_rows)) {
            A_try <- A_try / rowSums(A_try)
            init_try <- .ash_hmm_normalize(init_prob_new[keep])
            rho_try <- rho_new[keep, , drop = FALSE]
            mu_try <- mu_new[keep]
            transition_prior_try <-
              transition_prior[keep, keep, drop = FALSE]
            mixture_prior_try <- mixture_prior[keep, , drop = FALSE]
            init_prior_try <- init_prior[keep]
            fixed_pointmass_try <- match(fixed_pointmass_states, keep,
                                         nomatch = 0L)
            fixed_pointmass_try <- fixed_pointmass_try[fixed_pointmass_try > 0L]
            fixed_mean_try <- match(fixed_mean_states, keep, nomatch = 0L)
            fixed_mean_try <- fixed_mean_try[fixed_mean_try > 0L]
            mixture_active_try <- matrix(TRUE, length(keep), l)
            if (length(fixed_pointmass_try)) {
              mixture_active_try[fixed_pointmass_try, ] <- FALSE
            }
            component_log_emission_try <- NULL
            if (cache_component_emissions) {
              component_log_emission_try <-
                component_log_emission_new[, , keep, drop = FALSE]
              log_emission_try <- .ash_hmm_mix_component_cache(
                component_log_emission_try, rho_try)
            } else {
              log_emission_try <- .ash_hmm_log_emission(
                y, se, mu_try, prior_sd, rho_try, effect_support)
            }
            initial_probability_try <-
              fb_new$state_probability[, keep, drop = FALSE]
            initial_probability_try <- pmax(
              initial_probability_try,
              variational_probability_floor)
            dim(initial_probability_try) <- c(n, length(keep))
            initial_probability_try <- initial_probability_try /
              rowSums(initial_probability_try)
            fb_try <- infer_states(
              log_emission_try, A_try, init_try,
              initial_probability = initial_probability_try,
              mask_current = mask_try)
            loss <- fb_new$log_likelihood - fb_try$log_likelihood
            if (is.finite(loss) && loss <= prune_max_loglik_loss) {
              objective_try <- objective_value(
                fb_try$log_likelihood, A_try, rho_try, init_try,
                transition_prior_try, mask_try, mixture_prior_try,
                mixture_active_try, init_prior_try)
              accepted <- list(
                drop = drop, keep = keep, A = A_try, rho = rho_try,
                mu = mu_try, init_prob = init_try, mask = mask_try,
                transition_prior = transition_prior_try,
                mixture_prior = mixture_prior_try,
                init_prior = init_prior_try,
                mixture_active = mixture_active_try,
                fixed_pointmass_states = fixed_pointmass_try,
                fixed_mean_states = fixed_mean_try,
                component_log_emission = component_log_emission_try,
                log_emission = log_emission_try, fb = fb_try,
                objective = objective_try)
            }
          }
          number_to_try <- floor(number_to_try / 2L)
        }

        if (!is.null(accepted)) {
          removed <- accepted$drop
          reason <- ifelse(
            removed %in% low_candidates & removed %in% close_candidates,
            "low occupancy and near duplicate",
            ifelse(removed %in% close_candidates,
                   "near duplicate", "low occupancy"))
          pruning_history <- rbind(
            pruning_history,
            data.frame(iteration = rep(iter, length(removed)),
                       state_id = state_id[removed],
                       grid_index = grid_index[removed],
                       anchor_mu = mean_anchor[removed],
                       fitted_mu = mu_new[removed],
                       occupancy = occupancy_new[removed],
                       reason = reason, stringsAsFactors = FALSE))
          pruned_grid_index <- c(pruned_grid_index, grid_index[removed])
          state_id <- state_id[accepted$keep]
          grid_index <- grid_index[accepted$keep]
          mean_anchor <- mean_anchor[accepted$keep]
          mean_state_enabled <- mean_state_enabled[accepted$keep]
          A_new <- accepted$A
          rho_new <- accepted$rho
          mu_new <- accepted$mu
          init_prob_new <- accepted$init_prob
          transition_mask <- accepted$mask
          transition_prior <- accepted$transition_prior
          mixture_prior <- accepted$mixture_prior
          init_prior <- accepted$init_prior
          mixture_active <- accepted$mixture_active
          fixed_pointmass_states <- accepted$fixed_pointmass_states
          fixed_mean_states <- accepted$fixed_mean_states
          component_log_emission_new <-
            accepted$component_log_emission
          log_emission_new <- accepted$log_emission
          fb_new <- accepted$fb
          objective_new <- accepted$objective
          states_pruned_iter <- length(removed)
          m <- length(mu_new)
        }
      }
    }

    history <- rbind(
      history,
      data.frame(iteration = iter,
                 log_likelihood = fb_new$log_likelihood,
                 objective = objective_new,
                 criterion = if (inference == "exact") {
                   "log_likelihood"
                 } else "variational_bound",
                 states = m,
                 means_updated = sum(mean_update$updated),
                 states_pruned = states_pruned_iter))
    mean_history[[length(mean_history) + 1L]] <- data.frame(
      iteration = iter, state_id = state_id, grid_index = grid_index,
      anchor_mu = mean_anchor, mu = mu_new,
      occupancy = colSums(fb_new$state_probability))

    if (verbose) {
      message(sprintf(
        paste0("iteration %d: %s = %.10f, objective = %.10f, ",
               "states = %d, means moved = %d, pruned = %d"),
        iter,
        if (inference == "exact") "logLik" else "ELBO",
        fb_new$log_likelihood, objective_new, m,
        sum(mean_update$updated), states_pruned_iter))
    }

    A <- A_new
    rho <- rho_new
    mu <- mu_new
    init_prob <- init_prob_new
    log_emission <- log_emission_new
    component_log_emission <- component_log_emission_new
    fb <- fb_new
    objective <- objective_new

    minimum_iterations <- max(
      if (learn_state_means) mean_update_start else 1L,
      if (prune_states) prune_start else 1L)
    pruning_checkpoint_complete <- !prune_states || pruning_due
    if (iter >= minimum_iterations && pruning_checkpoint_complete &&
        states_pruned_iter == 0L &&
        abs(delta_em) <= tolerance * (1 + abs(objective))) {
      converged <- TRUE
      break
    }
  }

  training_fb <- fb
  training_objective <- objective
  final_exact_used <- inference == "exact" ||
    variational_final_exact || null_model_selection == "bic"
  if (inference == "variational" && final_exact_used) {
    fb <- .ash_hmm_forward_backward(
      log_emission, A, init_prob, transition_mask, sequence_id,
      engine = forward_backward_engine)
    objective <- objective_value(
      fb$log_likelihood, A, rho, init_prob)
  }
  model_log_likelihood <- if (final_exact_used) {
    fb$log_likelihood
  } else {
    NA_real_
  }

  log_null <- sum(stats::dnorm(y, mean = 0, sd = se, log = TRUE))
  initial_selected_mu <- grid_selection$selected_mu
  initial_selected_index <- grid_selection$selected_index
  learned_mean_states <- setdiff(which(mean_state_enabled), fixed_mean_states)
  effective_parameters <- .ash_hmm_effective_parameter_count(
    A, transition_mask, rho, fixed_pointmass_states,
    learned_mean_states, estimate_init, init_prob, parameter_tolerance)
  candidate_count <- max(1L, length(grid_selection$full_mu) - 1L)
  full_bic <- -2 * model_log_likelihood +
    effective_parameters * log(max(2, n)) +
    2 * null_selection_gamma * log(candidate_count)
  null_bic <- -2 * log_null
  collapse_to_null <- null_model_selection == "bic" &&
    null_bic + null_bic_margin <= full_bic
  model_selection <- list(
    method = null_model_selection,
    selected = if (collapse_to_null) {
      "strict_null"
    } else if (effect_support == "nonnegative") {
      "one_sided_ash_hmm"
    } else {
      "signed_ash_hmm"
    },
    collapsed_to_null = collapse_to_null,
    full_log_likelihood = model_log_likelihood,
    training_criterion = training_fb$log_likelihood,
    training_criterion_type = if (inference == "exact") {
      "log_likelihood"
    } else "mean_field_elbo",
    null_log_likelihood = log_null,
    effective_parameters = effective_parameters,
    candidate_count = candidate_count,
    full_bic = full_bic,
    null_bic = null_bic,
    grid_penalty_gamma = null_selection_gamma,
    margin = null_bic_margin)

  if (collapse_to_null) {
    removed <- if (m > 1L) 2L:m else integer()
    occupancy_before_null <- colSums(fb$state_probability)
    pruning_history <- rbind(
      pruning_history,
      data.frame(iteration = rep(max(history$iteration) + 1L, length(removed)),
                 state_id = state_id[removed],
                 grid_index = grid_index[removed],
                 anchor_mu = mean_anchor[removed],
                 fitted_mu = mu[removed],
                 occupancy = occupancy_before_null[removed],
                 reason = rep("global BIC selected strict null", length(removed)),
                 stringsAsFactors = FALSE))
    pruned_grid_index <- c(pruned_grid_index, grid_index[removed])

    state_id <- state_id[1L]
    grid_index <- grid_index[1L]
    mean_anchor <- 0
    mu <- 0
    m <- 1L
    rho <- matrix(0, 1L, l)
    rho[1L, 1L] <- 1
    A <- matrix(1, 1L, 1L)
    transition_mask <- matrix(TRUE, 1L, 1L)
    init_prob <- 1
    transition_prior <- matrix(transition_prior[1L, 1L], 1L, 1L)
    mixture_prior <- matrix(mixture_prior[1L, ], 1L, l)
    init_prior <- init_prior[1L]
    mixture_active <- matrix(FALSE, 1L, l)
    fixed_pointmass_states <- 1L
    fixed_mean_states <- 1L
    mean_state_enabled <- FALSE
    log_emission <- .ash_hmm_log_emission(
      y, se, mu, prior_sd, rho, effect_support)
    fb <- .ash_hmm_forward_backward(
      log_emission, A, init_prob, transition_mask, sequence_id,
      engine = forward_backward_engine)
    objective <- objective_value(
      fb$log_likelihood, A, rho, init_prob,
      transition_prior, transition_mask, mixture_prior,
      mixture_active, init_prior)
    history <- rbind(
      history,
      data.frame(iteration = max(history$iteration) + 1L,
                 log_likelihood = fb$log_likelihood,
                 objective = objective,
                 criterion = "log_likelihood",
                 states = 1L,
                 means_updated = 0L,
                 states_pruned = length(removed)))
    mean_history[[length(mean_history) + 1L]] <- data.frame(
      iteration = max(history$iteration), state_id = state_id,
      grid_index = grid_index, anchor_mu = 0, mu = 0,
      occupancy = n)
    converged <- TRUE
    model_log_likelihood <- fb$log_likelihood
  }

  posterior <- .ash_hmm_posterior_summary(
    y, se, mu, prior_sd, rho, log_emission, fb$state_probability,
    effect_support)
  viterbi_state <- .ash_hmm_viterbi(log_emission, A, init_prob,
                                    transition_mask, sequence_id)
  penalized_state <- .ash_hmm_penalized_viterbi(
    log_emission, transition_mask, sequence_id, step_penalty)
  step_summary <- .ash_hmm_step_summary(
    penalized_state, sequence_id, mu)
  step_summary$state <- penalized_state
  step_summary$penalty <- step_penalty
  grid_selection$initial_selected_mu <- initial_selected_mu
  grid_selection$initial_selected_index <- initial_selected_index
  grid_selection$selected_mu <- mu
  grid_selection$selected_index <- grid_index
  grid_selection$removed_index <- sort(unique(c(
    grid_selection$removed_index, pruned_grid_index)))
  grid_selection$pruning_history <- pruning_history

  ans <- list(call = call,
              state_probability = fb$state_probability,
              posterior = posterior,
              viterbi_state = viterbi_state,
              penalized_state = penalized_state,
              step_selection = step_summary,
              boundary_probability = fb$boundary_probability,
              fitted = list(mu = mu,
                            prior_sd = prior_sd,
                            mixture_weight = rho,
                            transition = A,
                            transition_mask = transition_mask,
                            forward_backward_engine =
                              forward_backward_engine,
                            inference = inference,
                            variational_final_exact =
                              variational_final_exact,
                            cache_component_emissions =
                              cache_component_emissions,
                            init_prob = init_prob,
                            shared_mixture = shared_mixture,
                            null_state = if (1L %in% fixed_pointmass_states) {
                              "pointmass"
                            } else "adaptive",
                            fixed_pointmass_states = fixed_pointmass_states,
                            fixed_mean_states = fixed_mean_states,
                            state_id = state_id,
                            mean_anchor = mean_anchor,
                            state_occupancy = colSums(fb$state_probability),
                            effect_support = effect_support,
                            learn_state_means = learn_state_means,
                            mean_bounds = mean_bounds,
                            mean_state_enabled = mean_state_enabled,
                            mean_update_start = mean_update_start,
                            mean_min_effective_count = mean_min_effective_count,
                            mean_min_pointmass_weight = mean_min_pointmass_weight,
                            mean_min_self_transition = mean_min_self_transition,
                            mean_damping = mean_damping,
                            nonnegative_state_means = nonnegative_state_means,
                            positive_state_means = nonnegative_state_means,
                            positive_mean_floor = positive_mean_floor),
              grid = c(grid_selection,
                       list(automatic_prior_sd = automatic_prior_sd,
                            selected_prior_sd = prior_sd)),
              log_likelihood = model_log_likelihood,
              variational_bound = if (inference == "variational") {
                training_fb$variational_bound
              } else NA_real_,
              training_objective = training_objective,
              inference = list(
                method = inference,
                final_exact_used = final_exact_used,
                variational_maxiter = variational_maxiter,
                variational_tolerance = variational_tolerance,
                variational_probability_floor =
                  variational_probability_floor,
                final_variational_iterations =
                  if (inference == "variational") {
                    training_fb$variational_iterations
                  } else NA_integer_,
                final_variational_converged =
                  if (inference == "variational") {
                    training_fb$variational_converged
                  } else NA),
              log_null = log_null,
              log_evidence_ratio = model_log_likelihood - log_null,
              model_selection = model_selection,
              history = history,
              mean_history = do.call(rbind, mean_history),
              pruning_history = pruning_history,
              converged = converged,
              iterations = nrow(history) - 1L)
  class(ans) <- "ash_hmm_fit"
  ans
}


#' Backward-compatible ash-HMM entry point
#'
#' Calls [fit_ash_hmm()] using the original argument names `x` and `sd`.
#'
#' @param x Numeric vector of noisy effect estimates.
#' @param sd Numeric vector of known standard errors.
#' @param ... Additional arguments passed to [fit_ash_hmm()].
#' @return An object of class `ash_hmm_fit`; see [fit_ash_hmm()].
#' @export
fit_hmm <- function(x, sd, ...) {
  fit_ash_hmm(y = x, se = sd, ...)
}

#' Fit the exact two-state binary Markov model
#'
#' Fits the point-state special case with a symmetric transition matrix
#' `A(q) = matrix(c(1-q, q, q, 1-q), 2, 2, byrow = TRUE)`. The scalar `q` is
#' selected by likelihood maximization over the requested interval.
#'
#' @param y Numeric vector of noisy observations.
#' @param se Numeric vector of known, strictly positive standard errors.
#' @param state_means Two finite point-state means.
#' @param sequence_id Vector identifying independent contiguous sequences.
#' @param q_interval Increasing length-two search interval contained in `[0,1]`.
#' @param grid_size Number of initial likelihood evaluations over `q_interval`.
#' @param init_prob Nonnegative length-two initial state distribution.
#' @return A list containing the estimated flip probability `q`, transition
#'   matrix, log likelihood, smoothed state probabilities, and boundary
#'   probabilities.
#' @export
fit_binary_markov <- function(y, se, state_means = c(0, 1),
                              sequence_id = rep(1L, length(y)),
                              q_interval = c(0, 0.5),
                              grid_size = 501L,
                              init_prob = c(0.5, 0.5)) {
  .ash_hmm_check_data(y, se, sequence_id)
  if (!is.numeric(state_means) || length(state_means) != 2L ||
      anyNA(state_means) || any(!is.finite(state_means))) {
    stop("'state_means' must contain two finite values.", call. = FALSE)
  }
  if (length(q_interval) != 2L || q_interval[1L] < 0 ||
      q_interval[2L] > 1 || q_interval[1L] >= q_interval[2L]) {
    stop("'q_interval' must be an increasing subinterval of [0, 1].",
         call. = FALSE)
  }
  if (!is.numeric(init_prob) || length(init_prob) != 2L || anyNA(init_prob) ||
      any(!is.finite(init_prob)) || any(init_prob < 0) || sum(init_prob) <= 0) {
    stop("'init_prob' must be a nonnegative length-two vector with positive mass.",
         call. = FALSE)
  }
  grid_size <- as.integer(grid_size)
  if (grid_size < 5L) stop("'grid_size' must be at least five.", call. = FALSE)
  init_prob <- .ash_hmm_normalize(init_prob)
  mask <- matrix(TRUE, 2L, 2L)
  log_emission <- cbind(
    stats::dnorm(y, state_means[1L], se, log = TRUE),
    stats::dnorm(y, state_means[2L], se, log = TRUE))

  evaluate <- function(q, details = FALSE) {
    A <- matrix(c(1 - q, q, q, 1 - q), 2L, 2L, byrow = TRUE)
    if (!details) {
      return(.ash_hmm_forward_log_likelihood(
        log_emission, A, init_prob, sequence_id, mask))
    }
    fb <- .ash_hmm_forward_backward(
      log_emission, A, init_prob, mask, sequence_id)
    list(A = A, fb = fb)
  }

  grid <- seq(q_interval[1L], q_interval[2L], length.out = grid_size)
  values <- vapply(grid, evaluate, numeric(1L))
  local <- which(values[2L:(grid_size - 1L)] >= values[1L:(grid_size - 2L)] &
                   values[2L:(grid_size - 1L)] >= values[3L:grid_size]) + 1L
  candidates <- data.frame(q = c(grid[1L], grid[grid_size]),
                           log_likelihood = c(values[1L], values[grid_size]))
  if (length(local)) {
    for (i in local) {
      refined <- stats::optimize(function(q) -evaluate(q),
                                 interval = c(grid[i - 1L], grid[i + 1L]))
      candidates <- rbind(candidates,
                          data.frame(q = refined$minimum,
                                     log_likelihood = -refined$objective))
    }
  }
  best <- candidates[which.max(candidates$log_likelihood), ]
  details <- evaluate(best$q, details = TRUE)
  viterbi <- .ash_hmm_viterbi(log_emission, details$A, init_prob, mask, sequence_id)

  list(q = best$q,
       expected_run_length = if (best$q > 0) 1 / best$q else Inf,
       state_probability = details$fb$state_probability,
       boundary_probability = details$fb$boundary_probability,
       viterbi_state = viterbi,
       transition = details$A,
       init_prob = init_prob,
       log_likelihood = details$fb$log_likelihood,
       search = list(grid = grid, values = values, candidates = candidates))
}



.ash_hmm_version <- "3.0.0-variational"

#' Return the implementation version.
#' @export
ash_hmm_version <- function() .ash_hmm_version

.ash_hmm_logsumexp <- function(x) {
  if (!length(x)) return(-Inf)
  z <- max(x)
  if (!is.finite(z)) return(z)
  z + log(sum(exp(x - z)))
}

.ash_hmm_row_logsumexp <- function(x) {
  if (is.null(dim(x))) return(.ash_hmm_logsumexp(x))
  z <- apply(x, 1L, max)
  ans <- rep(-Inf, nrow(x))
  ok <- is.finite(z)
  if (any(ok)) {
    ans[ok] <- z[ok] + log(rowSums(exp(x[ok, , drop = FALSE] - z[ok])))
  }
  ans
}

.ash_hmm_safe_log <- function(x) {
  ans <- rep(-Inf, length(x))
  ans[x > 0] <- log(x[x > 0])
  dim(ans) <- dim(x)
  ans
}

.ash_hmm_normalize <- function(x, fallback = NULL) {
  s <- sum(x)
  if (is.finite(s) && s > 0) return(x / s)
  if (is.null(fallback)) stop("Cannot normalize a zero-mass vector.", call. = FALSE)
  fallback / sum(fallback)
}

.ash_hmm_check_data <- function(y, se, sequence_id) {
  if (!is.numeric(y) || !is.numeric(se) || length(y) != length(se)) {
    stop("'y' and 'se' must be numeric vectors of equal length.", call. = FALSE)
  }
  if (!length(y)) stop("At least one observation is required.", call. = FALSE)
  if (anyNA(y) || anyNA(se) || any(!is.finite(y)) || any(!is.finite(se))) {
    stop("'y' and 'se' must be finite and non-missing.", call. = FALSE)
  }
  if (any(se <= 0)) stop("Every standard error must be strictly positive.", call. = FALSE)
  if (length(sequence_id) != length(y) || anyNA(sequence_id)) {
    stop("'sequence_id' must contain one non-missing value per observation.", call. = FALSE)
  }
  invisible(TRUE)
}

# Stephens-style geometric grid, including the point-mass component at zero.
ash_sigma_grid <- function(y, se, multiplier = sqrt(2), sigma_min = NULL,
                           sigma_max = NULL) {
  .ash_hmm_check_data(y, se, rep(1L, length(y)))
  if (!is.numeric(multiplier) || length(multiplier) != 1L || multiplier <= 1) {
    stop("'multiplier' must be a scalar greater than one.", call. = FALSE)
  }
  if (is.null(sigma_min)) sigma_min <- min(se) / 10
  if (is.null(sigma_max)) {
    signal_var <- max(y^2 - se^2)
    sigma_max <- if (signal_var > 0) 2 * sqrt(signal_var) else 8 * sigma_min
  }
  if (sigma_min <= 0 || sigma_max < sigma_min) {
    stop("Require 0 < sigma_min <= sigma_max.", call. = FALSE)
  }
  positive <- sigma_min * multiplier^(0:ceiling(log(sigma_max / sigma_min) /
                                                  log(multiplier)))
  positive[length(positive)] <- sigma_max
  sort(unique(c(0, positive)))
}

# Symmetric state-mean grid with the zero/hub state first.
ash_mu_grid <- function(y, half_grid = 100L, shape = 3, expansion = 1.5,
                        max_abs = NULL, positive_only = FALSE) {
  if (!is.numeric(y) || !length(y) || anyNA(y) || any(!is.finite(y))) {
    stop("'y' must be a finite numeric vector.", call. = FALSE)
  }
  half_grid <- as.integer(half_grid)
  if (half_grid < 1L || shape <= 0 || expansion <= 0) {
    stop("Require half_grid >= 1, shape > 0, and expansion > 0.", call. = FALSE)
  }
  if (is.null(max_abs)) max_abs <- max(abs(y))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  positive <- seq(0, 1, length.out = half_grid)^(1 / shape) * expansion * max_abs
  if (isTRUE(positive_only)) return(positive)
  c(positive, -positive[-1L])
}

# Screen a candidate state-mean grid before any HMM forward-backward pass.
# Each observation is assigned soft iid weights under point emissions centered
# on the candidate means. Their column sums are effective observation counts.
.ash_hmm_prefilter_mu <- function(y, se, mu,
                                  min_state_count = 0.5,
                                  min_state_fraction = 1e-5,
                                  screening_prior_sd = 0,
                                  block_size = 5000L) {
  if (!is.numeric(mu) || !length(mu) || anyNA(mu) || any(!is.finite(mu))) {
    stop("'mu' must be a finite numeric vector.", call. = FALSE)
  }
  if (length(min_state_count) != 1L || !is.finite(min_state_count) ||
      min_state_count < 0 || length(min_state_fraction) != 1L ||
      !is.finite(min_state_fraction) || min_state_fraction < 0 ||
      min_state_fraction > 1) {
    stop("Grid-screening count and fraction thresholds must be valid scalars.",
         call. = FALSE)
  }
  if (length(screening_prior_sd) != 1L || !is.finite(screening_prior_sd) ||
      screening_prior_sd < 0) {
    stop("'screening_prior_sd' must be a finite nonnegative scalar.",
         call. = FALSE)
  }
  block_size <- as.integer(block_size)
  if (block_size < 1L) stop("'block_size' must be positive.", call. = FALSE)

  n <- length(y)
  m <- length(mu)
  effective_count <- numeric(m)
  starts <- seq.int(1L, n, by = block_size)

  for (first in starts) {
    last <- min(n, first + block_size - 1L)
    index <- first:last
    log_score <- matrix(-Inf, length(index), m)
    for (state in seq_len(m)) {
      log_score[, state] <- stats::dnorm(
        y[index], mean = mu[state],
        sd = sqrt(se[index]^2 + screening_prior_sd^2), log = TRUE)
    }
    log_normalizer <- .ash_hmm_row_logsumexp(log_score)
    effective_count <- effective_count +
      colSums(exp(log_score - log_normalizer))
  }

  cutoff <- max(min_state_count, n * min_state_fraction)
  keep <- which(effective_count >= cutoff)
  # State 1 is the zero/hub state for automatically constructed grids and is
  # retained even when the data are entirely nonzero.
  keep <- sort(unique(c(1L, keep)))

  list(full_mu = mu,
       selected_mu = mu[keep],
       selected_index = keep,
       removed_index = setdiff(seq_len(m), keep),
       effective_count = effective_count,
       mean_weight = effective_count / n,
       cutoff = cutoff,
       settings = list(min_state_count = min_state_count,
                       min_state_fraction = min_state_fraction,
                       screening_prior_sd = screening_prior_sd,
                       block_size = block_size))
}

# Construct the original symmetric power grid and remove unsupported states
# before fitting the HMM. This is also useful for inspecting the automatic
# choice without fitting a model.
ash_hmm_select_grid <- function(y, se,
                                half_grid = 100L,
                                shape = 3,
                                expansion = 1.5,
                                max_abs = NULL,
                                min_state_count = 0.5,
                                min_state_fraction = 1e-5,
                                screening_prior_sd = 0,
                                block_size = 5000L,
                                positive_only = FALSE) {
  .ash_hmm_check_data(y, se, rep(1L, length(y)))
  candidate <- ash_mu_grid(y, half_grid = half_grid, shape = shape,
                           expansion = expansion, max_abs = max_abs,
                           positive_only = positive_only)
  ans <- .ash_hmm_prefilter_mu(
    y, se, candidate,
    min_state_count = min_state_count,
    min_state_fraction = min_state_fraction,
    screening_prior_sd = screening_prior_sd,
    block_size = block_size)
  ans$automatic <- TRUE
  ans
}

ash_hmm_transition_mask <- function(n_states, topology = c("hub", "full")) {
  topology <- match.arg(topology)
  n_states <- as.integer(n_states)
  if (n_states < 1L) stop("'n_states' must be positive.", call. = FALSE)
  if (topology == "full") return(matrix(TRUE, n_states, n_states))

  # Hub-and-spoke: the zero/hub state is state 1. A non-hub state may
  # persist or return to the hub, but may not jump directly to another spoke.
  mask <- diag(TRUE, n_states)
  mask[1L, ] <- TRUE
  mask[, 1L] <- TRUE
  mask
}

.ash_hmm_initial_transition <- function(mask, stay_probability = 0.95) {
  m <- nrow(mask)
  stay_probability <- rep(stay_probability, length.out = m)
  if (any(stay_probability < 0 | stay_probability > 1)) {
    stop("'stay_probability' must lie in [0, 1].", call. = FALSE)
  }
  A <- matrix(0, m, m)
  for (i in seq_len(m)) {
    allowed <- which(mask[i, ])
    if (length(allowed) == 1L) {
      A[i, allowed] <- 1
    } else {
      A[i, i] <- stay_probability[i]
      other <- setdiff(allowed, i)
      A[i, other] <- (1 - stay_probability[i]) / length(other)
    }
  }
  A
}

.ash_hmm_stationary <- function(A, tolerance = 1e-13, maxiter = 100000L) {
  p <- rep(1 / nrow(A), nrow(A))
  for (iter in seq_len(maxiter)) {
    p_new <- drop(p %*% A)
    if (max(abs(p_new - p)) < tolerance) return(.ash_hmm_normalize(p_new))
    p <- p_new
  }
  warning("Stationary-distribution iteration did not converge; using its last iterate.")
  .ash_hmm_normalize(p)
}

.ash_hmm_validate_transition <- function(A, mask, tolerance = 1e-10) {
  m <- nrow(mask)
  if (!is.matrix(A) || any(dim(A) != c(m, m)) || anyNA(A) ||
      any(!is.finite(A)) || any(A < 0)) {
    stop("'init_transition' must be a finite, nonnegative M by M matrix.", call. = FALSE)
  }
  if (any(A[!mask] > tolerance)) {
    stop("'init_transition' has positive mass on a structurally forbidden edge.", call. = FALSE)
  }
  rs <- rowSums(A)
  if (any(abs(rs - 1) > tolerance)) {
    stop("Every row of 'init_transition' must sum to one.", call. = FALSE)
  }
  A[!mask] <- 0
  A / rowSums(A)
}

# Log marginal density for one ash component. In the one-sided model the
# continuous prior is N(mu, tau^2) truncated to [0, Inf). The identity
#
# p(y | beta >= 0) = p_untruncated(y) *
#   Pr(beta >= 0 | y) / Pr(beta >= 0)
#
# gives an exact and stable normal-CDF expression for the convolution with
# N(beta, se^2) observation noise.
.ash_hmm_component_log_emission <- function(
    y, se, mu, tau, effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  if (tau == 0) {
    if (effect_support == "nonnegative" && mu < 0) {
      return(rep(-Inf, length(y)))
    }
    return(stats::dnorm(y, mean = mu, sd = se, log = TRUE))
  }

  observation_variance <- se^2 + tau^2
  ans <- stats::dnorm(
    y, mean = mu, sd = sqrt(observation_variance), log = TRUE)
  if (effect_support == "real") return(ans)

  posterior_mean <- (se^2 * mu + tau^2 * y) / observation_variance
  posterior_sd <- tau * se / sqrt(observation_variance)
  ans + stats::pnorm(
    posterior_mean / posterior_sd, log.p = TRUE) -
    stats::pnorm(mu / tau, log.p = TRUE)
}

.ash_hmm_component_log_emission_matrix <- function(
    y, se, mu, prior_sd, effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  ans <- matrix(-Inf, length(y), length(prior_sd))
  for (component in seq_along(prior_sd)) {
    ans[, component] <- .ash_hmm_component_log_emission(
      y, se, mu, prior_sd[component], effect_support)
  }
  ans
}

.ash_hmm_log_emission <- function(
    y, se, mu, prior_sd, rho,
    effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  ans <- matrix(-Inf, n, m)
  log_rho <- .ash_hmm_safe_log(rho)
  for (state in seq_len(m)) {
    terms <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, effect_support)
    terms <- sweep(terms, 2L, log_rho[state, ], "+")
    ans[, state] <- .ash_hmm_row_logsumexp(terms)
  }
  ans
}

# Cache every component log emission. The state-specific mixture weights change
# at every M-step, but these component densities change only when a state center
# moves. Keeping this array avoids evaluating the same normal densities once
# for the mixture likelihood and again for the component responsibilities.
.ash_hmm_component_log_emission_array <- function(
    y, se, mu, prior_sd,
    effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  n <- length(y)
  l <- length(prior_sd)
  m <- length(mu)
  answer <- array(-Inf, c(n, l, m))
  for (state in seq_len(m)) {
    answer[, , state] <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, effect_support)
  }
  answer
}

.ash_hmm_mix_component_cache <- function(component_log_emission, rho) {
  dimensions <- dim(component_log_emission)
  n <- dimensions[1L]
  l <- dimensions[2L]
  m <- dimensions[3L]
  if (!is.matrix(rho) || any(dim(rho) != c(m, l))) {
    stop("Cached component emissions and mixture weights are incompatible.",
         call. = FALSE)
  }
  answer <- matrix(-Inf, n, m)
  log_rho <- .ash_hmm_safe_log(rho)
  for (state in seq_len(m)) {
    component <- matrix(
      component_log_emission[, , state], n, l)
    answer[, state] <- .ash_hmm_row_logsumexp(
      sweep(component, 2L, log_rho[state, ], "+"))
  }
  answer
}

.ash_hmm_refresh_component_cache <- function(
    component_log_emission, y, se, old_mu, new_mu, prior_sd,
    effect_support = c("real", "nonnegative"),
    tolerance = 0) {
  effect_support <- match.arg(effect_support)
  changed <- which(abs(new_mu - old_mu) > tolerance)
  if (!length(changed)) return(component_log_emission)
  n <- length(y)
  l <- length(prior_sd)
  for (state in changed) {
    component_log_emission[, , state] <-
      .ash_hmm_component_log_emission_matrix(
        y, se, new_mu[state], prior_sd, effect_support)
  }
  dim(component_log_emission) <- c(n, l, length(new_mu))
  component_log_emission
}

.ash_hmm_component_statistics_cached <- function(
    y, se, mu, prior_sd, rho, log_emission, gamma,
    component_log_emission,
    effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  if (!identical(dim(component_log_emission), c(n, l, m))) {
    stop("The component-emission cache has incompatible dimensions.",
         call. = FALSE)
  }
  counts <- matrix(0, m, l)
  mean_numerator <- mean_precision <- numeric(m)
  log_rho <- .ash_hmm_safe_log(rho)
  if (effect_support == "real") {
    precision <- 1 / outer(se^2, prior_sd^2, "+")
    weighted_y <- y * precision
  }
  for (state in seq_len(m)) {
    component <- matrix(
      component_log_emission[, , state], n, l)
    conditional <- exp(
      sweep(component, 2L, log_rho[state, ], "+") -
        log_emission[, state])
    joint <- conditional * gamma[, state]
    counts[state, ] <- colSums(joint)
    if (effect_support == "real") {
      mean_numerator[state] <- sum(joint * weighted_y)
      mean_precision[state] <- sum(joint * precision)
    }
  }
  list(
    counts = counts,
    occupancy = rowSums(counts),
    mean_numerator = mean_numerator,
    mean_precision = mean_precision)
}

.ash_hmm_segments <- function(sequence_id) {
  starts <- c(1L, which(sequence_id[-1L] != sequence_id[-length(sequence_id)]) + 1L)
  ends <- c(starts[-1L] - 1L, length(sequence_id))
  list(starts = starts, ends = ends)
}

.ash_hmm_forward_backward_log <- function(log_emission, A, init_prob, mask,
                                          sequence_id) {
  n <- nrow(log_emission)
  m <- ncol(log_emission)
  log_A <- .ash_hmm_safe_log(A)
  log_pi <- .ash_hmm_safe_log(init_prob)
  edge <- which(mask, arr.ind = TRUE)
  edge_from <- edge[, 1L]
  edge_to <- edge[, 2L]
  incoming <- lapply(seq_len(m), function(j) which(edge_to == j))
  outgoing <- lapply(seq_len(m), function(i) which(edge_from == i))
  segments <- .ash_hmm_segments(sequence_id)

  log_alpha <- matrix(-Inf, n, m)
  log_beta <- matrix(-Inf, n, m)
  gamma <- matrix(0, n, m)
  transition_counts <- matrix(0, m, m)
  boundary_probability <- rep(NA_real_, max(0L, n - 1L))
  total_log_likelihood <- 0

  for (segment in seq_along(segments$starts)) {
    first <- segments$starts[segment]
    last <- segments$ends[segment]
    log_alpha[first, ] <- log_pi + log_emission[first, ]

    if (first < last) {
      for (t in (first + 1L):last) {
        for (state in seq_len(m)) {
          e <- incoming[[state]]
          log_alpha[t, state] <- log_emission[t, state] +
            .ash_hmm_logsumexp(log_alpha[t - 1L, edge_from[e]] +
                                 log_A[cbind(edge_from[e], edge_to[e])])
        }
      }
    }

    segment_ll <- .ash_hmm_logsumexp(log_alpha[last, ])
    if (!is.finite(segment_ll)) {
      stop("The HMM assigns zero probability to an observed sequence.", call. = FALSE)
    }
    total_log_likelihood <- total_log_likelihood + segment_ll
    log_beta[last, ] <- 0

    if (first < last) {
      for (t in (last - 1L):first) {
        for (state in seq_len(m)) {
          e <- outgoing[[state]]
          log_beta[t, state] <- .ash_hmm_logsumexp(
            log_A[cbind(edge_from[e], edge_to[e])] +
              log_emission[t + 1L, edge_to[e]] +
              log_beta[t + 1L, edge_to[e]])
        }
      }
    }

    for (t in first:last) {
      lg <- log_alpha[t, ] + log_beta[t, ] - segment_ll
      lg <- lg - .ash_hmm_logsumexp(lg)
      gamma[t, ] <- exp(lg)
    }

    if (first < last) {
      for (t in first:(last - 1L)) {
        log_xi <- log_alpha[t, edge_from] +
          log_A[cbind(edge_from, edge_to)] +
          log_emission[t + 1L, edge_to] +
          log_beta[t + 1L, edge_to] - segment_ll
        log_xi <- log_xi - .ash_hmm_logsumexp(log_xi)
        xi <- exp(log_xi)
        for (e in seq_along(xi)) {
          transition_counts[edge_from[e], edge_to[e]] <-
            transition_counts[edge_from[e], edge_to[e]] + xi[e]
        }
        boundary_probability[t] <- sum(xi[edge_from != edge_to])
      }
    }
  }

  list(log_likelihood = total_log_likelihood,
       state_probability = gamma,
       transition_counts = transition_counts,
       boundary_probability = boundary_probability,
       initial_counts = colSums(gamma[segments$starts, , drop = FALSE]),
       starts = segments$starts)
}

# Fast scaled-probability forward-backward recursion. Subtracting the largest
# log emission separately at every position prevents overflow, and the usual
# forward scale factors prevent underflow across long sequences. Matrix
# multiplication moves the state recursion out of interpreted R loops.
#
# The function returns NULL in the exceptionally ill-conditioned case where
# all reachable scaled emissions underflow at a position. The dispatcher then
# falls back to the fully log-domain implementation above.
.ash_hmm_forward_backward_scaled <- function(log_emission, A, init_prob, mask,
                                             sequence_id) {
  n <- nrow(log_emission)
  m <- ncol(log_emission)
  segments <- .ash_hmm_segments(sequence_id)
  row_shift <- apply(log_emission, 1L, max)
  if (any(!is.finite(row_shift))) return(NULL)
  emission <- exp(log_emission - row_shift)

  A_work <- A
  A_work[!mask] <- 0
  alpha <- matrix(0, n, m)
  beta <- matrix(0, n, m)
  gamma <- matrix(0, n, m)
  scale <- numeric(n)
  transition_counts <- matrix(0, m, m)
  boundary_probability <- rep(NA_real_, max(0L, n - 1L))
  total_log_likelihood <- 0

  for (segment in seq_along(segments$starts)) {
    first <- segments$starts[segment]
    last <- segments$ends[segment]

    forward <- init_prob * emission[first, ]
    scale[first] <- sum(forward)
    if (!is.finite(scale[first]) || scale[first] <= 0) return(NULL)
    alpha[first, ] <- forward / scale[first]
    total_log_likelihood <- total_log_likelihood +
      row_shift[first] + log(scale[first])

    if (first < last) {
      for (t in (first + 1L):last) {
        forward <- drop(alpha[t - 1L, ] %*% A_work) * emission[t, ]
        scale[t] <- sum(forward)
        if (!is.finite(scale[t]) || scale[t] <= 0) return(NULL)
        alpha[t, ] <- forward / scale[t]
        total_log_likelihood <- total_log_likelihood +
          row_shift[t] + log(scale[t])
      }
    }

    beta[last, ] <- 1
    if (first < last) {
      for (t in (last - 1L):first) {
        next_weight <- emission[t + 1L, ] * beta[t + 1L, ]
        beta[t, ] <- drop(A_work %*% next_weight) / scale[t + 1L]
      }
    }

    index <- first:last
    smoothed <- alpha[index, , drop = FALSE] *
      beta[index, , drop = FALSE]
    smoothed_sum <- rowSums(smoothed)
    if (any(!is.finite(smoothed_sum)) || any(smoothed_sum <= 0)) return(NULL)
    gamma[index, ] <- smoothed / smoothed_sum

    if (first < last) {
      from <- first:(last - 1L)
      to <- from + 1L
      left <- alpha[from, , drop = FALSE]
      right <- emission[to, , drop = FALSE] *
        beta[to, , drop = FALSE]
      pair_mass <- rowSums((left %*% A_work) * right)
      if (any(!is.finite(pair_mass)) || any(pair_mass <= 0)) return(NULL)
      right <- right / pair_mass

      # sum_t xi_t = A * sum_t alpha_t right_t^T. This replaces T
      # allocations of an M by M matrix by one BLAS cross-product.
      transition_counts <- transition_counts +
        A_work * crossprod(left, right)
      same_state <- rowSums(
        left * sweep(right, 2L, diag(A_work), "*"))
      boundary_probability[from] <- pmin(1, pmax(0, 1 - same_state))
    }
  }

  list(log_likelihood = total_log_likelihood,
       state_probability = gamma,
       transition_counts = transition_counts,
       boundary_probability = boundary_probability,
       initial_counts = colSums(gamma[segments$starts, , drop = FALSE]),
       starts = segments$starts)
}

# The dense scaled engine is substantially faster for the full topology. The
# sparse log engine remains useful for large hub-and-spoke masks and is also a
# numerical fallback. Advanced users can override the automatic choice with
# options(ashHMM.forward_backward = "scaled") or "log".
.ash_hmm_forward_backward <- function(
    log_emission, A, init_prob, mask, sequence_id,
    engine = getOption("ashHMM.forward_backward", "auto")) {
  engine <- match.arg(engine, c("auto", "scaled", "log"))
  if (engine == "auto") {
    density <- sum(mask) / length(mask)
    engine <- if (nrow(mask) <= 4L || density >= 0.25) "scaled" else "log"
  }
  if (engine == "log") {
    return(.ash_hmm_forward_backward_log(
      log_emission, A, init_prob, mask, sequence_id))
  }
  ans <- .ash_hmm_forward_backward_scaled(
    log_emission, A, init_prob, mask, sequence_id)
  if (!is.null(ans)) return(ans)
  .ash_hmm_forward_backward_log(
    log_emission, A, init_prob, mask, sequence_id)
}

.ash_hmm_softmax_rows <- function(score, probability_floor = 0) {
  shift <- apply(score, 1L, max)
  if (any(!is.finite(shift))) {
    stop("A variational state update has no finite support.", call. = FALSE)
  }
  probability <- exp(score - shift)
  probability <- probability / rowSums(probability)
  if (probability_floor > 0) {
    probability <- pmax(probability, probability_floor)
    dim(probability) <- dim(score)
    probability <- probability / rowSums(probability)
  }
  probability
}

.ash_hmm_variational_bound <- function(
    log_emission, A, init_prob, state_probability, sequence_id) {
  n <- nrow(log_emission)
  segments <- .ash_hmm_segments(sequence_id)
  starts <- segments$starts
  within <- if (n > 1L) {
    which(sequence_id[-n] == sequence_id[-1L])
  } else integer()
  log_A <- log(A)
  expected_emission <- sum(state_probability * log_emission)
  expected_initial <- sum(
    state_probability[starts, , drop = FALSE] *
      matrix(log(init_prob), length(starts), ncol(log_emission),
             byrow = TRUE))
  expected_transition <- if (length(within)) {
    sum((state_probability[within, , drop = FALSE] %*% log_A) *
          state_probability[within + 1L, , drop = FALSE])
  } else 0
  use <- state_probability > 0
  entropy <- -sum(state_probability[use] * log(state_probability[use]))
  expected_emission + expected_initial + expected_transition + entropy
}

# Fully factorized variational HMM update q(Q)=prod_t q_t(Q_t).
#
# A chain graph is bipartite. Conditional on all even positions, the odd
# positions form an independent coordinate block, and conversely. Alternating
# these two vectorized checkerboard blocks is therefore genuine coordinate
# ascent on the mean-field ELBO, rather than a synchronous fixed-point
# heuristic. Hard structural zeros are intentionally disallowed because a
# product distribution would otherwise assign positive mass to forbidden
# adjacent pairs and have ELBO -Inf.
.ash_hmm_variational_states <- function(
    log_emission, A, init_prob, mask, sequence_id,
    initial_probability = NULL,
    maxiter = 50L, tolerance = 1e-6,
    probability_floor = 1e-12) {
  n <- nrow(log_emission)
  m <- ncol(log_emission)
  if (!all(mask)) {
    stop(paste(
      "Variational state inference requires a full transition mask;",
      "use inference='exact' for hub or custom structural-zero topologies."),
      call. = FALSE)
  }
  if (any(A <= 0) || any(init_prob <= 0)) {
    stop(paste(
      "Variational state inference requires positive transition and initial",
      "probabilities."), call. = FALSE)
  }
  maxiter <- as.integer(maxiter)
  if (length(maxiter) != 1L || is.na(maxiter) || maxiter < 1L ||
      !is.finite(tolerance) || tolerance <= 0 ||
      !is.finite(probability_floor) || probability_floor < 0 ||
      probability_floor >= 1 / m) {
    stop("Invalid variational state-inference controls.", call. = FALSE)
  }
  if (is.null(initial_probability)) {
    state_probability <- .ash_hmm_softmax_rows(
      log_emission, probability_floor)
  } else {
    if (!is.matrix(initial_probability) ||
        any(dim(initial_probability) != c(n, m)) ||
        anyNA(initial_probability) ||
        any(!is.finite(initial_probability)) ||
        any(initial_probability < 0) ||
        any(rowSums(initial_probability) <= 0)) {
      stop("'initial_probability' must be a nonnegative T by M matrix.",
           call. = FALSE)
    }
    state_probability <- initial_probability /
      rowSums(initial_probability)
    state_probability <- pmax(state_probability, probability_floor)
    dim(state_probability) <- c(n, m)
    state_probability <- state_probability /
      rowSums(state_probability)
  }

  segments <- .ash_hmm_segments(sequence_id)
  is_start <- rep(FALSE, n)
  is_end <- rep(FALSE, n)
  is_start[segments$starts] <- TRUE
  is_end[segments$ends] <- TRUE
  log_A <- log(A)
  log_init <- log(init_prob)
  converged <- FALSE
  previous_bound <- -Inf
  maximum_change <- Inf

  for (iteration in seq_len(maxiter)) {
    previous <- state_probability
    for (parity in 0:1) {
      index <- which(seq_len(n) %% 2L == parity)
      if (!length(index)) next
      score <- log_emission[index, , drop = FALSE]

      at_start <- is_start[index]
      if (any(at_start)) {
        score[at_start, ] <- score[at_start, , drop = FALSE] +
          matrix(log_init, sum(at_start), m, byrow = TRUE)
      }
      has_left <- !at_start
      if (any(has_left)) {
        score[has_left, ] <- score[has_left, , drop = FALSE] +
          state_probability[index[has_left] - 1L, , drop = FALSE] %*%
          log_A
      }
      has_right <- !is_end[index]
      if (any(has_right)) {
        score[has_right, ] <- score[has_right, , drop = FALSE] +
          state_probability[index[has_right] + 1L, , drop = FALSE] %*%
          t(log_A)
      }
      state_probability[index, ] <- .ash_hmm_softmax_rows(
        score, probability_floor)
    }

    maximum_change <- max(abs(state_probability - previous))
    bound <- .ash_hmm_variational_bound(
      log_emission, A, init_prob, state_probability, sequence_id)
    if (iteration > 1L &&
        maximum_change <= sqrt(tolerance) &&
        abs(bound - previous_bound) <=
        tolerance * (1 + abs(previous_bound))) {
      converged <- TRUE
      break
    }
    previous_bound <- bound
  }

  within <- if (n > 1L) {
    which(sequence_id[-n] == sequence_id[-1L])
  } else integer()
  transition_counts <- if (length(within)) {
    crossprod(
      state_probability[within, , drop = FALSE],
      state_probability[within + 1L, , drop = FALSE])
  } else matrix(0, m, m)
  boundary_probability <- rep(NA_real_, max(0L, n - 1L))
  if (length(within)) {
    boundary_probability[within] <- 1 -
      rowSums(
        state_probability[within, , drop = FALSE] *
          state_probability[within + 1L, , drop = FALSE])
  }
  list(
    log_likelihood = bound,
    variational_bound = bound,
    state_probability = state_probability,
    transition_counts = transition_counts,
    boundary_probability = boundary_probability,
    initial_counts = colSums(
      state_probability[segments$starts, , drop = FALSE]),
    starts = segments$starts,
    variational_iterations = iteration,
    variational_converged = converged,
    variational_maximum_change = maximum_change)
}

# Likelihood-only scaled forward recursion. Optimization routines should use
# this path when smoothing probabilities and transition counts are unnecessary.
.ash_hmm_forward_log_likelihood <- function(
    log_emission, A, init_prob, sequence_id, mask = A > 0) {
  n <- nrow(log_emission)
  row_shift <- apply(log_emission, 1L, max)
  if (any(!is.finite(row_shift))) {
    return(.ash_hmm_forward_backward_log(
      log_emission, A, init_prob, mask, sequence_id)$log_likelihood)
  }
  emission <- exp(log_emission - row_shift)
  A_work <- A
  A_work[!mask] <- 0
  segments <- .ash_hmm_segments(sequence_id)
  total_log_likelihood <- 0

  for (segment in seq_along(segments$starts)) {
    first <- segments$starts[segment]
    last <- segments$ends[segment]
    forward <- init_prob * emission[first, ]
    scale <- sum(forward)
    if (!is.finite(scale) || scale <= 0) {
      return(.ash_hmm_forward_backward_log(
        log_emission, A, init_prob, mask, sequence_id)$log_likelihood)
    }
    forward <- forward / scale
    total_log_likelihood <- total_log_likelihood +
      row_shift[first] + log(scale)

    if (first < last) {
      for (t in (first + 1L):last) {
        forward <- drop(forward %*% A_work) * emission[t, ]
        scale <- sum(forward)
        if (!is.finite(scale) || scale <= 0) {
          return(.ash_hmm_forward_backward_log(
            log_emission, A, init_prob, mask, sequence_id)$log_likelihood)
        }
        forward <- forward / scale
        total_log_likelihood <- total_log_likelihood +
          row_shift[t] + log(scale)
      }
    }
  }
  total_log_likelihood
}

.ash_hmm_component_statistics <- function(
    y, se, mu, prior_sd, rho, log_emission, gamma,
    effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  counts <- matrix(0, m, l)
  mean_numerator <- mean_precision <- numeric(m)
  log_rho <- .ash_hmm_safe_log(rho)
  for (state in seq_len(m)) {
    terms <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, effect_support)
    terms <- sweep(terms, 2L, log_rho[state, ], "+")
    conditional <- exp(terms - log_emission[, state])
    joint <- conditional * gamma[, state]
    counts[state, ] <- colSums(joint)
    if (effect_support == "real") {
      for (component in seq_len(l)) {
        observation_variance <- se^2 + prior_sd[component]^2
        mean_numerator[state] <- mean_numerator[state] +
          sum(joint[, component] * y / observation_variance)
        mean_precision[state] <- mean_precision[state] +
          sum(joint[, component] / observation_variance)
      }
    }
  }
  list(counts = counts,
       occupancy = rowSums(counts),
       mean_numerator = mean_numerator,
       mean_precision = mean_precision)
}

.ash_hmm_component_counts <- function(
    y, se, mu, prior_sd, rho, log_emission, gamma,
    effect_support = c("real", "nonnegative")) {
  .ash_hmm_component_statistics(y, se, mu, prior_sd, rho,
                                log_emission, gamma,
                                effect_support)$counts
}

# Nonoverlapping intervals around the current set of anchor centers. Restricting
# each learned center to its own interval prevents label crossing and stops a
# dense collection of states from collapsing onto the same plateau. When states
# are pruned, recomputing these intervals automatically gives the survivors more
# room to move.
.ash_hmm_voronoi_bounds <- function(anchor) {
  m <- length(anchor)
  lower <- rep(-Inf, m)
  upper <- rep(Inf, m)
  if (m <= 1L) return(cbind(lower = lower, upper = upper))
  if (anyDuplicated(anchor)) {
    stop("State-center anchors must be distinct when using Voronoi bounds.",
         call. = FALSE)
  }
  order_mu <- order(anchor)
  sorted <- anchor[order_mu]
  midpoint <- (sorted[-m] + sorted[-1L]) / 2
  lower_sorted <- c(-Inf, midpoint)
  upper_sorted <- c(midpoint, Inf)
  lower[order_mu] <- lower_sorted
  upper[order_mu] <- upper_sorted
  cbind(lower = lower, upper = upper)
}

.ash_hmm_update_means <- function(mu, anchor, statistics,
                                  learn, fixed, minimum_count = 2,
                                  eligible = rep(TRUE, length(mu)),
                                  damping = 1,
                                  minimum = rep(-Inf, length(mu)),
                                  bounds = c("voronoi", "none")) {
  bounds <- match.arg(bounds)
  m <- length(mu)
  if (length(eligible) != m || anyNA(eligible)) {
    stop("'eligible' must contain one non-missing value per state.",
         call. = FALSE)
  }
  eligible <- as.logical(eligible)
  damping <- rep(damping, length.out = m)
  if (anyNA(damping) || any(!is.finite(damping)) ||
      any(damping < 0 | damping > 1)) {
    stop("Every mean damping value must lie in [0, 1].", call. = FALSE)
  }
  minimum <- rep(minimum, length.out = m)
  if (anyNA(minimum)) {
    stop("Mean lower bounds must not be missing.", call. = FALSE)
  }
  if (!learn || !m) {
    return(list(mu = mu, raw = mu, updated = rep(FALSE, m),
                lower = rep(-Inf, m), upper = rep(Inf, m)))
  }
  interval <- if (bounds == "voronoi") {
    .ash_hmm_voronoi_bounds(anchor)
  } else {
    cbind(lower = rep(-Inf, m), upper = rep(Inf, m))
  }
  interval[, "lower"] <- pmax(interval[, "lower"], minimum)
  raw <- mu
  estimable <- eligible & statistics$occupancy >= minimum_count &
    is.finite(statistics$mean_precision) & statistics$mean_precision > 0
  estimable[fixed] <- FALSE
  raw[estimable] <- statistics$mean_numerator[estimable] /
    statistics$mean_precision[estimable]
  target <- pmax(interval[, "lower"], pmin(interval[, "upper"], raw))
  updated <- estimable & is.finite(target) & damping > 0
  ans <- mu
  ans[updated] <- mu[updated] +
    damping[updated] * (target[updated] - mu[updated])
  ans[fixed] <- mu[fixed]
  list(mu = ans, raw = raw, updated = updated,
       lower = interval[, "lower"], upper = interval[, "upper"])
}

# Conditional maximization for centers in the one-sided ash model. Component
# responsibilities are held at the current parameter snapshot, so each state's
# criterion is a one-dimensional weighted sum of exact truncated-normal
# component log emissions. Candidate moves are accepted only when this
# conditional criterion does not decrease; damping is backtracked if needed.
.ash_hmm_update_truncated_means <- function(
    y, se, mu, anchor, prior_sd, rho, log_emission, gamma, statistics,
    learn, fixed, minimum_count = 2,
    eligible = rep(TRUE, length(mu)), damping = 1,
    minimum = rep(0, length(mu)), bounds = c("voronoi", "none")) {
  bounds <- match.arg(bounds)
  m <- length(mu)
  if (length(eligible) != m || anyNA(eligible)) {
    stop("'eligible' must contain one non-missing value per state.",
         call. = FALSE)
  }
  eligible <- as.logical(eligible)
  damping <- rep(damping, length.out = m)
  if (anyNA(damping) || any(!is.finite(damping)) ||
      any(damping < 0 | damping > 1)) {
    stop("Every mean damping value must lie in [0, 1].", call. = FALSE)
  }
  minimum <- rep(minimum, length.out = m)
  if (anyNA(minimum) || any(!is.finite(minimum))) {
    stop("One-sided mean lower bounds must be finite.", call. = FALSE)
  }

  interval <- if (bounds == "voronoi") {
    .ash_hmm_voronoi_bounds(anchor)
  } else {
    cbind(lower = rep(-Inf, m), upper = rep(Inf, m))
  }
  interval[, "lower"] <- pmax(interval[, "lower"], minimum)

  # stats::optimize requires finite intervals. This upper cap is deliberately
  # conservative relative to both the starting grid and plausible Gaussian
  # observations; it only affects the outermost anchor cell.
  scale_pad <- max(1e-8, stats::median(se))
  search_upper <- max(c(mu, anchor, y + 8 * se), na.rm = TRUE)
  search_upper <- max(search_upper, max(interval[, "lower"]) + scale_pad)
  infinite_upper <- !is.finite(interval[, "upper"])
  interval[infinite_upper, "upper"] <- pmax(
    search_upper, interval[infinite_upper, "lower"] + scale_pad)

  raw <- ans <- mu
  updated <- rep(FALSE, m)
  criterion_before <- criterion_after <- rep(NA_real_, m)
  if (!learn || !m) {
    return(list(mu = ans, raw = raw, updated = updated,
                lower = interval[, "lower"], upper = interval[, "upper"],
                criterion_before = criterion_before,
                criterion_after = criterion_after))
  }

  estimable <- eligible & statistics$occupancy >= minimum_count & damping > 0
  estimable[fixed] <- FALSE
  log_rho <- .ash_hmm_safe_log(rho)

  for (state in which(estimable)) {
    component_log_current <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, "nonnegative")
    terms <- sweep(
      component_log_current, 2L, log_rho[state, ], "+")
    component_probability <- exp(terms - log_emission[, state])
    joint <- component_probability * gamma[, state]
    positive_weight <- joint > 0
    if (!any(positive_weight)) next

    conditional_objective <- function(center) {
      component_log <- .ash_hmm_component_log_emission_matrix(
        y, se, center, prior_sd, "nonnegative")
      sum(joint[positive_weight] * component_log[positive_weight])
    }

    lower <- interval[state, "lower"]
    upper <- interval[state, "upper"]
    if (!is.finite(lower) || !is.finite(upper) || upper <= lower) next
    current <- min(upper, max(lower, mu[state]))
    current_value <- conditional_objective(current)
    optimum <- stats::optimize(
      function(center) -conditional_objective(center),
      interval = c(lower, upper),
      tol = sqrt(.Machine$double.eps))
    candidates <- unique(c(current, lower, optimum$minimum, upper))
    values <- vapply(candidates, conditional_objective, numeric(1L))
    target <- candidates[which.max(values)]
    target_value <- max(values)
    raw[state] <- target
    criterion_before[state] <- current_value

    if (!is.finite(target_value) ||
        target_value <= current_value + 1e-10 * (1 + abs(current_value))) {
      criterion_after[state] <- current_value
      next
    }

    fraction <- damping[state]
    proposed <- current + fraction * (target - current)
    proposed_value <- conditional_objective(proposed)
    while (fraction > 2^-20 &&
           proposed_value < current_value - 1e-10 * (1 + abs(current_value))) {
      fraction <- fraction / 2
      proposed <- current + fraction * (target - current)
      proposed_value <- conditional_objective(proposed)
    }
    if (is.finite(proposed_value) && proposed_value >= current_value -
        1e-10 * (1 + abs(current_value))) {
      ans[state] <- proposed
      updated[state] <- abs(proposed - mu[state]) > 1e-12
      criterion_after[state] <- proposed_value
    } else {
      criterion_after[state] <- current_value
    }
  }

  ans[fixed] <- mu[fixed]
  list(mu = ans, raw = raw, updated = updated,
       lower = interval[, "lower"], upper = interval[, "upper"],
       criterion_before = criterion_before,
       criterion_after = criterion_after)
}

# Effective dimension used by the optional BIC/extended-BIC comparison with
# the strict all-null model. Only numerically active simplex coordinates count;
# this is a pragmatic diagnostic because ordinary BIC regularity fails on
# mixture boundaries.
.ash_hmm_effective_parameter_count <- function(
    A, mask, rho, fixed_pointmass_states = integer(),
    learned_mean_states = integer(), estimate_init = FALSE,
    init_prob = NULL, tolerance = 1e-6) {
  active_transition <- mask & A > tolerance
  transition_df <- sum(pmax(rowSums(active_transition) - 1L, 0L))
  mixture_df <- 0L
  free <- setdiff(seq_len(nrow(rho)), fixed_pointmass_states)
  if (length(free)) {
    mixture_df <- sum(pmax(rowSums(rho[free, , drop = FALSE] > tolerance) - 1L, 0L))
  }
  mean_df <- length(unique(learned_mean_states))
  init_df <- 0L
  if (estimate_init && !is.null(init_prob)) {
    init_df <- max(sum(init_prob > tolerance) - 1L, 0L)
  }
  as.integer(transition_df + mixture_df + mean_df + init_df)
}

# Decode a state path with an explicit l0 penalty for every state change. This
# separates step-count selection from the fitted Markov transition matrix.
.ash_hmm_penalized_viterbi <- function(log_emission, mask, sequence_id,
                                       change_penalty) {
  n <- nrow(log_emission)
  m <- ncol(log_emission)
  if (!is.numeric(change_penalty) || length(change_penalty) != 1L ||
      !is.finite(change_penalty) || change_penalty < 0) {
    stop("'change_penalty' must be a finite nonnegative scalar.", call. = FALSE)
  }
  segments <- .ash_hmm_segments(sequence_id)
  path <- integer(n)
  for (segment in seq_along(segments$starts)) {
    first <- segments$starts[segment]
    last <- segments$ends[segment]
    len <- last - first + 1L
    delta <- matrix(-Inf, len, m)
    back <- matrix(NA_integer_, len, m)
    delta[1L, ] <- log_emission[first, ]
    if (len > 1L) {
      for (tt in 2L:len) {
        for (state in seq_len(m)) {
          from <- which(mask[, state])
          value <- delta[tt - 1L, from] -
            change_penalty * as.numeric(from != state)
          winner <- which.max(value)
          delta[tt, state] <- log_emission[first + tt - 1L, state] + value[winner]
          back[tt, state] <- from[winner]
        }
      }
    }
    path[last] <- which.max(delta[len, ])
    if (len > 1L) {
      for (t in last:(first + 1L)) {
        path[t - 1L] <- back[t - first + 1L, path[t]]
      }
    }
  }
  path
}

.ash_hmm_step_summary <- function(path, sequence_id, mu = NULL) {
  n <- length(path)
  segments <- .ash_hmm_segments(sequence_id)
  rows <- vector("list", 0L)
  per_sequence <- data.frame(
    sequence = character(), steps = integer(), changes = integer(),
    stringsAsFactors = FALSE)
  for (s in seq_along(segments$starts)) {
    first <- segments$starts[s]
    last <- segments$ends[s]
    local_change <- if (first < last) {
      which(path[(first + 1L):last] != path[first:(last - 1L)]) + first
    } else integer()
    starts <- c(first, local_change)
    ends <- c(local_change - 1L, last)
    for (j in seq_along(starts)) {
      state <- path[starts[j]]
      rows[[length(rows) + 1L]] <- data.frame(
        sequence = as.character(sequence_id[first]),
        start = starts[j], end = ends[j], length = ends[j] - starts[j] + 1L,
        state = state,
        state_mean = if (is.null(mu)) NA_real_ else mu[state],
        stringsAsFactors = FALSE)
    }
    per_sequence <- rbind(
      per_sequence,
      data.frame(sequence = as.character(sequence_id[first]),
                 steps = length(starts), changes = length(starts) - 1L,
                 stringsAsFactors = FALSE))
  }
  segment_table <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(segments = segment_table,
       per_sequence = per_sequence,
       step_count = sum(per_sequence$steps),
       change_count = sum(per_sequence$changes),
       occupied_state_count = length(unique(path)))
}

# Mark the lower-occupancy member of each very close pair as a possible pruning
# candidate. Actual deletion is accepted only after the caller checks the HMM
# marginal likelihood.
.ash_hmm_close_state_candidates <- function(mu, occupancy, distance,
                                            protected = integer()) {
  if (length(mu) < 2L || !is.finite(distance) || distance <= 0) {
    return(integer())
  }
  ord <- order(mu)
  close <- which(diff(mu[ord]) <= distance)
  if (!length(close)) return(integer())
  drop <- integer()
  for (j in close) {
    pair <- ord[c(j, j + 1L)]
    eligible <- setdiff(pair, protected)
    if (!length(eligible)) next
    if (length(eligible) == 1L) {
      drop <- c(drop, eligible)
    } else {
      drop <- c(drop, eligible[which.min(occupancy[eligible])])
    }
  }
  unique(drop)
}

.ash_hmm_expand_prior <- function(x, nr, nc, name) {
  if (length(x) == 1L) x <- matrix(x, nr, nc)
  else if (length(x) == nc) x <- matrix(rep(x, each = nr), nr, nc)
  else if (is.matrix(x) && all(dim(x) == c(nr, nc))) x <- x
  else stop(sprintf("'%s' must be scalar, length %d, or a %d by %d matrix.",
                    name, nc, nr, nc), call. = FALSE)
  if (anyNA(x) || any(!is.finite(x)) || any(x < 1)) {
    stop(sprintf("Every '%s' entry must be finite and at least one.", name),
         call. = FALSE)
  }
  x
}

.ash_hmm_dirichlet_penalty <- function(probability, alpha, active = NULL) {
  if (is.null(active)) active <- matrix(TRUE, nrow(probability), ncol(probability))
  use <- active & alpha > 1
  if (!any(use)) return(0)
  if (any(probability[use] <= 0)) return(-Inf)
  sum((alpha[use] - 1) * log(probability[use]))
}

.ash_hmm_viterbi <- function(log_emission, A, init_prob, mask, sequence_id) {
  n <- nrow(log_emission)
  m <- ncol(log_emission)
  log_A <- .ash_hmm_safe_log(A)
  log_pi <- .ash_hmm_safe_log(init_prob)
  edge <- which(mask, arr.ind = TRUE)
  incoming <- lapply(seq_len(m), function(j) which(edge[, 2L] == j))
  segments <- .ash_hmm_segments(sequence_id)
  path <- integer(n)

  for (segment in seq_along(segments$starts)) {
    first <- segments$starts[segment]
    last <- segments$ends[segment]
    delta <- matrix(-Inf, last - first + 1L, m)
    back <- matrix(NA_integer_, last - first + 1L, m)
    delta[1L, ] <- log_pi + log_emission[first, ]
    if (first < last) {
      for (tt in 2L:(last - first + 1L)) {
        for (state in seq_len(m)) {
          e <- incoming[[state]]
          from <- edge[e, 1L]
          values <- delta[tt - 1L, from] + log_A[cbind(from, rep(state, length(from)))]
          winner <- which.max(values)
          delta[tt, state] <- log_emission[first + tt - 1L, state] + values[winner]
          back[tt, state] <- from[winner]
        }
      }
    }
    path[last] <- which.max(delta[nrow(delta), ])
    if (first < last) {
      for (t in last:(first + 1L)) {
        path[t - 1L] <- back[t - first + 1L, path[t]]
      }
    }
  }
  path
}

# Mean and variance of a standard normal truncated below at -z. For very
# negative z, direct inverse-Mills calculations suffer catastrophic
# cancellation; the displayed tail expansions retain accuracy.
.ash_hmm_lower_truncated_standard_moments <- function(z) {
  mean_shift <- variance_factor <- numeric(length(z))
  extreme <- z < -10

  if (any(!extreme)) {
    zz <- z[!extreme]
    inverse_mills <- exp(
      stats::dnorm(zz, log = TRUE) - stats::pnorm(zz, log.p = TRUE))
    mean_shift[!extreme] <- inverse_mills
    variance_factor[!extreme] <- pmax(
      0, 1 - zz * inverse_mills - inverse_mills^2)
  }

  if (any(extreme)) {
    a <- -z[extreme]
    inverse_a <- 1 / a
    # lambda(-a) - a and Var{Z | Z >= a}; terms through a^-7/a^-6.
    residual <- inverse_a - 2 * inverse_a^3 +
      10 * inverse_a^5 - 74 * inverse_a^7
    mean_shift[extreme] <- a + residual
    variance_factor[extreme] <- pmax(
      0, inverse_a^2 - 6 * inverse_a^4 + 50 * inverse_a^6)
  }

  list(mean_shift = mean_shift, variance_factor = variance_factor)
}

.ash_hmm_posterior_summary <- function(
    y, se, mu, prior_sd, rho, log_emission, gamma,
    effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  posterior_mean <- posterior_second <- prob_ge_zero <- prob_le_zero <-
    prob_zero <- numeric(n)
  log_rho <- .ash_hmm_safe_log(rho)

  for (state in seq_len(m)) {
    terms <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, effect_support)
    terms <- sweep(terms, 2L, log_rho[state, ], "+")
    component_probability <- exp(terms - log_emission[, state])
    state_mean <- state_second <- state_ge <- state_le <- state_zero <- numeric(n)

    for (component in seq_len(l)) {
      tau <- prior_sd[component]
      if (tau == 0) {
        component_mean <- rep(mu[state], n)
        component_variance <- rep(0, n)
        ge <- rep(as.numeric(mu[state] >= 0), n)
        le <- rep(as.numeric(mu[state] <= 0), n)
        p0 <- rep(as.numeric(mu[state] == 0), n)
      } else {
        component_variance <- tau^2 * se^2 / (tau^2 + se^2)
        component_mean <- (se^2 * mu[state] + tau^2 * y) / (tau^2 + se^2)
        component_sd <- sqrt(component_variance)
        if (effect_support == "nonnegative") {
          z <- component_mean / component_sd
          truncated <- .ash_hmm_lower_truncated_standard_moments(z)
          component_mean <- component_mean +
            component_sd * truncated$mean_shift
          component_variance <- component_variance *
            truncated$variance_factor
          component_mean <- pmax(0, component_mean)
          ge <- rep(1, n)
          le <- rep(0, n)
        } else {
          ge <- stats::pnorm(component_mean / component_sd)
          le <- stats::pnorm(-component_mean / component_sd)
        }
        p0 <- rep(0, n)
      }
      w <- component_probability[, component]
      state_mean <- state_mean + w * component_mean
      state_second <- state_second + w * (component_variance + component_mean^2)
      state_ge <- state_ge + w * ge
      state_le <- state_le + w * le
      state_zero <- state_zero + w * p0
    }

    g <- gamma[, state]
    posterior_mean <- posterior_mean + g * state_mean
    posterior_second <- posterior_second + g * state_second
    prob_ge_zero <- prob_ge_zero + g * state_ge
    prob_le_zero <- prob_le_zero + g * state_le
    prob_zero <- prob_zero + g * state_zero
  }

  if (effect_support == "nonnegative") {
    posterior_mean <- pmax(0, posterior_mean)
    prob_ge_zero <- rep(1, n)
  }
  posterior_variance <- pmax(0, posterior_second - posterior_mean^2)
  list(mean = posterior_mean,
       sd = sqrt(posterior_variance),
       probability_ge_zero = pmin(1, pmax(0, prob_ge_zero)),
       probability_le_zero = pmin(1, pmax(0, prob_le_zero)),
       probability_zero = pmin(1, pmax(0, prob_zero)),
       lfsr = pmin(prob_ge_zero, prob_le_zero))
}




#' @export
print.ash_hmm_fit <- function(x, ...) {
  cat("Adaptive-shrinkage HMM fit\n")
  cat("  observations:", nrow(x$state_probability), "\n")
  cat("  states:", ncol(x$state_probability))
  if (isTRUE(x$grid$automatic)) {
    cat(" (retained from", length(x$grid$full_mu), "automatic candidates)")
  }
  cat("\n")
  if (!is.null(x$fitted$effect_support)) {
    cat("  effect support:", x$fitted$effect_support, "\n")
  }
  if (!is.null(x$inference$method)) {
    cat("  training inference:", x$inference$method)
    if (isTRUE(x$inference$final_exact_used) &&
        identical(x$inference$method, "variational")) {
      cat(" (exact final smoother)")
    }
    cat("\n")
  }
  if (isTRUE(x$fitted$learn_state_means)) {
    moved <- sum(abs(x$fitted$mu - x$fitted$mean_anchor) > 1e-10)
    cat("  learned state means:", moved, "of", length(x$fitted$mu), "\n")
  }
  if (!is.null(x$pruning_history) && nrow(x$pruning_history)) {
    cat("  dynamically pruned states:", nrow(x$pruning_history), "\n")
  }
  if (!is.null(x$model_selection) &&
      isTRUE(x$model_selection$collapsed_to_null)) {
    cat("  global model selection: strict null\n")
  }
  if (!is.null(x$step_selection)) {
    cat("  penalized steps:", x$step_selection$step_count,
        "(changes:", x$step_selection$change_count, ")\n")
  }
  cat("  iterations:", x$iterations,
      if (x$converged) "(converged)" else "(maximum reached)", "\n")
  cat("  log likelihood:", format(x$log_likelihood, digits = 8), "\n")
  if (is.finite(x$variational_bound)) {
    cat("  variational training bound:",
        format(x$variational_bound, digits = 8), "\n")
  }
  invisible(x)
}

# Byte-compile the pure-R hot paths when this file is sourced outside an
# installed byte-compiled package.
if (requireNamespace("compiler", quietly = TRUE)) {
  .ash_hmm_forward_backward_scaled <- compiler::cmpfun(
    .ash_hmm_forward_backward_scaled)
  .ash_hmm_forward_backward_log <- compiler::cmpfun(
    .ash_hmm_forward_backward_log)
  .ash_hmm_variational_states <- compiler::cmpfun(
    .ash_hmm_variational_states)
  .ash_hmm_component_statistics <- compiler::cmpfun(
    .ash_hmm_component_statistics)
  .ash_hmm_component_statistics_cached <- compiler::cmpfun(
    .ash_hmm_component_statistics_cached)
  .ash_hmm_posterior_summary <- compiler::cmpfun(
    .ash_hmm_posterior_summary)
}
