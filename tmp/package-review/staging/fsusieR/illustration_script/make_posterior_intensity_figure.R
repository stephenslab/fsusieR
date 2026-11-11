# Publication figure: recovery of the structured log-intensity
#
# Run this script from the fsusieR package root. It writes:
#   plots/posterior_intensity/posterior_intensity_calibration.pdf
#   plots/posterior_intensity/posterior_intensity_calibration.png
#   plots/posterior_intensity/posterior_intensity_calibration_metrics.csv
#   plots/posterior_intensity/posterior_intensity_calibration_caption.txt

options(stringsAsFactors = FALSE)

if (!requireNamespace("fsusieR", quietly = TRUE)) {
  stop("Install fsusieR before running this script.")
}
if (!requireNamespace("susieR", quietly = TRUE)) {
  stop("Install susieR to access the N3finemapping example data.")
}

## Configuration ---------------------------------------------------------

sigma_grid <- c(0, 0.5, 1, 2)

N <- 100L
P_requested <- 100L
T_grid <- 2^7
L_true <- 2L
L_fit <- 3L

# The original code added 0.3 once for each of two true effects. Writing the
# resulting intercept explicitly makes the data-generating model unambiguous.
intercept <- 0.6

seed_design <- 1L
seed_counts <- 20260826L

maxit_gaussian <- 100L
maxit_outer <- 100L
maxit_inner <- 100L
sigma2_subcycles <- 10L
fit_tolerance <- 1e-6

# A common y-axis is used in all panels. Very extreme observed counts are
# clipped from display only; calibration metrics always use every value.
observed_display_quantile <- 0.995

# Keep the outputs in the package-level plots folder, including when the script
# is run from inside illustration_script.
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project_root) == "illustration_script") {
  project_root <- dirname(project_root)
}
if (!file.exists(file.path(project_root, "DESCRIPTION"))) {
  stop("Run this script from the fsusieR package root or illustration_script.")
}

output_dir <- file.path(project_root, "plots", "posterior_intensity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_pdf <- file.path(output_dir, "posterior_intensity_calibration.pdf")
figure_png <- file.path(output_dir, "posterior_intensity_calibration.png")
metrics_csv <- file.path(
  output_dir,
  "posterior_intensity_calibration_metrics.csv"
)
caption_txt <- file.path(
  output_dir,
  "posterior_intensity_calibration_caption.txt"
)

## Simulate one shared structured intensity ------------------------------

set.seed(seed_design)
data("N3finemapping", package = "susieR")

X_raw <- N3finemapping$X
if (nrow(X_raw) < N) {
  stop("N3finemapping$X has fewer rows than requested N.")
}

P_initial <- min(P_requested, ncol(X_raw))
sampled_rows <- sample.int(nrow(X_raw), size = N, replace = FALSE)
genotype <- X_raw[
  sampled_rows,
  seq_len(P_initial),
  drop = FALSE
]

nonconstant <- apply(genotype, 2L, stats::var) >= 1e-15
genotype <- genotype[, nonconstant, drop = FALSE]
if (ncol(genotype) < L_true) {
  stop("Too few nonconstant genotype columns to select the true effects.")
}

scale_denominator <- 0.5 * max(genotype)
if (!is.finite(scale_denominator) || scale_denominator == 0) {
  stop("The genotype scaling denominator is not finite and nonzero.")
}

# This preserves the scaling used in the exploratory script.
X <- (genotype - 0.99 * min(genotype)) / scale_denominator
G <- X

true_positions <- sample.int(ncol(G), size = L_true, replace = FALSE)
G[, true_positions] <- sweep(
  G[, true_positions, drop = FALSE],
  2L,
  apply(G[, true_positions, drop = FALSE], 2L, min),
  FUN = "-"
)

true_effects <- matrix(0, nrow = L_true, ncol = T_grid)
true_effects[1L, 10:70] <- 2.5 * sin((10:70) / 50)
true_effects[2L, 80:120] <- 2.5 * sin((80:120) / 50)

B_true <- matrix(intercept, nrow = N, ncol = T_grid)
for (l in seq_len(L_true)) {
  B_true <- B_true + tcrossprod(G[, true_positions[l]], true_effects[l, ])
}

# Reusing Z makes B_true and the standardized latent deviations identical in
# every scenario. Only the residual standard deviation sigma changes.
Z <- matrix(stats::rnorm(N * T_grid), nrow = N, ncol = T_grid)

## Fit both methods and collect calibration statistics ------------------

mixsqp_control <- list(
  verbose = FALSE,
  eps = 1e-6,
  numiter.em = 40L
)

calibration_metrics <- function(truth,
                                estimate,
                                sigma,
                                method,
                                converged = NA,
                                iterations = NA_integer_,
                                estimated_sigma2 = NA_real_) {
  truth <- as.numeric(truth)
  estimate <- as.numeric(estimate)

  calibration_fit <- stats::lm(estimate ~ truth)
  coefficients <- stats::coef(calibration_fit)

  data.frame(
    sigma = sigma,
    method = method,
    intercept = unname(coefficients[[1L]]),
    slope = unname(coefficients[[2L]]),
    bias = mean(estimate - truth),
    rmse = sqrt(mean((estimate - truth)^2)),
    mae = mean(abs(estimate - truth)),
    correlation = stats::cor(estimate, truth),
    n_values = length(truth),
    converged = converged,
    iterations = iterations,
    estimated_sigma2 = estimated_sigma2,
    check.names = FALSE
  )
}

scenario_results <- vector("list", length(sigma_grid))
metric_results <- vector("list", 2L * length(sigma_grid))

for (k in seq_along(sigma_grid)) {
  sigma <- sigma_grid[k]
  message("Fitting latent log-intensity SD sigma = ", sigma)

  latent_log_intensity <- B_true + sigma * Z
  lambda <- exp(latent_log_intensity)

  # Resetting this seed supplies a common conditional Poisson random-number
  # stream across scenarios and makes results independent of loop order.
  set.seed(seed_counts)
  Y <- matrix(
    stats::rpois(length(lambda), lambda = as.numeric(lambda)),
    nrow = N,
    ncol = T_grid
  )

  gaussian_fit <- fsusieR::susiF(
    Y = log1p(Y),
    X = X,
    L = L_fit,
    maxit = maxit_gaussian,
    tol = fit_tolerance,
    control_mixsqp = mixsqp_control,
    verbose = FALSE
  )

  poisson_fit <- fsusieR::Pois_fSuSiE(
    Y = Y,
    X = X,
    L = L_fit,
    maxit_outer = maxit_outer,
    maxit_inner = maxit_inner,
    tol = fit_tolerance,
    control_mixsqp = mixsqp_control,
    sigma2_subcycles = sigma2_subcycles,
    verbose = FALSE,
    diagnostic_plot = FALSE
  )

  if (!isTRUE(poisson_fit$converged)) {
    warning(
      "Poisson fSuSiE did not converge for sigma = ",
      sigma,
      " after ",
      poisson_fit$n_iter,
      " iterations."
    )
  }

  poisson_estimate <- as.matrix(poisson_fit$B_pm)
  gaussian_estimate <- as.matrix(gaussian_fit$ind_fitted_func)

  if (!identical(dim(poisson_estimate), dim(B_true)) ||
      !identical(dim(gaussian_estimate), dim(B_true))) {
    stop("At least one fitted intensity has the wrong dimensions.")
  }

  scenario_results[[k]] <- list(
    sigma = sigma,
    observed = log1p(Y),
    poisson = poisson_estimate,
    gaussian_log1p = gaussian_estimate
  )

  metric_results[[2L * k - 1L]] <- calibration_metrics(
    truth = B_true,
    estimate = poisson_estimate,
    sigma = sigma,
    method = "Poisson fSuSiE",
    converged = poisson_fit$converged,
    iterations = poisson_fit$n_iter,
    estimated_sigma2 = poisson_fit$sigma2
  )

  metric_results[[2L * k]] <- calibration_metrics(
    truth = B_true,
    estimate = gaussian_estimate,
    sigma = sigma,
    method = "fSuSiE on log1p counts"
  )
}

metrics <- do.call(rbind, metric_results)
metrics$method <- factor(
  metrics$method,
  levels = c("Poisson fSuSiE", "fSuSiE on log1p counts")
)
metrics <- metrics[order(metrics$sigma, metrics$method), ]
row.names(metrics) <- NULL

utils::write.csv(metrics, metrics_csv, row.names = FALSE)

## Publication figure ----------------------------------------------------

method_colors <- c(
  "Poisson fSuSiE" = "#0072B2",
  "fSuSiE on log1p counts" = "#D55E00"
)
observed_color <- grDevices::adjustcolor("#333333", alpha.f = 0.20)
poisson_point_color <- grDevices::adjustcolor(
  method_colors[["Poisson fSuSiE"]],
  alpha.f = 0.58
)
gaussian_point_color <- grDevices::adjustcolor(
  method_colors[["fSuSiE on log1p counts"]],
  alpha.f = 0.58
)

truth_vector <- as.numeric(B_true)
all_observed <- unlist(
  lapply(scenario_results, function(result) as.numeric(result$observed)),
  use.names = FALSE
)
all_estimates <- unlist(
  lapply(
    scenario_results,
    function(result) c(
      as.numeric(result$poisson),
      as.numeric(result$gaussian_log1p)
    )
  ),
  use.names = FALSE
)

x_padding <- 0.04 * diff(range(truth_vector))
x_limits <- range(truth_vector) + c(-x_padding, x_padding)

observed_upper <- as.numeric(stats::quantile(
  all_observed,
  probs = observed_display_quantile,
  names = FALSE
))
y_limits <- c(
  min(-0.25, all_estimates, truth_vector),
  max(observed_upper, all_estimates, truth_vector)
)
y_padding <- 0.025 * diff(y_limits)
y_limits <- y_limits + c(-y_padding, y_padding)

displayed_observation_fraction <- mean(
  all_observed >= y_limits[1L] & all_observed <= y_limits[2L]
)

format_sigma <- function(x) format(x, trim = TRUE, scientific = FALSE)

draw_figure <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::layout(
    matrix(c(1:4, rep(5, 4)), nrow = 2L, byrow = TRUE),
    heights = c(1, 0.16)
  )
  graphics::par(
    oma = c(3.5, 4.0, 0.2, 0.3),
    family = "sans",
    las = 1,
    xaxs = "i",
    yaxs = "i"
  )

  for (k in seq_along(scenario_results)) {
    result <- scenario_results[[k]]
    sigma <- result$sigma

    panel_metrics <- metrics[metrics$sigma == sigma, , drop = FALSE]
    poisson_metrics <- panel_metrics[
      panel_metrics$method == "Poisson fSuSiE",
      ,
      drop = FALSE
    ]
    gaussian_metrics <- panel_metrics[
      panel_metrics$method == "fSuSiE on log1p counts",
      ,
      drop = FALSE
    ]

    panel_title <- sprintf(
      "(%s)   \u03c3 = %s",
      letters[k],
      format_sigma(sigma)
    )

    graphics::par(mar = c(2.3, if (k == 1L) 2.8 else 1.0, 2.4, 0.6))
    graphics::plot(
      truth_vector,
      as.numeric(result$observed),
      type = "n",
      xlim = x_limits,
      ylim = y_limits,
      xlab = "",
      ylab = "",
      axes = FALSE,
      main = panel_title,
      cex.main = 0.96,
      font.main = 2
    )

    x_ticks <- graphics::axTicks(1L)
    y_ticks <- graphics::axTicks(2L)
    graphics::abline(
      v = x_ticks,
      h = y_ticks,
      col = "#E6E6E6",
      lwd = 0.65
    )

    graphics::points(
      truth_vector,
      as.numeric(result$observed),
      pch = 16,
      cex = 0.42,
      col = observed_color
    )
    graphics::abline(a = 0, b = 1, col = "#222222", lty = 2, lwd = 1.0)

    graphics::points(
      truth_vector,
      as.numeric(result$poisson),
      pch = 16,
      cex = 0.62,
      col = poisson_point_color
    )
    graphics::points(
      truth_vector,
      as.numeric(result$gaussian_log1p),
      pch = 17,
      cex = 0.66,
      col = gaussian_point_color
    )

    graphics::abline(
      a = poisson_metrics$intercept,
      b = poisson_metrics$slope,
      col = method_colors[["Poisson fSuSiE"]],
      lwd = 1.35
    )
    graphics::abline(
      a = gaussian_metrics$intercept,
      b = gaussian_metrics$slope,
      col = method_colors[["fSuSiE on log1p counts"]],
      lwd = 1.35
    )

    annotation <- c(
      sprintf("%-8s %9s | %5s | %4s", "", "intercept", "slope", "RMSE"),
      sprintf(
        "%-8s %9.2f | %5.2f | %4.2f",
        "Poisson",
        poisson_metrics$intercept,
        poisson_metrics$slope,
        poisson_metrics$rmse
      ),
      sprintf(
        "%-8s %9.2f | %5.2f | %4.2f",
        "log1p",
        gaussian_metrics$intercept,
        gaussian_metrics$slope,
        gaussian_metrics$rmse
      )
    )

    annotation_x <- x_limits[1L] + 0.025 * diff(x_limits)
    annotation_y <- y_limits[2L] - c(0.035, 0.095, 0.155) * diff(y_limits)
    graphics::rect(
      xleft = annotation_x - 0.008 * diff(x_limits),
      ybottom = annotation_y[3L] - 0.035 * diff(y_limits),
      xright = annotation_x + 0.73 * diff(x_limits),
      ytop = annotation_y[1L] + 0.025 * diff(y_limits),
      border = NA,
      col = grDevices::adjustcolor("white", alpha.f = 0.88)
    )
    graphics::text(
      x = rep(annotation_x, 3L),
      y = annotation_y,
      labels = annotation,
      adj = c(0, 0.5),
      family = "mono",
      font = c(2, 1, 1),
      col = c(
        "#222222",
        method_colors[["Poisson fSuSiE"]],
        method_colors[["fSuSiE on log1p counts"]]
      ),
      cex = 0.60
    )

    graphics::axis(1L, at = x_ticks, cex.axis = 0.78, tck = -0.018)
    graphics::axis(
      2L,
      at = y_ticks,
      labels = if (k == 1L) y_ticks else FALSE,
      cex.axis = 0.78,
      tck = -0.018
    )
    graphics::box(col = "#444444", lwd = 0.8)
  }

  # Shared legend below the four panels.
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = c(
      "Observed log1p count",
      "Poisson fSuSiE",
      "fSuSiE on log1p counts",
      "Identity"
    ),
    col = c(
      grDevices::adjustcolor("#333333", alpha.f = 0.65),
      method_colors[["Poisson fSuSiE"]],
      method_colors[["fSuSiE on log1p counts"]],
      "#222222"
    ),
    pch = c(16, 16, 17, NA),
    lty = c(NA, 1, 1, 2),
    lwd = c(NA, 1.35, 1.35, 1.0),
    horiz = TRUE,
    bty = "n",
    cex = 0.86,
    pt.cex = c(1.05, 1.15, 1.20, NA),
    seg.len = 2.2,
    x.intersp = 0.7
  )

  graphics::mtext(
    "True structured log-intensity",
    side = 1,
    outer = TRUE,
    line = 2.05,
    cex = 0.94
  )
  graphics::mtext(
    "Observed or estimated log-intensity",
    side = 2,
    outer = TRUE,
    line = 2.55,
    las = 0,
    cex = 0.92
  )
}

if (isTRUE(capabilities("cairo"))) {
  grDevices::cairo_pdf(
    filename = figure_pdf,
    width = 14.50,
    height = 4.35,
    family = "Arial",
    symbolfamily = "Arial",
    pointsize = 10
  )
} else {
  grDevices::pdf(
    file = figure_pdf,
    width = 14.50,
    height = 4.35,
    family = "sans",
    useDingbats = FALSE,
    pointsize = 10
  )
}
draw_figure()
grDevices::dev.off()

png_arguments <- list(
  filename = figure_png,
  width = 14.50,
  height = 4.35,
  units = "in",
  res = 320,
  pointsize = 10,
  bg = "white"
)
if (isTRUE(capabilities("cairo"))) {
  png_arguments$type <- "cairo"
}
do.call(grDevices::png, png_arguments)
draw_figure()
grDevices::dev.off()

## Data-driven figure caption --------------------------------------------

sigma_text <- paste(vapply(sigma_grid, format_sigma, character(1L)), collapse = ", ")

lowest_sigma <- min(sigma_grid)
highest_sigma <- max(sigma_grid)

get_metric <- function(sigma, method, variable) {
  metrics[
    metrics$sigma == sigma & metrics$method == method,
    variable,
    drop = TRUE
  ]
}

calibration_sentence <- sprintf(
  paste0(
    "From sigma = %s to sigma = %s, the calibration slope changed from ",
    "%.2f to %.2f for Poisson fSuSiE and from %.2f to %.2f for fSuSiE ",
    "on log1p counts; the corresponding calibration intercepts changed ",
    "from %.2f to %.2f and from %.2f to %.2f, respectively."
  ),
  format_sigma(lowest_sigma),
  format_sigma(highest_sigma),
  get_metric(lowest_sigma, "Poisson fSuSiE", "slope"),
  get_metric(highest_sigma, "Poisson fSuSiE", "slope"),
  get_metric(lowest_sigma, "fSuSiE on log1p counts", "slope"),
  get_metric(highest_sigma, "fSuSiE on log1p counts", "slope"),
  get_metric(lowest_sigma, "Poisson fSuSiE", "intercept"),
  get_metric(highest_sigma, "Poisson fSuSiE", "intercept"),
  get_metric(lowest_sigma, "fSuSiE on log1p counts", "intercept"),
  get_metric(highest_sigma, "fSuSiE on log1p counts", "intercept")
)

display_sentence <- if (displayed_observation_fraction < 0.999999) {
  sprintf(
    paste0(
      "The common plotting range contains %.1f%% of observed log1p counts; ",
      "clipping is for display only, and all values were used to compute ",
      "the calibration statistics."
    ),
    100 * displayed_observation_fraction
  )
} else {
  "All observed log1p counts fall within the common plotting range."
}

caption <- paste(
  sprintf(
    paste0(
      "Recovery of the structured log-intensity under increasing latent ",
      "overdispersion. Counts were generated for %d individuals at %d grid ",
      "points according to Y_it | eta_it ~ Poisson[exp(eta_it)], with ",
      "eta_it = B_it + sigma*z_it, z_it ~ N(0,1), and sigma in {%s}. ",
      "The same genotype matrix, %d causal variables, structured ",
      "log-intensity B, and standardized latent deviations z were used in ",
      "every panel."
    ),
    N,
    T_grid,
    sigma_text,
    L_true
  ),
  paste0(
    "Gray points are the observed log1p counts, blue circles are posterior ",
    "means of B from Poisson fSuSiE, and orange triangles are fitted values ",
    "from Gaussian fSuSiE applied to log1p counts. The dashed line is the ",
    "identity line; colored solid lines are least-squares calibration lines. ",
    "Panel annotations report the calibration intercept, slope, and root ",
    "mean squared error, computed over all individual-grid pairs."
  ),
  calibration_sentence,
  display_sentence
)

writeLines(enc2utf8(caption), con = caption_txt, useBytes = TRUE)

message("Wrote figure:  ", normalizePath(figure_pdf, winslash = "/"))
message("Wrote preview: ", normalizePath(figure_png, winslash = "/"))
message("Wrote metrics: ", normalizePath(metrics_csv, winslash = "/"))
message("Wrote caption: ", normalizePath(caption_txt, winslash = "/"))
message("\nFigure caption:\n", caption)
