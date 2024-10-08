% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_susiF}
\alias{plot_susiF}
\alias{plot.susiF}
\alias{plot_susiF_pip}
\alias{plot_susiF_effect}
\title{fSuSiE Plots}
\usage{
plot_susiF(
  obj,
  which_plot = c("both", "pip", "effect"),
  effect = "all",
  cred_band = TRUE,
  show_affected_region = TRUE,
  show_outing_grid = ifelse(diff(range(diff(obj$outing_grid))) < 1e-06, FALSE, TRUE),
  lfsr_curve = TRUE,
  line_width = 0.35,
  point_size = 1.25,
  dot_size = 0.5,
  pos_SNP,
  point_shape,
  font_size = 10,
  title = "",
  ...
)

\method{plot}{susiF}(x, ...)

plot_susiF_pip(
  obj,
  title = "",
  pos_SNP,
  point_shape,
  point_size = 1.25,
  font_size = 10
)

plot_susiF_effect(
  obj,
  effect = "all",
  title = "",
  cred_band = TRUE,
  show_affected_region = TRUE,
  show_outing_grid = ifelse(diff(range(diff(obj$outing_grid))) < 1e-06, FALSE, TRUE),
  lfsr_curve = TRUE,
  line_width = 0.35,
  dot_size = 0.5,
  font_size = 10
)
}
\arguments{
\item{obj}{Output of the susiF function.}

\item{which_plot}{Which plots to return; a PIP plot, effect plot,
or both (in which case the return value is a list containing the
two plots.}

\item{effect}{The indices of the effects to be plotted, or use
\code{effect = "all"} to plot all effects.}

\item{cred_band}{logical. If \code{TRUE}, plot credible bands if
the fSuSiE model was fitted with wavelet regression.}

\item{show_affected_region}{If \code{show_affected_region = TRUE},
the regions in which the credible bands cross zero are also shown.}

\item{show_outing_grid}{If \code{show_outing_grid = TRUE}, show
the grid positions at which the effects were estimated. By default,
this option is set to \code{TRUE} only when the grid poisitions are
uneven.}

\item{lfsr_curve}{Logical. If \code{TRUE}, plot estimated lfsr of the
effect at each base pair if obj fitted with HMM regression. This
has no effect unless the \code{\link{susiF}} option
\code{post_processing = "HMM"} was used.}

\item{line_width}{Numeric. Width of the plotted lines.}

\item{point_size}{numeric, size of the points.}

\item{dot_size}{numeric, size of the points in the effect plot.}

\item{pos_SNP}{vector, containing the base pair of the SNPs}

\item{point_shape}{vector, containing the shape of dots}

\item{font_size}{Passed as the \dQuote{ont_size} argument to
\code{\link[cowplot]{theme_cowplot}}.}

\item{title}{The title of the plot.}

\item{x}{Output of the susiF function.}

\item{ldots}{additional arguments}
}
\description{
Various visualizations of fSuSiE results.
}
