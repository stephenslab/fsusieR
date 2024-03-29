% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_on_susiF_obj.R
\name{plot_susiF}
\alias{plot_susiF}
\title{Plot susiF object}
\usage{
plot_susiF(
  obj,
  title = "",
  cred.band = TRUE,
  lfsr.curve = TRUE,
  size_line = 2,
  size_point = 4,
  pos_SNP,
  pip_only = FALSE,
  point_shape,
  ...
)
}
\arguments{
\item{obj}{output of the susiF function}

\item{title}{character}

\item{cred.band}{logical, if TRUE, plot credible bands if obj fitted with wavelets regression. Set as TRUE by default}

\item{lfsr.curve}{logical, if TRUE, plot estimated lfsr of the effect at each base pair  if obj fitted with HMM regression. Set as TRUE by default}

\item{size_line}{numeric, width of the plotted lines}

\item{size_point}{numeric, size of the point}

\item{pos_SNP}{vector, containing the base pair of the SNPs}

\item{pip_only}{logical, if TRUE only ouput the PIP plot}

\item{point_shape}{vector, containing the shape of dots}

\item{\dots}{Other arguments..}
}
\description{
Plot susiF object
}
