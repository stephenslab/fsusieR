% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ELBO.R
\name{get_objective}
\alias{get_objective}
\alias{get_objective.susiF}
\title{Get objective function from data and susiF object}
\usage{
get_objective(obj, Y, X, D, C, indx_lst, ...)

\method{get_objective}{susiF}(obj, Y, X, D, C, indx_lst, ...)
}
\arguments{
\item{obj}{a susiF object defined by  init_susiF_obj  function}

\item{Y}{Matrix of outcomes}

\item{X}{matrix of covariate}

\item{D}{matrix of wavelet D coefficients from the original input data (Y)}

\item{C}{vector of wavelet scaling coefficient from the original input data (Y)}

\item{indx_lst}{list generated by gen_wavelet_indx for the given level of resolution}
}
\value{
objective function value
}
\description{
Get objective function from data and susiF object
}
\keyword{internal}
