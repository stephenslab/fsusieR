% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM.R
\name{L_mixsq}
\alias{L_mixsq}
\alias{L_mixsq.mixture_normal}
\alias{L_mixsq.mixture_normal_per_scale}
\title{Compute likelihood matrix for mixsqp}
\usage{
L_mixsq(G_prior, Bhat, Shat, indx_lst, ...)

\method{L_mixsq}{mixture_normal}(G_prior, Bhat, Shat, indx_lst, is.EBmvFR = FALSE, ...)

\method{L_mixsq}{mixture_normal_per_scale}(G_prior, Bhat, Shat, indx_lst, is.EBmvFR = FALSE, ...)
}
\arguments{
\item{G_prior}{mixture normal prior}

\item{Bhat}{matrix pxJ regression coefficient, Bhat[j,t] corresponds to regression coefficient of Y[,t] on X[,j]}

\item{Shat}{matrix pxJ standard error, Shat[j,t] corresponds to standard error of the regression coefficient of Y[,t] on X[,j]}

\item{indx_lst}{list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class  mixture_normal_per_scale}

\item{\dots}{Other arguments.}
}
\value{
See L argument mixsqp package mixsqp function
}
\description{
Compute likelihood matrix for mixsqp
}
\keyword{internal}
