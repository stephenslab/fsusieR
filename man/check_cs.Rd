% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_on_susiF_obj.R
\name{check_cs}
\alias{check_cs}
\alias{check_cs.susiF}
\title{Check purity credible sets}
\usage{
check_cs(obj, min_purity = 0.5, X, ...)

\method{check_cs}{susiF}(obj, min_purity = 0.5, X, ...)
}
\arguments{
\item{obj}{a susif object defined by init_susiF_obj function}

\item{min_purity}{minimal purity within a CS}

\item{X}{matrix of covariates}
}
\value{
a obj without "dummy" credible s
}
\description{
Check purity credible sets
}
\keyword{internal}
