% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_on_susiF_obj.R
\name{update_cal_cs}
\alias{update_cal_cs}
\alias{update_cal_cs.susiF}
\title{Update susiF by computing credible sets}
\usage{
update_cal_cs(obj, cov_lev = 0.95, ...)

\method{update_cal_cs}{susiF}(obj, cov_lev = 0.95, l, ...)
}
\arguments{
\item{obj}{a susiF object defined by  init_susiF_obj  function}

\item{cov_lev}{numeric between 0 and 1, corresponding to the expected level of coverage of the cs if not specified set to 0.95}
}
\value{
susiF object
}
\description{
Update susiF by computing credible sets
}
\keyword{internal}
