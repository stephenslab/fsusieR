% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_on_EBmvFR_obj.R,
%   R/operation_on_susiF_obj.R
\name{out_prep.EBmvFR}
\alias{out_prep.EBmvFR}
\alias{out_prep}
\alias{out_prep.susiF}
\title{Preparing output of main susiF function}
\usage{
\method{out_prep}{EBmvFR}(obj, Y, X, indx_lst, outing_grid, ...)

out_prep(obj, Y, X, indx_lst, outing_grid, ...)

\method{out_prep}{susiF}(
  obj,
  Y,
  X,
  indx_lst,
  outing_grid,
  filter_cs,
  filter.number = 10,
  family = "DaubLeAsymm",
  post_processing = "TI",
  tidx = NULL,
  names_colX = NULL,
  pos,
  ...
)
}
\arguments{
\item{obj}{a susiF object defined by init_susiF_obj function}

\item{Y}{functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelets that assume that J is of the form J^2. If J is not a power of 2, susiF internally remaps the data into a grid of length 2^J}

\item{X}{matrix of size N by p}

\item{indx_lst}{list generated by gen_wavelet_indx for the given level of resolution}

\item{outing_grid}{grid use to fit fsusie}

\item{filter_cs}{logical, if TRUE filter the credible set (removing low purity cs and cs with estimated prior equal to 0)}

\item{pos}{the original position of the Y column}
}
\value{
susiF object
}
\description{
Preparing output of main susiF function
}
\keyword{internal}
