############# Define custom data class
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @title Define sufficient stat data format
#'
#' @param Bhat A p by t matrix of estimated effects.
#'
#' @param Shat A p by t matrix of standard errors.
#'
#' @param R A p by p correlation matrix. It should be estimated from
#'   the same samples used to compute \code{bhat} and \code{shat}. Using
#'   an out-of-sample matrix may produce unreliable results.
#'
#'  @param N The sample size.
#'
#' @param var_y a T vector of The sample variance of Y at each time point , defined as \eqn{y'y/(n-1)}.
#'   When the sample variance cannot be provided, the coefficients
#'   (returned from \code{coef}) are computed on the "standardized" X, y
#'   scale.
#'
#' @param XtX A p by p matrix \eqn{X'X} in which the columns of X
#'   are centered to have mean zero.
#'
#' @param Xty A p by T matrix \eqn{X'y} in which y and the columns of X are
#'   centered to have mean zero.
#'
#' @param yty A T by T scalar \eqn{y'y} in which y is centered to have mean
#'   zero.
#'
make_data_suff_stat <- function(Bhat , Shat , R , N  , var_y , XtX , Xty , yty )
{
  out <-  list( Bhat=Bhat,
                Shat=Shat,
                R=R,
                N=N,
                var_y=var_y,
                XtX=XtX,
                Xty=Xty,
                yty=yty,
                exp_residual= Xty,
                part_exp_residual= Xty)
  class(out)  <- c("suff_stat","list")
  return(out)
}
