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
#'   the same samples used to compute \code{Bhat} and \code{Shat}. Using
#'   an out-of-sample matrix may produce unreliable results.
#'
#' @param N The sample size.
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
#' @param wav_trans logical, if true the algorithm will consider that the summary statistics based on wavelet transformed data (\code{Bhat} and \code{Shat}).
#'  if False, the algorithm will rescale \code{Bhat} and \code{Shat} to obtain summary statistics from wavelet regression. Default set as FALSE .
#'
make_data_suff_stat <- function(Bhat , Shat , R , N  , var_y , XtX , Xty , yty , wav_trans=  FALSE)
{


  if(wav_trans){
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
  }else{


    wav_suff_stat <- trans_sum_stat_wreg(Bhat,Shat)


    out <-  list( Bhat=wav_suff_stat$Bhat,
                  Shat=wav_suff_stat$Shat,
                  R=R,
                  N=N,
                  var_y=var_y,
                  XtX=XtX,
                  Xty=Xty,
                  yty=yty,
                  exp_residual= Xty,
                  part_exp_residual= Xty)
  }

  class(out)  <- c("suff_stat","list")
  return(out)
}

#' @title Transform raw regression coefficient to wavelet regression coefficients
#'
#' @param Bhat matrix of regression coefficient
#'
#' @param Shat matrix of standard error
#'
#'
#' @return list of two
#'
#' \item{Bhat}{ matrix pxJ wavelet  regression coefficient}
#' \item{Shat}{ matrix pxJ  wavelet  regression coefficient standard error}
#'@export

trans_sum_stat_wreg <- function( Bhat, Shat){

  W1 <- (GenW(n=ncol(Bhat)  , filter.number = 10, family = "DaubLeAsymm"))
  wBhat <-  Bhat_recov(Bhat)
  wShat <-  Shat_recov(Shat,W1)
  out  <- list( Bhat = wBhat,
                Shat = wShat)
  return(out)
}


#' @title Transform raw regression coefficient to wavelet regression coefficients
#'
#' @param Bhat matrix of regression coefficient
#'
#' @param W1 matrix associated with a wavelet transform
#'
#' @return Wavelet regression coefficient using wavelet transform from W1
#'@export
Bhat_recov <- function(Bhat, W1)
{
  W <- DWT2(Bhat )
  t_Bhat <-   cbind( W$D,W$C)
  return(t_Bhat)
}


#' @title Transform raw regression coefficient standard error to wavelet regression coefficients standard error
#'
#' @param Shat matrix of standard error
#'
#' @param W1 matrix associated with a wavelet transform
#'
#' @return Matrix of standard error of the wavelet regression coefficient  using wavelet transform from W1
#'@export
Shat_recov <- function(Shat, W1)
{
  t_Shat<-  t(apply(Shat,1,function(x)  Shatwc_recov_vec(x,W1)  ))

  return(t_Shat)
}






#' @title Transform a vector of raw regression coefficient standard error to wavelet regression coefficients standard error
#'
#' @param sd_vec matrix of standard error (must be of power of two)
#'
#' @param W1 matrix associated with a wavelet transform
#'
#' @return Matrix of standard error of the wavelet regression coefficient  using wavelet transform from W1
#'@export


Shatwc_recov_vec <- function(sd_vec,W1)
{


  var_wc <-  c()
  for ( i in 1:length(sd_vec))
  {
    var_wc  <- c(  var_wc ,sum((W1[,i]^2)*diag(sd_vec^2)))
  }
  tt <- sqrt(var_wc)
  out<- c(shifter(tt[-1],-1),tt[1])
  return(out )
}

