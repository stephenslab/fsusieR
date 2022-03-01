#'
#' @title Initialize a susiF object using regression coefficients
#'
#' @param L number of non zero coefficients An L-vector containing the indices of the
#'   nonzero coefficients.
#'
#' @param G_prior prior object defined by init_prior function
#'
#' @param data  an object of the class suff_stat define by function \code{\link{make_data_suff_stat}}
#'
#' @export
#' @return A list with elements
#' \item{fitted_wc}{ list of length L, each element contains the fitted wavelet coefficients of effect l}
#' \item{fitted_wc2}{list of length L, each element contains the variance of the fitted wavelet coefficients of effect l}
#' \item{alpha_hist}{ history of the fitted alpha value}
#' \item{N}{ number of indidivual in the study}
#' \item{sigma2}{residual variance}
#' \item{n_wac}{number of wavelet coefficients}
#' \item{ind_fitted_func}{fitted curves of each individual }
#' \item{cs}{credible set}
#' \item{pip}{Posterior inclusion probabilites}
#' \item{G_prior}{a G_prior of the same class as the input G_prior, used for internal calculation}
#'
init_susiF_ss_obj <- function(L, G_prior, data )
{


  fitted_wc       <-  list()
  fitted_wc2      <-  list()
  alpha           <-  list()
  alpha_hist      <-  list()
  ind_fitted_func <-  matrix(NA, nrow = dim(Y)[1], ncol=dim(Y)[2]  )
  cs              <-  list()
  pip             <-  rep(0, dim(X)[2])
  est_pi          <-  list()
  est_sd          <-  list()
  L               <-  L
  G_prior         <-  G_prior
  N               <- data$N
  n_wac           <- ncol(data$Bhat)
  P               <- nrow(data$Bhat)
  sigma2          <- 1
  lBF             <- list()
  KL              <- rep(NA,L)
  ELBO            <- c()
  for ( l in 1:L )
  {
    fitted_wc[[l]]        <-  matrix(NA, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    fitted_wc2[[l]]       <-  matrix(NA, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    alpha [[l]]           <-  rep(NA, dim(X)[2])
    cs[[l]]               <-  list()
    est_pi [[l]]          <-  get_pi_G_prior(G_prior)
    est_sd [[l]]          <-  get_sd_G_prior(G_prior)
    lBF[[l]]              <- rep(NA, ncol(X))
  }
  obj <- list( fitted_wc       = fitted_wc,
               fitted_wc2      = fitted_wc2,
               lBF             = lBF,
               KL              = KL,
               ELBO            = ELBO,
               ind_fitted_func = ind_fitted_func,
               G_prior         = G_prior,
               alpha_hist      = alpha_hist,
               N               = N,
               n_wac           = n_wac,
               sigma2          = sigma2,
               P               = P,
               alpha           = alpha,
               cs              = cs,
               pip             = pip,
               est_pi          = est_pi,
               est_sd          = est_sd,
               L               = L)

  class(obj) <- "susiF_ss"
  return(obj)
}




#' @rdname get_pi
#'
#' @method get_pi susif_ss
#'
#' @export get_pi.susif_ss
#'
#' @export
#'
get_pi.susiF_ss <- function(susiF_ss.obj, l, ...)
{

  if( l >  length(susiF_ss.obj$est_pi))
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  out <- susiF_ss.obj$est_pi[[l]]
  return(out)
}




#' @rdname get_pi
#'
#' @method get_pi susiF_ss
#'
#' @export get_pi.susiF_ss
#'
#' @export
#'
get_lBF.susiF_ss <- function(susiF_ss.obj, l, ...)
{

  if( l >   susiF_ss.obj$L)
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  out <- susiF_ss.obj$lBF[[l]]
  return(out)
}

#' @rdname get_G_prior
#'
#' @method get_G_prior susiF_ss
#'
#' @export get_G_prior.susiF_ss
#'
#' @export
#'
get_G_prior.susiF_ss <- function(susiF_ss.obj, ...)
{
  out <- susiF_ss.obj$G_prior
  return(out)
}




#' @rdname update_pi
#'
#' @method update_pi susiF_ss
#'
#' @export update_pi.susiF_ss
#'
#' @export
#'
update_pi.susiF_ss <- function( susiF_ss.obj, l, tpi, ...)
{

  if( l > length(susiF_ss.obj$est_pi))
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  if( class(tpi)%!in% c("pi_mixture_normal" , "pi_mixture_normal_per_scale"))
  {
    stop("Error tpi should be of one of the follwoing class:\n
          pi_mixture_normal \n pi_mixture_normal_per_scale")
  }
  susiF_ss.obj$est_pi[[l]] <- tpi
  out <-susiF_ss.obj
  class(out) <- "susiF_ss"
  return(out)
}




#' @rdname update_alpha
#'
#' @method update_alpha susiF_ss
#'
#' @export update_alpha.susiF_ss
#'
#' @export
#'
update_alpha.susiF_ss <-  function(susiF_ss.obj, l, alpha, ... )
{
  susiF_ss.obj$alpha[[l]] <- alpha
  susiF_ss.obj$alpha_hist[[ (length(susiF_ss.obj$alpha_hist)+1)  ]] <- alpha
  return( susiF_ss.obj)
}



#' @rdname get_alpha
#'
#' @method get_alpha susiF_ss
#'
#' @export get_alpha.susiF_ss
#'
#' @export
#'

get_alpha.susiF_ss <-  function(susiF_ss.obj, l, ...  )
{
  out <- susiF_ss.obj$alpha[[l]]
  return( out)
}


#' @rdname update_lBF
#'
#' @method update_lBF susiF_ss
#'
#' @export update_lBF.susiF_ss
#'
#' @export
#'

update_lBF.susiF_ss <- function(susiF_ss.obj,l, lBF, ...)
{
  if(l> susiF_ss.obj$L)
  {
    stop("Error: trying to update more effects than the number of specified effect")
  }

  susiF_ss.obj$lBF[[l]] <- lBF
  return(susiF_ss.obj)
}



#' @rdname update_ELBO
#'
#' @method update_ELBO susiF_ss
#'
#' @export update_ELBO.susiF_ss
#'
#' @export
#'

update_ELBO.susiF_ss <- function (susiF_ss.obj,ELBO, ...)
{

  susiF_ss.obj$ELBO <- c(susiF_ss.obj$ELBO,ELBO)
  return(susiF_ss.obj)
}





#' @rdname update_cal_pip
#'
#' @method update_cal_pip susiF_ss
#'
#' @export update_cal_pip.susiF_ss
#'
#' @export
#'

update_cal_pip.susiF_ss <- function (susiF_ss.obj, ...)
{
  if(sum( is.na(unlist(susiF_ss.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  tpip <- list()
  for ( l in 1:susiF_ss.obj$L)
  {
    tpip[[l]] <- rep(1, lengths(susiF_ss.obj$alpha)[[l]])-susiF_ss.obj$alpha[[l]]
  }
  susiF_ss.obj$pip <- 1-  apply( do.call(rbind,tpip),2, prod)
  return(susiF_ss.obj)
}


