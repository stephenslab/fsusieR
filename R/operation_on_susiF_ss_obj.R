################################## Operations on susiF_ss object ############################
#'
#'






#' @title Compute expected residuals for susiF_ss.obj
#'
#' @description Compute expected residuals for susiF_ss.obj
#'
#' @param  susiF_ss.obj  susiF_ss.obj object
#'
#' @param data  an object of the class suff_stat define by function \code{\link{make_data_suff_stat}}
#'
#' @return an update suff_stat object
#'
#' @export
#'
cal_expected_residual <- function (susiF_ss.obj, data, ...)
  UseMethod("cal_expected_residual")

#' @rdname cal_expected_residual
#'
#' @method cal_expected_residual susiF_ss
#'
#' @export cal_expected_residual.susiF_ss
#'
#' @export
#'
cal_expected_residual.susiF_ss <- function( susiF_ss.obj, data)
{

  data$exp_residual <- data$Xty - data$XtX%*%get_post_F(susiF_ss.obj)

  return( data)
}





#' @rdname cal_Bhat_Shat
#'
#' @method cal_Bhat_Shat  susiF_ss
#'
#' @export cal_Bhat_Shat.susiF_ss
#'
#' @export
#'
cal_Bhat_Shat.susiF_ss <- function( susiF_ss.obj,data , partial=TRUE )
{


  if (partial)
  {
    rho <-data$part_exp_residual
  }else{
    rho <-data$Xty
  }
  Bhat  <- rho /diag(XtX)

  Shat  <- sqrt( susiF_ss.obj$sigma2 * matrix(1, ncol = ncol(Bhat), nrow= nrow(Bhat))/ diag(data$XtX))

  out  <- list( Bhat = Bhat,
                Shat = Shat)
  return(out)
}


#' @rdname estimate_residual_variance
#'
#' @method estimate_residual_variance susiF_ss
#'
#' @export estimate_residual_variance.susiF_ss
#'
#' @export
estimate_residual_variance.susiF_ss <- function(susiF_ss.obj,data, ... )
{
  out <-  (1/(data$N*ncol(data$Bhat)))*get_ER2 (susiF_ss.obj,data  )
  return(out)
}


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
    fitted_wc[[l]]        <-  matrix(0, nrow = nrow(data$Bhat), ncol=ncol(data$Bhat)  )
    fitted_wc2[[l]]       <-  matrix(1, nrow = nrow(data$Bhat), ncol=ncol(data$Bhat)  )
    alpha [[l]]           <-  rep(1, P )/P
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
#' @method get_pi susiF_ss
#'
#' @export get_pi.susiF_ss
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
#' @method get_lBF susiF_ss
#'
#' @export get_lBF.susiF_ss
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





#' @title Compute partial residuals for susiF_ss.obj
#'
#' @description Compute partial residuals for susiF_ss.obj
#'
#' @param  susiF_ss.obj  susiF_ss.obj object
#'
#' @param data  an object of the class suff_stat define by function \code{\link{make_data_suff_stat}}
#'
#' @param l, optional effect to update
#'
#' @return an update suff_stat object
#'
#' @export
#'
get_partial_residual <- function (susiF_ss.obj, data, l, ...)
  UseMethod("get_partial_residual")

#' @rdname get_partial_residual
#'
#' @method get_partial_residual susiF_ss
#'
#' @export get_partial_residual.susiF_ss
#'
#' @export
#'
get_partial_residual.susiF_ss <- function( susiF_ss.obj, data,l)
{

  data$part_exp_residual <- data$exp_residual + data$XtX%*%get_post_F(susiF_ss.obj,l)

  return( data)
}






#' @rdname get_post_F
#'
#' @method get_post_F susiF_ss
#'
#' @export get_post_F.susiF_ss
#'
#' @export
#'

get_post_F.susiF_ss <- function(susiF_ss.obj,l,...)
{
  if(missing(l))
  {
    out <-  Reduce("+",lapply(1:susiF_ss.obj$L, FUN=function(l) susiF_ss.obj$alpha[[l]] * susiF_ss.obj$fitted_wc[[l]]))
  }else{
    out <-   susiF_ss.obj$alpha[[l]] * susiF_ss.obj$fitted_wc[[l]]
  }

  return(out)
}



#' @rdname get_post_F2
#'
#' @method get_post_F2 susiF_ss
#'
#' @export get_post_F2.susiF_ss
#'
#' @export
#'
get_post_F2.susiF_ss <- function(susiF_ss.obj, l,...)
{
  if(missing(l))
  {
    out <-  Reduce("+",lapply(1:susiF_ss.obj$L, FUN=function(l) susiF_ss.obj$alpha[[l]] *(susiF_ss.obj$sigma2 + susiF_ss.obj$fitted_wc2[[l]])))
  }else{
    out <-   susiF_ss.obj$alpha[[l]] *(susiF_ss.obj$sigma2 + susiF_ss.obj$fitted_wc2[[l]])
  }

  return(out)
}



#' @rdname get_ER2
#'
#' @method get_ER2 susiF_ss
#'
#' @export get_ER2.susiF_ss
#'
#' @export

get_ER2.susiF_ss = function (  susiF_ss.obj,data ) {
  postF <- get_post_F( susiF_ss.obj )# J by N matrix
  postF2 <- get_post_F2( susiF_ss.obj ) # Posterior second moment.
  return( sum(data$yty) -2*sum(t(postF)%*%Xty) +  sum(( t(postF) %*% data$XtX %*% postF))  -sum(postF)^2 + sum( postF2))

}


#' @rdname out_prep
#'
#' @method out_prep susiF_ss
#'
#' @export out_prep.susiF_ss
#'
#' @export
#'

out_prep.susiF_ss <- function(susiF_ss.obj, ...)
{
  susiF_ss.obj <-  update_cal_pip(susiF_ss.obj)
  susiF_ss.obj <-  update_cal_cs(susiF_ss.obj)
  return(susiF_ss.obj)
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


#' @rdname update_cal_cs
#'
#' @method update_cal_cs susiF_ss
#'
#' @export update_cal_cs.susiF_ss
#'
#' @export
#'

update_cal_cs.susiF_ss <- function(susiF_ss.obj, cov_lev=0.95)
{
  if(sum( is.na(unlist(susiF_ss.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  for ( l in 1:susiF_ss.obj$L)
  {
    temp        <- susiF_ss.obj$alpha[[l]]
    temp_cumsum <- cumsum( temp[order(temp, decreasing =TRUE)])
    max_indx_cs <- min(which( temp_cumsum >cov_lev ))
    susiF_ss.obj$cs[[l]]  <- order(temp, decreasing = TRUE)[1:max_indx_cs ]

  }

  return(susiF_ss.obj)
}


#' @rdname update_susiF_obj
#'
#' @method update_susiF_obj  susiF_ss
#'
#' @export update_susiF_obj.susiF_ss
#'
#' @export
#'

update_susiF_obj.susiF_ss <- function(susiF_ss.obj, l, EM_pi,Bhat, Shat , indx_lst, ...)
{

  if( l > length(susiF_ss.obj$est_pi))
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  if(  "EM_pi"  %!in%  class(EM_pi)  )
  {
    stop("Error EM_pi should be of the class EM_pi")
  }
  susiF_ss.obj         <-   update_pi(susiF_ss.obj   = susiF_ss.obj   ,
                                      l              = l ,
                                      tpi            =  EM_pi$tpi_k
  )

  susiF_ss.obj$G_prior          <-  update_prior(get_G_prior(susiF_ss.obj  ) , EM_pi$tpi_k  )
  susiF_ss.obj$fitted_wc[[l]]   <-  post_mat_mean(get_G_prior(susiF_ss.obj) , Bhat=Bhat, Shat=Shat,indx_lst= indx_lst )
  susiF_ss.obj$fitted_wc2[[l]]  <-  post_mat_sd  (get_G_prior(susiF_ss.obj) , Bhat=Bhat, Shat=Shat, indx_lst= indx_lst)^2


  new_alpha    <- cal_zeta(  EM_pi$lBF)
  susiF_ss.obj <- update_alpha(susiF_ss.obj, l, new_alpha)
  susiF_ss.obj <- update_lBF(susiF_ss.obj, l, EM_pi$lBF)
  return(susiF_ss.obj)
}





#' @rdname update_residual_variance
#'
#' @method update_residual_variance susiF_ss
#'
#' @export update_residual_variance.susiF_ss
#'
#' @export
#'

update_residual_variance.susiF_ss <- function(susiF_ss.obj,sigma2)
{
  susiF_ss.obj$sigma2 <- sigma2
  return(susiF_ss.obj)
}



#' @title Update expected residuals for susiF_ss.obj
#'
#' @description Update expected residuals for susiF_ss.obj for effect l
#'
#' @param  susiF_ss.obj  susiF_ss.obj object
#'
#' @param data  an object of the class suff_stat define by function \code{\link{make_data_suff_stat}}
#'
#' @param l, optional effect to update
#'
#' @return an updated suff_stat object
#'
#' @export
#'
update_expected_residual  <- function( susiF_ss.obj, data,l)
  UseMethod("update_expected_residual")



#' @rdname update_expected_residual
#'
#' @method update_expected_residual susiF_ss
#'
#' @export update_expected_residual.susiF_ss
#'
#' @export
#'
update_expected_residual.susiF_ss <- function( susiF_ss.obj, data,l)
{

  data$exp_residual <- data$part_exp_residual - data$XtX%*%get_post_F(susiF_ss.obj,l)

  return( data)
}






#' @rdname update_KL
#'
#' @method update_KL susiF_ss
#'
#' @export update_KL.susiF_ss
#'
#' @export
#'
update_KL.susiF_ss <- function(susiF_ss.obj, data , ...)
{
  susiF_ss.obj$KL <-  do.call(c,lapply(1:susiF_ss.obj$L,FUN=function(l) cal_KL_l(susiF_ss.obj, l, data  )))
  return( susiF_ss.obj)
}




