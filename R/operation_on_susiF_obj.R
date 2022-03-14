################################## Operations on susiF object ############################
#'
#' @title Initialize a susiF object using regression coefficients
#'
#' @param L number of non zero coefficients An L-vector containing the indices of the
#'   nonzero coefficients.
#'
#' @param G_prior prior object defined by init_prior function
#'
#' @param Y Matrix of outcomes
#'
#' @param X Matrix of covariates
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
init_susiF_obj <- function(L, G_prior, Y,X )
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
  N               <- dim(Y)[1]
  n_wac           <- dim(Y)[2]
  P               <- dim(X)[2]
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

  class(obj) <- "susiF"
  return(obj)
}

#'
#' @title Access susiF mixture proportion of effect l
#'
#' @param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'
#' @return a vector of of proportion
#'
#' @export
#'
get_pi  <- function(susiF.obj, l, ...)
  UseMethod("get_pi")

#' @rdname get_pi
#'
#' @method get_pi susiF
#'
#' @export get_pi.susiF
#'
#' @export
#'
get_pi.susiF <- function(susiF.obj, l, ...)
{

  if( l >  length(susiF.obj$est_pi))
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  out <- susiF.obj$est_pi[[l]]
  return(out)
}





#' @title Access susiF log Bayes factors of effect l
#'
#' @param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'
#' @return a vector of log Bayes Factors
#'
#' @export
#'
get_lBF  <- function(susiF.obj, l, ...)
  UseMethod("get_lBF")

#' @rdname get_pi
#'
#' @method get_pi susiF
#'
#' @export get_pi.susiF
#'
#' @export
#'
get_lBF.susiF <- function(susiF.obj, l, ...)
{

  if( l >   susiF.obj$L)
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  out <- susiF.obj$lBF[[l]]
  return(out)
}





#'
#'@title Access susiF internal prior
#'
#' @param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'
#' @return G_prior object
#'
#' @export
#'
#'
get_G_prior  <- function(susiF.obj, ...)
  UseMethod("get_G_prior")


#' @rdname get_G_prior
#'
#' @method get_G_prior susiF
#'
#' @export get_G_prior.susiF
#'
#' @export
#'
get_G_prior.susiF <- function(susiF.obj, ...)
{
  out <- susiF.obj$G_prior
  return(out)
}

#' @title Update mixture proportion of susiF mixture proportions of effect l
#'
#' @param susiF.obj a susif object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'
#' @param tpi an object of the class "pi_mixture_normal" or "pi_mixture_normal_per_scale"
#'
#' @return susiF object
#'
#' @export

update_pi <- function( susiF.obj, l, tpi, ...)
    UseMethod("update_pi")

#' @rdname update_pi
#'
#' @method update_pi susiF
#'
#' @export update_pi.susiF
#'
#' @export
#'
update_pi.susiF <- function( susiF.obj, l, tpi, ...)
{

  if( l > length(susiF.obj$est_pi))
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
  susiF.obj$est_pi[[l]] <- tpi
  out <- susiF.obj
  class(out) <- "susiF"
  return(out)
}

#' @title Update alpha   susiF mixture proportion of effect l
#'
#' @param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'
#' @param alpha  vector of p alpha values summing up to one
#'
#' @return susiF object
#'
#' @export
#'

update_alpha  <-  function(susiF.obj, l, alpha, ... )
        UseMethod("update_alpha")


#' @rdname update_alpha
#'
#' @method update_alpha susiF
#'
#' @export update_alpha.susiF
#'
#' @export
#'
update_alpha.susiF <-  function(susiF.obj, l, alpha, ... )
{
  susiF.obj$alpha[[l]] <- alpha
  susiF.obj$alpha_hist[[ (length(susiF.obj$alpha_hist)+1)  ]] <- alpha
  return( susiF.obj)
}

#' @title Update alpha  susiF mixture proportion of effect l
#'
#' @param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'
#' @param alpha  vector of p alpha values summing up to one
#'
#' @return susiF object
#'
#' @export
#'
#'
get_alpha  <-  function(susiF.obj, l, ...  )
  UseMethod("get_alpha")

#' @rdname get_alpha
#'
#' @method get_alpha susiF
#'
#' @export get_alpha.susiF
#'
#' @export
#'

get_alpha.susiF <-  function(susiF.obj, l, ...  )
{
  out <- susiF.obj$alpha[[l]]
  return( out)
}



#' @title Compute partial residual for effect l
#'
#' @param susiF.obj a susisF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'
#' @param X matrix of covariate
#'
#' @param D matrix of wavelet D coefficients from the original input data (Y)
#'
#' @param C vector of wavelet scaling coefficient from the original input data (Y)
#'
#' @param indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'
#' @return a matrix of size N by size J of partial residuals
#'
#' @export
cal_partial_resid  <- function( susiF.obj, l, X, D, C,  indx_lst, ... )
      UseMethod("cal_partial_resid")


#' @rdname cal_partial_resid
#'
#' @method cal_partial_resid susiF
#'
#' @export cal_partial_resid.susiF
#'
#' @export
#'

cal_partial_resid.susiF  <- function( susiF.obj, l, X, D, C,  indx_lst, ... )
{
   L <- susiF.obj$L
  if (L > 1){
    id_L <- (1:L)[ - ( (l%%L)+1) ]#Computing residuals R_{l+1} by removing all the effect except effect l+1

    if(class(get_G_prior(susiF.obj))=="mixture_normal_per_scale" )
    {
      update_D  <-  D - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% (susiF.obj$fitted_wc[[l]][,-indx_lst[[length(indx_lst)]]])   ) )
      update_C  <-  C - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% susiF.obj$fitted_wc[[l]][,indx_lst[[length(indx_lst)]]] ) )
      update_Y  <- cbind(  update_D, update_C)
    }
    if(class(get_G_prior(susiF.obj))=="mixture_normal" )
    {
      id_L <- (1:L)[ - ( (l%%L)+1) ]#Computing residuals R_{l+1} by removing all the effect except effect l+1
      update_D  <-  D - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% (susiF.obj$fitted_wc[[l]][,-dim(susiF.obj$fitted_wc[[l]])[2]])   ) )
      update_C  <-  C - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% susiF.obj$fitted_wc[[l]][,dim(susiF.obj$fitted_wc[[l]])[2]] ) )
      update_Y  <- cbind(  update_D, update_C)
    }
  }else{
    id_L <- 1

    if(class(get_G_prior(susiF.obj))=="mixture_normal_per_scale" )
    {
      update_D  <-  D - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% (susiF.obj$fitted_wc[[l]][,-indx_lst[[length(indx_lst)]]])   ) )
      update_C  <-  C - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% susiF.obj$fitted_wc[[l]][,indx_lst[[length(indx_lst)]]] ) )
      update_Y  <- cbind(  update_D, update_C)
    }
    if(class(get_G_prior(susiF.obj))=="mixture_normal" )
    {
      update_D  <-  D - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% (susiF.obj$fitted_wc[[l]][,-dim(susiF.obj$fitted_wc[[l]])[2]])   ) )
      update_C  <-  C - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% susiF.obj$fitted_wc[[l]][,dim(susiF.obj$fitted_wc[[l]])[2]] ) )
      update_Y  <- cbind(  update_D, update_C)
    }
  }


  return(update_Y)
}


#'@title Update  susiF object using the output of EM_pi
#'
#' @param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'
#' @param EM_pi an object of the class "EM_pi" generated by the function EM_pi
#'
#' @return susiF object
#'
#' @export

update_susiF_obj  <- function(susiF.obj, l, EM_pi, Bhat, Shat, indx_lst, ...)
      UseMethod("update_susiF_obj")

#' @rdname update_susiF_obj
#'
#' @method update_susiF_obj susiF
#'
#' @export update_susiF_obj.susiF
#'
#' @export
#'

update_susiF_obj.susiF <- function(susiF.obj, l, EM_pi, Bhat, Shat, indx_lst, ...)
{

  if( l > length(susiF.obj$est_pi))
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
  susiF.obj         <-   update_pi(susiF.obj = susiF.obj ,
                                   l = l ,
                                   tpi =  EM_pi$tpi_k)
  susiF.obj$G_prior <-   update_prior(get_G_prior(susiF.obj) , EM_pi$tpi_k  )

  susiF.obj$fitted_wc[[l]]   <- post_mat_mean(get_G_prior(susiF.obj) , Bhat, Shat,indx_lst= indx_lst )
  susiF.obj$fitted_wc2[[l]]  <- post_mat_sd  (get_G_prior(susiF.obj) , Bhat, Shat, indx_lst= indx_lst)^2


  new_alpha <- cal_zeta(  EM_pi$lBF)
  susiF.obj <- update_alpha(susiF.obj, l, new_alpha)
  susiF.obj <- update_lBF(susiF.obj, l, EM_pi$lBF)
  return(susiF.obj)
}


#'@title Update susiF log Bayes factor
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@param l effect to update
#'@param lBF vector of length p, containning the updated log Bayes factors
#'@return susiF object
#'@export

update_lBF  <- function    (susiF.obj, l, lBF, ...)
  UseMethod("update_lBF")

#' @rdname update_lBF
#'
#' @method update_lBF susiF
#'
#' @export update_lBF.susiF
#'
#' @export
#'

update_lBF.susiF <- function    (susiF.obj,l, lBF, ...)
{
  if(l> susiF.obj$L)
  {
    stop("Error: trying to update more effects than the number of specified effect")
  }

  susiF.obj$lBF[[l]] <- lBF
  return(susiF.obj)
}





#'@title Update susiF log Bayes factor
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@param  ELBO new ELBO value
#'@param lBF vector of length p, containning the updated log Bayes factors
#'@return susiF object
#'@export

update_ELBO  <- function    (susiF.obj,ELBO , ...)
  UseMethod("update_ELBO")

#' @rdname update_ELBO
#'
#' @method update_ELBO susiF
#'
#' @export update_ELBO.susiF
#'
#' @export
#'

update_ELBO.susiF <- function    (susiF.obj,ELBO, ...)
{

  susiF.obj$ELBO <- c(susiF.obj$ELBO,ELBO)
  return(susiF.obj)
}

#'@title Update susiF by computing PiP
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@return susiF object
#'@export

update_cal_pip  <- function (susiF.obj, ...)
  UseMethod("update_cal_pip")

#' @rdname update_cal_pip
#'
#' @method update_cal_pip susiF
#'
#' @export update_cal_pip.susiF
#'
#' @export
#'

update_cal_pip.susiF <- function (susiF.obj, ...)
{
  if(sum( is.na(unlist(susiF.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  tpip <- list()
  for ( l in 1:susiF.obj$L)
  {
    tpip[[l]] <- rep(1, lengths(susiF.obj$alpha)[[l]])-susiF.obj$alpha[[l]]
  }
  susiF.obj$pip <- 1-  apply( do.call(rbind,tpip),2, prod)
  return(susiF.obj)
}


#'@title Update susiF by computing credible sets
#'
#' @param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'
#' @param cov_lev numeric between 0 and 1, corresponding to the expected level of coverage of the cs if not specified set to 0.95
#'
#' @return susiF object
#'
#' @export

update_cal_cs  <- function(susiF.obj, cov_lev=0.95, ...)
      UseMethod("update_cal_cs")

#' @rdname update_cal_cs
#'
#' @method update_cal_cs susiF
#'
#' @export update_cal_cs.susiF
#'
#' @export
#'

update_cal_cs.susiF <- function(susiF.obj, cov_lev=0.95)
{
  if(sum( is.na(unlist(susiF.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  for ( l in 1:susiF.obj$L)
  {
    temp        <- susiF.obj$alpha[[l]]
    temp_cumsum <- cumsum( temp[order(temp, decreasing =TRUE)])
    max_indx_cs <- min(which( temp_cumsum >cov_lev ))
    susiF.obj$cs[[l]]  <- order(temp, decreasing = TRUE)[1:max_indx_cs ]

  }

  return(susiF.obj)
}

#'@title Update susiF by computing predicted curves
#'
#'@param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'@param Y functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelet which assume that J is of the form J^2. If J not a power of 2, susiF internally remaps the data into grid of length 2^J
#'@param X matrix of size N by p
#'#'@param indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'@return susiF object
#'@export
#'

update_cal_indf <- function(susiF.obj, Y, X, indx_lst, ...)
      UseMethod("update_cal_indf")

#' @rdname update_cal_indf
#'
#' @method update_cal_indf susiF
#'
#' @export update_cal_indf.susiF
#'
#' @importFrom wavethresh wr
#'
#' @importFrom wavethresh wd
#'
#' @export
#'

update_cal_indf.susiF <- function(susiF.obj, Y, X, indx_lst, ...)
{
  if(sum( is.na(unlist(susiF.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  temp <- wd(rep(0, susiF.obj$n_wac)) #create dummy wd object


  if(class(get_G_prior(susiF.obj))=="mixture_normal_per_scale" )
  {
    for ( i in 1:susiF.obj$N)
    {
      susiF.obj$ind_fitted_func[i,]  <- rep(0,susiF.obj$n_wac)#fitted_baseline future implementation
      for ( l in 1:susiF.obj$L)
      {
        #add wavelet coefficient
        temp$D                         <-    (susiF.obj$alpha[[l]] *X[i,])%*%susiF.obj$fitted_wc[[l]][,-indx_lst[[length(indx_lst)]]]
        temp$C[length(temp$C)]         <-    (susiF.obj$alpha[[l]] *X[i,])%*%susiF.obj$fitted_wc[[l]][,indx_lst[[length(indx_lst)]]]
        #transform back
        susiF.obj$ind_fitted_func[i,]  <-  susiF.obj$ind_fitted_func[i,]+wr(temp)
      }
    }
  }
  if(class(get_G_prior(susiF.obj))=="mixture_normal" )
  {
    for ( i in 1:susiF.obj$N)
    {
      susiF.obj$ind_fitted_func[i,]  <- rep(0,susiF.obj$n_wac)#fitted_baseline
      for ( l in 1:susiF.obj$L)
      {
        #add wavelet coefficient
        temp$D                         <-    (susiF.obj$alpha[[l]] *X[i,])%*%susiF.obj$fitted_wc[[l]][,-dim(susiF.obj$fitted_wc[[l]])[2]]
        temp$C[length(temp$C)]         <-    (susiF.obj$alpha[[l]] *X[i,]) %*%susiF.obj$fitted_wc[[l]][,dim(susiF.obj$fitted_wc[[l]])[2]]
        #transform back
        susiF.obj$ind_fitted_func[i,]  <-  susiF.obj$ind_fitted_func[i,]+wr(temp)
      }
    }
  }
  return( susiF.obj)
}



#' @title Update susiF by computing posterior curves
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param Y functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelet which assume that J is of the form J^2. If J not a power of 2, susiF internally remaps the data into grid of length 2^J
#'
#' @param X matrix of size N by p
#'
#' @param indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'
#' @return susiF object
#'
#' @export
update_cal_fit_func  <- function(susiF.obj, indx_lst, ...)
    UseMethod("update_cal_fit_func")

#' @rdname update_cal_fit_func
#'
#' @method update_cal_fit_func susiF
#'
#' @export update_cal_fit_func.susiF
#'
#' @importFrom wavethresh wr
#'
#' @importFrom wavethresh wd
#'
#' @export
#'

update_cal_fit_func.susiF <- function(susiF.obj, indx_lst, ...)
{

  if(sum( is.na(unlist(susiF.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  temp <- wd(rep(0, susiF.obj$n_wac))

  if(class(get_G_prior(susiF.obj))=="mixture_normal_per_scale" )
  {
      for ( l in 1:susiF.obj$L)
    {
      temp$D                     <- (susiF.obj$alpha[[l]])%*%susiF.obj$fitted_wc[[l]][,-indx_lst[[length(indx_lst)]]]
      temp$C[length(temp$C)]     <- (susiF.obj$alpha[[l]])%*%susiF.obj$fitted_wc[[l]][,indx_lst[[length(indx_lst)]]]
      susiF.obj$fitted_func[[l]] <- wr(temp)
    }
  }
  if(class(get_G_prior(susiF.obj))=="mixture_normal" )
  {
     for ( l in 1:susiF.obj$L)
    {
      temp$D                     <- (susiF.obj$alpha[[l]])%*%susiF.obj$fitted_wc[[l]][,-dim(susiF.obj$fitted_wc[[l]])[2]]
      temp$C[length(temp$C)]     <- (susiF.obj$alpha[[l]])%*%susiF.obj$fitted_wc[[l]][,dim(susiF.obj$fitted_wc[[l]])[2]]
      susiF.obj$fitted_func[[l]] <- wr(temp)
    }
  }
  return(susiF.obj)
}

#' @title Update residual variance
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param sigma2 estimate residual variance
#'
#' @return updated susiF.obj




update_residual_variance  <- function(susiF.obj,sigma2, ...)
  UseMethod("update_residual_variance")

#' @rdname update_residual_variance
#'
#' @method update_residual_variance susiF
#'
#' @export update_residual_variance.susiF
#'
#' @export
#'

update_residual_variance.susiF <- function(susiF.obj,sigma2)
{
  susiF.obj$sigma2 <- sigma2
  return(susiF.obj)
}



#' @title Preparing output of main susiF function
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param Y functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelets that assume that J is of the form J^2. If J is not a power of 2, susiF internally remaps the data into a grid of length 2^J
#'
#' @param X matrix of size N by p
#'
#' @param indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'
#' @return susiF object
#'
#' @export
#'
out_prep <- function(susiF.obj,Y, X, indx_lst, ...)
      UseMethod("out_prep")

#' @rdname out_prep
#'
#' @method out_prep susiF
#'
#' @export out_prep.susiF
#'
#' @export
#'

out_prep.susiF <- function(susiF.obj,Y, X, indx_lst, ...)
{
  susiF.obj <-  update_cal_pip(susiF.obj)
  susiF.obj <-  update_cal_cs(susiF.obj)
  susiF.obj <-  update_cal_indf(susiF.obj, Y, X, indx_lst)
  susiF.obj <-  update_cal_fit_func(susiF.obj, indx_lst)
  return(susiF.obj)
}


#' @title Compute posterior mean of the fitted effect
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l , optional effect to update
#'
#' @return  A J by T matrix of posterior wavelet coefficient,
#' \item{if l not missng}{return effect specific posterior mean }
#' \item{if l not missng}{return effect specific posterior mean}
get_post_F <- function(susiF.obj,l,...)
  UseMethod("get_post_F")

#' @rdname get_post_F
#'
#' @method get_post_F susiF
#'
#' @export get_post_F.susiF
#'
#' @export
#'

get_post_F.default <- function(susiF.obj,l,...)
{
  if(missing(l))
  {
    out <-  Reduce("+",lapply(1:susiF.obj$L, FUN=function(l) susiF.obj$alpha[[l]] * susiF.obj$fitted_wc[[l]]))
  }else{
    out <-   susiF.obj$alpha[[l]] * susiF.obj$fitted_wc[[l]]
  }

   return(out)
}



#' @title Compute posterior second moment
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#' @param l , optional effect to update
#'
#' @return  A J by T matrix of posterior wavelet coefficient,
#' \item{if l not missng}{return effect specific posterior second moment }
#' \item{if l not missng}{return effect specific posterior second moment}

get_post_F2 <- function(susiF.obj,l,...)
  UseMethod("get_post_F2")

#' @rdname get_post_F2
#'
#' @method get_post_F2 susiF
#'
#' @export get_post_F2.susiF
#'
#' @export
#'
get_post_F2.default <- function(susiF.obj, l,...)
{
  if(missing(l))
  {
    out <-  Reduce("+",lapply(1:susiF.obj$L, FUN=function(l) susiF.obj$alpha[[l]] *(susiF.obj$s2+ susiF.obj$fitted_wc2[[l]])))
  }else{
    out <-   susiF.obj$alpha[[l]] *(susiF.obj$s2+ susiF.obj$fitted_wc2[[l]])
  }

  return(out)
}


#' @title Update residual variance
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param Y wavelet transformed  functional phenotype, matrix of size N by size J.
#'
#' @param X matrix of size N by p
#'
#'
#' @return estimated residual variance
#' @export
estimate_residual_variance <- function(susiF.obj,Y,X, ... )
  UseMethod("estimate_residual_variance")


#' @rdname estimate_residual_variance
#'
#' @method estimate_residual_variance susiF
#'
#' @export estimate_residual_variance.susiF
#'
#' @export
estimate_residual_variance.susiF <- function(susiF.obj,Y,X, ... )
{
  out <-  (1/(prod(dim(Y))))*get_ER2 (susiF.obj,Y, X  )
  return(out)
}




#' @titleCompute Epected sum of square
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param Y wavelet transformed  functional phenotype, matrix of size N by size J.
#'
#' @param X matrix of size N by p
#'
#'
#' @return estimated residual variance
#' @export
get_ER2 <- function(susiF.obj,Y,X, ... )
  UseMethod("estimate_residual_variance")


#' @rdname get_ER2
#'
#' @method get_ER2 susiF
#'
#' @export get_ER2.susiF
#'
#' @export

get_ER2.susiF = function (  susiF.obj,Y, X) {
  postF <- get_post_F(susiF.obj )# J by N matrix
  Xr_L = t(X%*% postF)
  postF2 <- get_post_F2(susiF.obj ) # Posterior second moment.
  return(sum((Y - X%*%postF )^2)  -sum(postF)^2 + sum(postF2))
}
