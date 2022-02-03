################################## Operations on susiF object ############################
#'
#'@title Initialize a susiF object using regression coefficients
#'
#'@param L number of non zero coefficients An L-vector containing the indices of the
#'   nonzero coefficients.
#'
#' @param G_prior prior object defined by init_prior function
#' @param Y Matrix of outcomes
#' @param X Matrix of covariate
#' @return A list with elements
#' fitted_wc2
#' fitted_wc2
#' alpha_hist
#' ind_fitted_func
#' cs (credible set)
#' pip Posterior inclusion probabilites
#' G_prior a G_prior of the same class as the input G_prior, used for internal calculation
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
  for ( l in 1:L )
  {
    fitted_wc[[l]]        <-  matrix(NA, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    fitted_wc2[[l]]       <-  matrix(NA, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    alpha [[l]]           <-  rep(NA, dim(X)[2])
    cs[[l]]               <-  list()
    est_pi [[l]]          <-  get_pi_G_prior(G_prior)
    est_sd [[l]]          <-  get_sd_G_prior(G_prior)
  }
  obj <- list( fitted_wc       = fitted_wc,
               fitted_wc2      = fitted_wc2,
               ind_fitted_func = ind_fitted_func,
               G_prior         = G_prior,
               alpha_hist      = alpha_hist,
               N               = N,
               n_wac           = n_wac,
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
#'@title Access susiF mixture proportion of effect l
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'@return susiF object
#'@export
get_pi.susiF <- function(susiF.obj, l)
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


#'
#'@title Access susiF internal prior
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@return G_prior object
#'@export
get_G_prior.susiF <- function(susiF.obj)
{
  out <- susiF.obj$G_prior
  return(out)
}

#'@title Update mixture proportion of susiF mixture proportions of effect l
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'@param tpi an object of the class "pi_mixture_normal" or "pi_mixture_normal_per_scale"
#'@return susiF object
#'@export
update_pi.susiF <- function( susiF.obj, l, tpi)
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

#'@title Update alpha   susiF mixture proportion of effect l
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'@param alpha  vector of p alpha values summing up to one
#'@return susiF object
#'@export
update_alpha.susiF <-  function(susiF.obj, l, alpha )
{
  susiF.obj$alpha[[l]] <- alpha
  susiF.obj$alpha_hist[[ (length(susiF.obj$alpha_hist)+1)  ]] <- alpha
  return( susiF.obj)
}

#'@title Update alpha   susiF mixture proportion of effect l
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'@param alpha  vector of p alpha values summing up to one
#'@return susiF object
#'@export
get_alpha.susiF <-  function(susiF.obj, l  )
{
  out <- susiF.obj$alpha[[l]]
  return( out)
}



#'@title Compute partial residual for effect l
#'
#'@param susiF.obj a susisF object defined by \code{\link{init_susiF_obj}} function
#'@param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'@param X matrix of covariate
#'@param D matrix of wavelet D coefficients from the original input data (Y)
#'@param C vector of wavelet scaling coefficient from the original input data (Y)
#'@param L total number of effect
#'@param indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'@return a matrix of size N by size J of partial residuals
#'@export



cal_partial_resid.susiF  <- function( susiF.obj, l, X, D, C, L, indx_lst )
{
  if (L > length(susiF.obj$alpha)){
    stop( "Error: trying to access more effect that initialized number of effect in susiF.obj")
  }
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
  }
  else{
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
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'@param EM_pi an object of the class "EM_pi" generated by the function EM_pi
#'@return susiF object
#'@export
update_susiF_obj.susiF <- function(susiF.obj, l, EM_pi, Bhat, Shat, indx_lst)
{

  if( l > length(susiF.obj$est_pi))
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  if( class(EM_pi)%!in% c("EM_pi"))
  {
    stop("Error EM_pi should be of the class EM_pi")
  }
  susiF.obj         <-   update_pi(susiF.obj = susiF.obj ,
                                   l = l ,
                                   tpi =  EM_pi$tpi_k)
  susiF.obj$G_prior <-   update_prior(get_G_prior(susiF.obj) , EM_pi$tpi_k  )

  susiF.obj$fitted_wc[[l]]   <- post_mat_mean(get_G_prior(susiF.obj) , Bhat, Shat,indx_lst= indx_lst )
  susiF.obj$fitted_wc2[[l]]  <- post_mat_sd  (get_G_prior(susiF.obj) , Bhat, Shat, indx_lst= indx_lst)


  new_alpha <- cal_zeta(  EM_pi$lBF)
  susiF.obj <- update_alpha.susiF(susiF.obj, l, new_alpha)
  return(susiF.obj)
}


#'@title Update susiF by computing PiP
#'
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@return susiF object
#'@export
update_cal_pip.susiF <- function    (susiF.obj)
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
#'@param susiF.obj a susisf object defined by \code{\link{init_susiF_obj}} function
#'@param cov_lev numeric between 0 and 1, corresponding to the expected level of coverage of the cs if not specified set to 0.95
#'@return susiF object
#'@export
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
#'@param X matrix of size n by p
#'#'@param indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'@return susiF object
#'@export
update_cal_indf.susiF <- function(susiF.obj, Y, X, indx_lst)
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



#'@title Update susiF by computing posterior curves
#'
#'@param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'@param Y functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelet which assume that J is of the form J^2. If J not a power of 2, susiF internally remaps the data into grid of length 2^J
#'@param X matrix of size n by p
#'#'@param indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'@return susiF object
#'@export
update_cal_fit_func.susiF <- function(susiF.obj, indx_lst)
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




#'@title Preparing output of main susiF function
#'
#'@param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'@param Y functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelets that assume that J is of the form J^2. If J is not a power of 2, susiF internally remaps the data into a grid of length 2^J
#'@param X matrix of size n by p
#'#'@param indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'@return susiF object
#'@export
#'
out_prep.susiF <- function(susiF.obj,Y, X, indx_lst)
{
  susiF.obj <-  update_cal_pip(susiF.obj)
  susiF.obj <-  update_cal_cs(susiF.obj)
  susiF.obj <-  update_cal_indf(susiF.obj, Y, X, indx_lst)
  susiF.obj <-  update_cal_fit_func(susiF.obj, indx_lst)
  return(susiF.obj)
}
