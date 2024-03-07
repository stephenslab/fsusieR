
cal_lik_EBmvFR <- function(Lmat, tpi_k){

  if( inherits(tpi_k ,"pi_mixture_normal")){


    out <- sum( Lmat %*%tpi_k )
  }

  if( inherits(tpi_k ,"pi_mixture_normal_per_scale")){
    out <-   sum(
      do.call( c,
               lapply(1:length(Lmat),
                      function( l )  Lmat[[l]]%*%tpi_k[[l]]
               )
      )
    )
  }
  return(out)
}


#  temp is a dummy wavethresh object
#  j effect to update

cal_post_effect <- function(obj,j, temp,indx_lst){
  temp$D                     <- obj$fitted_wc[[1]][j,-indx_lst[[length(indx_lst)]]]
  temp$D                     <- temp$D* 1/(obj$csd_X[j] )
  temp$C[length(temp$C)]     <- obj$fitted_wc[[1]][j,indx_lst[[length(indx_lst)]]]
  temp$C[length(temp$C)]     <- temp$C[length(temp$C)]*( 1/(obj$csd_X[j] ))
  return( wavethresh::wr(temp))
}


#' @rdname cal_partial_resid
#
#' @method cal_partial_resid EBmvFR
#
#' @export cal_partial_resid.EBmvFR
#
#' @export
#' @keywords internal

cal_partial_resid.EBmvFR  <- function(  obj, l, X, D, C,  indx_lst,... )
{

  obj <- obj

  update_D  <-  D -   X[,-l]%*%obj$fitted_wc[[1]][-l,-indx_lst[[length(indx_lst)]]]
  update_C  <-  C  -  as.vector(X[,-l]%*%obj$fitted_wc[[1]][-l,indx_lst[[length(indx_lst)]]])


  update_Y  <- cbind(  update_D, update_C)



  return(update_Y)
}



#' @rdname estimate_residual_variance
#
#' @method estimate_residual_variance EBmvFR
#
#' @export estimate_residual_variance.EBmvFR
#
#' @export
#
estimate_residual_variance.EBmvFR <- function(obj,Y,X,... )
{
  obj <- obj
  out <-  (1/(prod(dim(Y))))*get_ER2 (obj,Y, X  )
  #TODO: not correct, need to correct bottom part
  return(out)
}


fit_effect.EBmvFR <- function( obj, j,X,D, C,  indx_lst,lowc_wc){

  Y <- cal_partial_resid(
     obj = obj,
    l         = j,
    X         = X,
    D         = D,
    C         = C,
    indx_lst  = indx_lst
  )

  MLE_wc <- cal_Bhat_Shat(Y,X= matrix(X[,j], ncol=1),
                          resid_var  = obj$sigma2,
                          lowc_wc    = lowc_wc
  )
  obj <- update_effect.EBmvFR(obj,
                                     j          = j,
                                     MLE_wc     = MLE_wc,
                                     lowc_wc    = lowc_wc,
                                     indx_lst   = indx_lst)
  return(obj)
}

#' @export
#' @keywords internal


get_G_prior.EBmvFR <- function( obj, ...){

  out <- obj$G_prior

  return(out)
}



#' @rdname get_ER2
#
#' @method get_ER2 EBmvFR
#
#' @export get_ER2.EBmvFR
#
#' @export
#' @keywords internal

get_ER2.EBmvFR = function (   obj,Y, X,  ...) {


  obj <- obj
  postF  <- obj$fitted_wc[[1]]# J by N matrix
  postF2 <- obj$fitted_wc2[[1]] # Posterior second moment.

  return(sum(t((Y - X%*%postF ))%*%(Y - X%*%postF ) )  -sum(t(postF)%*%postF) + sum(    postF2))
}


#' @title Create an EBmvFR object
#' @details Create an EBmvFR object
#' @param G_prior a prior generated vi the init_prior function
#' @param Y matrix of wavelet coefficients
#' @param X matrix of covariates
#' @param \dots Other arguments.
#' @export

init_EBmvFR_obj <- function( G_prior, Y,X,... )
{


  MLE_wc          <- list()
  MLE_wc2         <- list()
  fitted_wc       <- list()
  fitted_wc2      <- list()
  ind_fitted_func <- matrix(0, nrow = dim(Y)[1], ncol=dim(Y)[2]  )
  est_pi          <- list()
  est_sd          <- list()
  G_prior         <- G_prior
  N               <- dim(Y)[1]
  n_wac           <- dim(Y)[2]
  P               <- dim(X)[2]
  sigma2          <- mean(apply(Y,2 ,var))
  KL              <- rep(NA,ncol(X))
  ELBO            <- c()
  mean_X          <- attr(X, "scaled:center")
  csd_X           <- attr(X, "scaled:scale")
  d               <- attr(X , "d")

  pi_hist         <- list()

    MLE_wc    [[1]]       <-  matrix(0, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    MLE_wc2   [[1]]       <-  matrix(1, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    fitted_wc [[1]]       <-  matrix(0, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    fitted_wc2[[1]]       <-  matrix(1, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    est_pi    [[1]]       <-  get_pi_G_prior(G_prior)
    est_sd    [[1]]       <-  get_sd_G_prior(G_prior)
    pi_hist   [[1]]       <- get_pi_G_prior(G_prior)

  obj <- list( MLE_wc          = MLE_wc,
               MLE_wc2         = MLE_wc2,
               fitted_wc       = fitted_wc,
               fitted_wc2      = fitted_wc2,
               KL              = KL,
               ELBO            = ELBO,
               ind_fitted_func = ind_fitted_func,
               G_prior         = G_prior,
               N               = N,
               n_wac           = n_wac,
               sigma2          = sigma2,
               P               = P,
               est_pi          = est_pi,
               est_sd          = est_sd,
               csd_X           = csd_X,
               d               = d,
               pi_hist         = pi_hist)

  class(obj) <- "EBmvFR"
  return(obj)
}


#' @rdname test_stop_cond
#
#' @method test_stop_cond EBmvFR
#
#' @export test_stop_cond.EBmvFR
#
#' @export
#' @keywords internal
test_stop_cond.EBmvFR <- function(obj, check, cal_obj, Y, X, D, C, indx_lst,...)
{


     obj <- obj

    if( cal_obj){

      stop("use cal_obj=FALSE, ELBO currently bieng implemented")
      obj <- update_KL(obj,
                             X,
                             D= D,
                             C= C , indx_lst)

      obj <- update_ELBO(obj,
                               get_objective( obj = obj,
                                              Y         = Y,
                                              X         = X,
                                              D         = D,
                                              C         = C,
                                              indx_lst  = indx_lst
                               )
      )

      if(length(obj$ELBO)>1    )#update parameter convergence,
      {
        check <- abs(diff(obj$ELBO)[(length( obj$ELBO )-1)])
        obj$check <- check
        return(obj)
      }else{
        obj$check <- check
        return(obj)
      }
    }
    else{
      len <- length( obj$pi_hist)
      if( len>1)#update parameter convergence, no ELBO for the moment
      {
        check <-0

        tpi_2 <- do.call( c, obj$pi_hist[[len ]])
        tpi_1 <- do.call( c, obj$pi_hist[[(len-1) ]])

        #print( sum(abs(tpi_2-tpi_1))/log(length(tpi_1)))
        check <- sum(abs(tpi_2-tpi_1))/log(length(tpi_1))#adding the log of length to limit computational time
        obj$check <- check
        return(obj)
        #print(check)
      }else{
        obj$check <- check
        return(obj)
      }
    }
  return(obj)

}


#' @rdname out_prep
#
#' @method out_prep EBmvFR
#
#' @export out_prep.EBmvFR
#
#' @export
#' @keywords internal

out_prep.EBmvFR <- function( obj,  Y,  X, indx_lst,    outing_grid,...)
{
  obj <- obj

  obj             <-  update_cal_fit_func(obj, indx_lst)
  obj             <-  update_cal_indf    (obj, X = X)
  obj$outing_grid <-  outing_grid
  return(obj)
}






update_pi_hist <- function(obj,tpi ){
  obj$pi_hist[[length(obj$pi_hist)+1]] <- tpi
  return(obj)
}







#' @rdname update_cal_fit_func
#
#' @method update_cal_fit_func EBmvFR
#
#' @export update_cal_fit_func.EBmvFR
#
#' @export
#' @keywords internal

update_cal_fit_func.EBmvFR <- function(obj, indx_lst,...)
{
  obj <- obj
  temp <- wavethresh::wd(rep(0, obj$n_wac))


       fit_func <-   do.call(rbind, lapply(1:obj$P,
                                   function( j) cal_post_effect(obj= obj,
                                                         j,
                                                         temp=temp,
                                                         indx_lst=indx_lst)
                                         )
                            )


    obj$fitted_func <-  fit_func

  return(obj)


}




#  @rdname update_cal_indf
#
#  @method update_cal_indf EBmvFR
#
#  @export update_cal_indf.EBmvFR
#
#' @export
#' @keywords internal
update_cal_indf.EBmvFR  <- function(obj,  Y, X, indx_lst, TI, ...){

  obj$ind_fitted_func <-   X%*%obj$fitted_func
  return(obj)
}




update_effect.EBmvFR <- function(obj, j, MLE_wc,lowc_wc,indx_lst){
  obj$MLE_wc[[1]][j,]     <- MLE_wc$Bhat
  obj$MLE_wc2[[1]][j,]    <- MLE_wc$Shat
  obj$fitted_wc[[1]][j,]  <- post_mat_mean( G_prior   = get_G_prior(obj),
                                                   Bhat     = MLE_wc$Bhat,
                                                   Shat     = MLE_wc$Shat,
                                                   lowc_wc  = lowc_wc,
                                                   indx_lst = indx_lst)
  obj$fitted_wc2[[1]][j,] <- post_mat_sd( G_prior    = get_G_prior(obj),
                                                 Bhat     = MLE_wc$Bhat,
                                                 Shat     = MLE_wc$Shat,
                                                 lowc_wc  = lowc_wc,
                                                 indx_lst = indx_lst)^2


  return(obj)

}

 
update_prior_EBmvFR <- function( obj,
                                 max_step = 100,
                                 espsilon = 0.0001,
                                 init_pi0_w =1,
                                 control_mixsqp,
                                 indx_lst,
                                 lowc_wc,
                                 nullweight  ){


  get_G_prior(obj)


  #static parameters
  Lmat  <-  L_mixsq(G_prior   = get_G_prior(obj),
                    Bhat      = obj$MLE_wc[[1]],
                    Shat      = sqrt(obj$MLE_wc2[[1]]),
                    indx_lst  = indx_lst,
                    is.EBmvFR =TRUE)

  J <- dim(obj$MLE_wc[[1]])[1]
  tsd_k <- get_sd_G_prior(get_G_prior( obj))

  #dynamic parameters
  tpi_k = get_pi_G_prior(get_G_prior( obj))
  oldloglik <-0
  newloglik <-1

  zeta <- rep(1/J,J) #assignation initial value
  k <- 1 #counting the number of iteration
  tpi_k <- get_pi_G_prior(get_G_prior(obj))

  while( k <=max_step &  abs(newloglik-oldloglik)>=espsilon)
  {
    # E step----
    oldloglik <- cal_lik_EBmvFR(Lmat,tpi_k)
    #print(paste("EM oldloglik is ",  oldloglik))

    # M step ----
    tpi_k   <- m_step(Lmat,zeta,indx_lst,
                      init_pi0_w     = init_pi0_w,
                      control_mixsqp = control_mixsqp,
                      nullweight     = nullweight,
                      is.EBmvFR = TRUE)

    newloglik <- cal_lik_EBmvFR(Lmat,tpi_k)
    #print(paste("EM  newloglik is ",   newloglik))
    k <- k+1

  }

  obj         <-  update_pi_hist(obj,
                                        tpi = tpi_k
  )
  obj$G_prior <-  update_prior(get_G_prior( obj) ,
                                      tpi = tpi_k
  )

  return( obj)

}





#' @rdname  update_residual_variance
#
#' @method  update_residual_variance EBmvFR
#
#' @export  update_residual_variance.EBmvFR
#
#' @export
#' @keywords internal

 update_residual_variance.EBmvFR <- function( obj, sigma2  ,... ){
   obj <- obj
   obj$sigma2 <- sigma2
   return(obj)
 }
