
# Expected squared residuals.
#
get_ER2 = function (  susiF.obj,Y, X) {
  postF <- get_post_F(susiF.obj )# J by N matrix
  Xr_L = t(X%*% postF)
  postF2 <- get_post_F2(susiF.obj ) # Posterior second moment.
  return(sum((Y - X%*%postF )^2)  -sum(postF)^2 + sum(postF2))
}



#' @title Compute KL divergence effect l
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param Y Matrix of outcomes
#'
#' @param X Matrix of covariates
#'
#' @return susiF object
#' @export
update_KL <- function    (susiF.obj, Y, X, ...)
  UseMethod("update_KL")



#' @rdname update_KL
#'
#' @method update_KL susiF
#'
#' @export update_KL.susiF
#'
#' @export
#'

update_KL.susiF <- function    (susiF.obj, Y, X, ...)
{
  susiF.obj$KL <-  do.call(c,lapply(1:susiF.obj$L,FUN=function(l) cal_KL_l(susiF.obj,l, Y, X )))
  return( susiF.obj)
}



#' @title Compute KL divergence effect l
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l effect to update
#'
#' @param Y Matrix of outcomes
#'
#' @param X Matrix of covariates
#'
#' @return susiF object
#' @export
cal_KL_l <- function    (susiF.obj, l, Y, X, ...)
  UseMethod("cal_KL_l")




#' @rdname cal_KL_l
#'
#' @method cal_KL_l susiF
#'
#' @export cal_KL_l.susiF
#'
#' @export
#'

cal_KL_l.susiF <- function(susiF.obj, l, Y, X,...)
{
  out <-  - loglik_SFR(susiF.obj, l,Y,X)+ loglik_SFR_post(susiF.obj, l,Y,X)#loglik_SFR_post should be using partial results insted of Y
  return(out)
}

#' @title Compute log likelihood of single function regression of effect l
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l effect to update
#'
#' @param Y Matrix of outcomes
#'
#' @param X Matrix of covariates
#'
#' @return The log-likelihood, \eqn{\log p(Y | X, V)}, where V is the prior parameters
#' @export
loglik_SFR <- function    (susiF.obj, l,  ...)
  UseMethod("loglik_SFR")




#' @rdname loglik_SFR
#'
#' @method loglik_SFR susiF
#'
#' @export loglik_SFR.susiF
#'
#' @export
#'

loglik_SFR.susiF <- function(susiF.obj, l,Y,X)
{
  lBF <- get_lBF(susiF.obj,l)
  prior_weights <- rep(1/ncol(X),ncol(X))
  maxlBF <- max(lBF)
  w = exp( lBF- maxlBF)
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)

  lBF_model = maxlBF + log(weighted_sum_w)
  loglik <- lBF_model + sum(dnorm(Y,0,sqrt(susiF.obj$sigma2),log = TRUE))

  return(loglik)
}




#' @title Compute posterior expected loglikelihood for  single function regression of effect l
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param l effect to update
#'
#' @param Y Matrix of outcomes
#'
#' @param X Matrix of covariates
#'
#' @return posterior expected loglikelihood for  single function regression of effect l
#' @export
loglik_SFR_post <- function    (susiF.obj, l,  ...)
  UseMethod("loglik_SFR_post")


#' @rdname loglik_SFR_post
#'
#' @method loglik_SFR_post susiF
#'
#' @export loglik_SFR_post.susiF
#'
#' @export
#'

loglik_SFR_post.susiF <- function(susiF.obj, l,Y,X)
{
  n <- nrow(Y)
  t <- nrow(Y)
  EF  <- get_post_F(susiF.obj,l)
  EF2 <- get_post_F2(susiF.obj,l)
  s2  <- susiF.obj$sigma2
  return(-0.5*n*t*log(2*pi*s2) - 0.5/s2*(sum(Y*Y)- 2*sum(Y*X%*%EF)+ sum(attr(X,"d") * EF2)))
}


#' @title Expected log likelihood for a   susiF   object
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param Y Matrix of outcomes
#'
#' @param X Matrix of covariates
#'
#' @return Expected log likelihood
#' @export
Eloglik <- function    (susiF.obj, Y, X,  ...)
  UseMethod("Eloglik")


#' @rdname Eloglik
#'
#' @method Eloglik susiF
#'
#' @export Eloglik.susiF
#' @export
Eloglik.susiF = function (susiF.obj,Y ,X) {
  n = nrow(Y)
  t = ncol(Y)

  return(-(n*t/2) * log(2*pi*susiF.obj$sigma2) - (1/(2*susiF.obj$sigma2)) * get_ER2( susiF.obj, Y, X))
}



#' @title Get objective function from data and susiF object
#'
#' @param susiF.obj a susiF object defined by \code{\link{init_susiF_obj}} function
#'
#' @param Y Matrix of outcomes
#'
#' @param X Matrix of covariates
#'
#' @return posterior expected log likelihood for  single function regression of effect l
#' @export
get_objective <- function    (susiF.obj, Y, X,  ...)
  UseMethod("get_objective")


#' @rdname get_objective
#'
#' @method get_objective susiF
#'
#' @export get_objective.susiF
#' @export
get_objective.susiF <- function    (susiF.obj, Y, X,  ...)
{
  susiF.obj <- update_KL(susiF.obj, Y, X)
  out <- Eloglik(susiF.obj, Y, X) - sum(susiF.obj$KL)
  return(out)

}
