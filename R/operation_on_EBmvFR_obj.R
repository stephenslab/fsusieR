
cal_lik_EBmvFR <- function(Lmat, tpi_k){

  if( inherits(tpi_k ,"pi_mixture_normal")){


    out <- sum(log(Lmat %*% tpi_k))
  }

  if( inherits(tpi_k ,"pi_mixture_normal_per_scale")){
    out <-   sum(
      do.call( c,
               lapply(1:length(Lmat),
                     function(l) log(Lmat[[l]] %*% tpi_k[[l]])
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

  # The workhorse maintains this cache with rank-one updates. Standalone
  # calls use the same partial residual without relying on a cached matrix.
  if (!is.null(obj$residual)) {
    return(obj$residual + tcrossprod(X[, l], obj$fitted_wc[[1]][l, ]))
  }
  cbind(D, C) - X[, -l, drop = FALSE] %*%
    obj$fitted_wc[[1]][-l, , drop = FALSE]
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
  active <- obj$active_wc
  if (is.null(active)) active <- seq_len(ncol(Y))
  floor <- obj$sigma2_floor
  if (is.null(floor)) floor <- .Machine$double.xmin
  max(floor, get_ER2(obj, Y, X) / (nrow(Y) * length(active)))
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
                          sigma2     = obj$sigma2,
                          lowc_wc    = lowc_wc
  )
  obj <- update_effect.EBmvFR(obj,
                                     j          = j,
                                     MLE_wc     = MLE_wc,
                                     lowc_wc    = lowc_wc,
                                     indx_lst   = indx_lst)
  if (!is.null(obj$residual)) {
    obj$residual <- Y - tcrossprod(X[, j], obj$fitted_wc[[1]][j, ])
  }
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


  active <- obj$active_wc
  if (is.null(active)) active <- seq_len(ncol(Y))
  residual <- obj$residual
  if (is.null(residual)) residual <- Y - X %*% obj$fitted_wc[[1]]
  # fitted_wc2 contains posterior VARIANCES, not second moments.
  sum(residual[, active, drop = FALSE]^2) +
    sum(colSums(X^2) * rowSums(obj$fitted_wc2[[1]][, active, drop = FALSE]))
}


#' @title Create an EBmvFR object
#' @details Create an EBmvFR object
#' @param G_prior a prior generated vi the init_prior function
#' @param Y matrix of wavelet coefficients
#' @param X matrix of covariates
#' @param lowc_wc Wavelet columns excluded from fitting.
#' @param \dots Other arguments.
#' @export

init_EBmvFR_obj <- function(G_prior, Y, X, ..., lowc_wc = NULL)
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
  active_wc       <- setdiff(seq_len(n_wac), lowc_wc)
  sigma2          <- mean(Y[, active_wc, drop = FALSE]^2)
  sigma2_floor    <- max(.Machine$double.xmin, .Machine$double.eps * sigma2)
  sigma2          <- max(sigma2, sigma2_floor)
  KL              <- rep(0, ncol(X))
  ELBO            <- c()
  mean_X          <- attr(X, "scaled:center")
  csd_X           <- attr(X, "scaled:scale")
  if (is.null(csd_X)) csd_X <- rep(1, P)
  d               <- colSums(X^2)
  prior_groups <- if (inherits(G_prior, "mixture_normal")) {
    list(active_wc)
  } else {
    lapply(gen_wavelet_indx(log2(n_wac)), setdiff, y = lowc_wc)
  }
  mixture_counts <- lapply(G_prior, function(g) {
    matrix(0, nrow = P, ncol = length(g$fitted_g$pi))
  })

  pi_hist         <- list()

    MLE_wc    [[1]]       <-  matrix(0, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    MLE_wc2   [[1]]       <-  matrix(1, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    fitted_wc [[1]]       <-  matrix(0, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    fitted_wc2[[1]]       <-  matrix(0, nrow = dim(X)[2], ncol=dim(Y)[2]  )
    est_pi    [[1]]       <-  get_pi_G_prior(G_prior)
    est_sd    [[1]]       <-  get_sd_G_prior(G_prior)
    pi_hist   [[1]]       <- get_pi_G_prior(G_prior)

  # Start from q = g, a valid variational distribution with zero KL.
  for (s in seq_along(G_prior)) {
    g <- G_prior[[s]]$fitted_g
    fitted_wc2[[1]][, prior_groups[[s]]] <- sum(g$pi * g$sd^2)
    mixture_counts[[s]] <- outer(rep(length(prior_groups[[s]]), P), g$pi)
  }

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
               sigma2_floor    = sigma2_floor,
               active_wc       = active_wc,
               prior_groups    = prior_groups,
               mixture_counts  = mixture_counts,
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
  if (cal_obj) {
    value <- .ebmvfr_objective(obj, Y, X)
    previous <- utils::tail(obj$ELBO, 1L)
    obj$ELBO <- c(obj$ELBO, value)
    if (length(previous)) {
      change <- value - previous
      if (change < -1e-8 * (1 + abs(previous))) {
        warning("EBmvFR ELBO decreased; inspect the fit for numerical problems")
      }
      check <- abs(change)
    }
  }
  # Without ELBO monitoring, the workhorse supplies the maximum change in
  # coefficients, mixture weights and residual variance, not weights alone.
  obj$check <- check
  obj
}

.ebmvfr_objective <- function(obj, Y, X) {
  nobs <- nrow(Y) * length(obj$active_wc)
  -nobs / 2 * log(2 * pi * obj$sigma2) -
    get_ER2(obj, Y, X) / (2 * obj$sigma2) - sum(obj$KL)
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
update_cal_indf.EBmvFR <- function(obj, Y, X, indx_lst, TI, ...) {
  # fitted_func has already been divided by the original predictor SDs.
  obj$ind_fitted_func <- sweep(X, 2, obj$csd_X, "*") %*% obj$fitted_func
  obj
}




update_effect.EBmvFR <- function(obj, j, MLE_wc, lowc_wc, indx_lst) {
  obj$MLE_wc[[1]][j, ] <- MLE_wc$Bhat
  # Preserve the existing field's standard-error units; no extra sqrt.
  obj$MLE_wc2[[1]][j, ] <- MLE_wc$Shat
  obj$fitted_wc[[1]][j, ] <- 0
  obj$fitted_wc2[[1]][j, ] <- 0
  obj$KL[j] <- 0
  for (s in seq_along(obj$G_prior)) {
    idx <- obj$prior_groups[[s]]
    if (!length(idx)) next
    post <- .ebmvfr_normal_means(MLE_wc$Bhat[1, idx],
                                MLE_wc$Shat[1, idx],
                                obj$G_prior[[s]]$fitted_g)
    obj$fitted_wc[[1]][j, idx] <- post$mean
    obj$fitted_wc2[[1]][j, idx] <- post$variance
    obj$mixture_counts[[s]][j, ] <- post$counts
    obj$KL[j] <- obj$KL[j] + post$KL
  }
  obj
}

# Normal-means posterior for a zero-centred normal mixture (including a
# point mass at zero). Compute moments, responsibilities and augmented KL
# together. Only this predictor's coefficient-by-component matrix is held
# temporarily; retain counts per predictor/scale, not a p-by-T-by-K array.
.ebmvfr_normal_means <- function(bhat, shat, g) {
  keep <- which(g$pi > 0)
  prior_var <- g$sd[keep]^2
  se2 <- shat^2
  total_var <- outer(se2, prior_var, "+")
  log_pi <- log(g$pi[keep])
  log_weight <- sweep(-0.5 * (log(2 * pi * total_var) + bhat^2 / total_var),
                      2, log_pi, "+")
  offset <- matrixStats::rowMaxs(log_weight)
  weight <- exp(log_weight - offset)
  normalizer <- rowSums(weight)
  weight <- weight / normalizer
  log_post <- log_weight - (offset + log(normalizer))

  shrink <- sweep(1 / total_var, 2, prior_var, "*")
  component_mean <- bhat * shrink
  mean <- rowSums(weight * component_mean)
  variance <- rowSums(weight * (se2 * shrink + (component_mean - mean)^2))
  counts <- numeric(length(g$pi))
  counts[keep] <- colSums(weight)

  # Gaussian KL = 0 for the degenerate zero component. This expression
  # avoids dividing by its zero prior variance or taking log(0).
  gaussian_KL <- 0.5 * (log1p(outer(1 / se2, prior_var, "*")) - shrink +
                         bhat^2 / total_var * shrink)
  KL <- sum(weight * (sweep(log_post, 2, log_pi, "-") + gaussian_KL))
  list(mean = mean, variance = variance, counts = counts, KL = KL)
}

 
update_prior_EBmvFR <- function(obj,
                               max_step = 100,
                               espsilon = 0.0001,
                               init_pi0_w = 1,
                               control_mixsqp,
                               indx_lst,
                               lowc_wc,
                               nullweight) {
  # Kim et al.'s M-step holds q fixed and averages its responsibilities.
  # The legacy optimizer arguments remain accepted for call compatibility.
  for (s in seq_along(obj$G_prior)) {
    counts <- obj$mixture_counts[[s]]
    totals <- colSums(counts)
    if (sum(totals) == 0) next  # A completely masked scale has no data.
    old_pi <- obj$G_prior[[s]]$fitted_g$pi
    new_pi <- totals / sum(totals)
    keep <- which(totals > 0)
    # q is unchanged: only its categorical prior contribution to KL moves.
    obj$KL <- obj$KL + drop(counts[, keep, drop = FALSE] %*%
                            (log(old_pi[keep]) - log(new_pi[keep])))
    obj$G_prior[[s]]$fitted_g$pi <- new_pi
  }
  tpi <- get_pi_G_prior(obj$G_prior)
  obj$est_pi[[1]] <- tpi
  update_pi_hist(obj, tpi)
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
