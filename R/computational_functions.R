################################## susiF computational FUNCTIONS ############################


#' @title Compute correlation between credible sets
#
#' @param obj a susiF object
#' @param  X matrix of covariates
#'
#' @importFrom stats median
#' @importFrom stats cor
#'
#' @export
#' @keywords internal
cal_cor_cs <- function(obj,X){

  if(length(obj$cs)==1)
  {return(obj)}else{

    mat_cor <- matrix(NA, ncol = length(obj$cs),
                      nrow= length(obj$cs))

    for ( l in 1:length(obj$cs)){

      for ( k in 1:length(obj$cs)){

        temp <- max(do.call( c,
                             sapply( 1:length(obj$cs[[l]]),
                                     function(l1)
                                       sapply( 1:length(obj$cs[[k]]),
                                               function( k1)  cor(X[,c(obj$cs[[l]][l1],
                                                                       obj$cs[[k]][k1])])[2,1],
                                               simplify=FALSE)
                             )
        )
        )
        mat_cor[l,k] <-temp
        mat_cor[k,l] <-temp

      }
    }
  }
  obj$cs_cor <- mat_cor
  return(obj)
}
# @title Regress Y nxJ on X nxp
#
# @description regression coefficients (and sd) of the column wise regression
#
# @param Y  functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelet which assume that J is of the form J^2. If J not a power of 2, susif internally remaps the data into grid of length 2^J
#
# @param X matrix of size n by p in
#
# @param v1 vector of ones of length n
#
# @param lowc_wc wavelet coefficient with low count to be discarded
#
# @param ind_analysis, optional, specify index for the individual to be analysied, allow analyis data with different entry with NA
# if a vector is provided, then we assume that the entry of Y have NA at the same place, if a list is provide
# @return list of two
#
# \item{Bhat}{ matrix pxJ regression coefficient, Bhat[j,t] corresponds to regression coefficient of Y[,t] on X[,j] }
#
# \item{Shat}{ matrix pxJ standard error, Shat[j,t] corresponds to standard error of the regression coefficient of Y[,t] on X[,j] }
#
# @export
#' @importFrom Rfast colsums
#' @importFrom Rfast colVars
#' @importFrom Rfast cova
cal_Bhat_Shat   <- function(Y,
                            X ,
                            v1 ,
                            resid_var=1,
                            lowc_wc=NULL,
                            ind_analysis=NULL,  ...  ){

    p <- ncol(X)
    J <- ncol(Y)

    Bhat <- matrix(0.0, p, J)
    Shat <- matrix(0.0, p, J)

    ## ============================================================
    ## Full data: vectorized (fastest possible in pure R)
    ## ============================================================
    if (is.null(ind_analysis)) {

      n <- nrow(Y)
      d <- colSums(X^2)

      Bhat <- crossprod(X, Y) / d

      Shat <- do.call(
        cbind,
        lapply(seq_len(J), function(j)
          Rfast::colVars(
            Y[, j] - sweep(X, 2, Bhat[, j], "*")
          )
        )
      )

      Shat <- sqrt(pmax(Shat, 1e-64)) / sqrt(n - 1)

      ## ============================================================
      ## Per-response subsets (ind_analysis is a list)
      ## ============================================================
    } else if (is.list(ind_analysis)) {

      for (j in seq_len(J)) {
        idx <- ind_analysis[[j]]
        nj  <- length(idx)

        Xj <- X[idx, , drop = FALSE]
        yj <- Y[idx, j]

        d <- colSums(Xj^2)
        Bhat[, j] <- crossprod(Xj, yj) / d

        scale <- 1 / sqrt(nj - 1)

        for (k in seq_len(p)) {
          xk <- Xj[, k]
          bk <- Bhat[k, j]
          r  <- yj - xk * bk
          Shat[k, j] <- sqrt(sum(r^2) / d[k]) * scale
        }
      }

      ## ============================================================
      ## Common subset (ind_analysis is a vector)
      ## ============================================================
    } else {

      idx <- ind_analysis
      ni  <- length(idx)

      Xi <- X[idx, , drop = FALSE]
      Yi <- Y[idx, , drop = FALSE]

      d <- colSums(Xi^2)
      Bhat <- crossprod(Xi, Yi) / d

      scale <- 1 / sqrt(ni - 1)

      for (j in seq_len(J)) {
        yj <- Yi[, j]
        for (k in seq_len(p)) {
          xk <- Xi[, k]
          bk <- Bhat[k, j]
          r  <- yj - xk * bk
          Shat[k, j] <- sqrt(sum(r^2) / d[k]) * scale
        }
      }
    }

    ## ============================================================
    ## Mask low-confidence coefficients if requested
    ## ============================================================
    if (!is.null(lowc_wc)) {
      Bhat[, lowc_wc] <- 0
      Shat[, lowc_wc] <- 1
    }
    Shat <- (pmax(Shat, 1e-32))
    out= list(Bhat = Bhat, Shat = Shat)
    return(out)
}




#' @title Compute likelihood for the weighted ash problem
#
#' @description Add description here.
#
#' @param lBF vector of log Bayes factors
#'
#' @param zeta assignment probabilities
#'
#' @return Likelihood value
#'
#' @export
#' @keywords internal
cal_lik <- function(lBF,zeta)
{
  out <- sum( zeta*exp(lBF - max(lBF ) ))
  return(out)
}

#' @title Compute conditional local false sign rate
#
#' @description Compute conditional local false sign rate
#' @param G_prior mixture normal prior
#
#' @param Bhat  matrix pxJ regression coefficient, Bhat[j,t] corresponds to regression coefficient of Y[,t] on X[,j]
#
#' @param Shat matrix pxJ standard error, Shat[j,t] corresponds to standard error of the regression coefficient of Y[,t] on X[,j]
#
#' @param indx_lst list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class mixture_normal_per_scale
#' @param \dots Other arguments.
#
#
#' @return pxJ matrix of local false sign rate
#
#' @export
#' @keywords internal


cal_clfsr <- function (G_prior, Bhat, Shat,...)
  UseMethod("cal_clfsr")


#' @rdname cal_clfsr
#
#' @method cal_clfsr mixture_normal
#
#' @export cal_clfsr.mixture_normal
#
#
#
#' @importFrom ashr set_data
#' @importFrom ashr get_fitted_g
#' @export
#' @keywords internal
#
cal_clfsr.mixture_normal  <- function(G_prior,
                                      Bhat,
                                      Shat,...)
{
  t_col_post <- function(t){
    m <- G_prior [[1]]
    data <-  ashr::set_data(Bhat[t,] ,Shat[t,] )
    return(calc_lfsr( m ,data))
  }

  out <- lapply(1:(dim(Bhat)[1] ),t_col_post )
  out <- t(Reduce("cbind", out))
  class(out) <- "clfsr_wc"
  return(out)
}


#' @rdname cal_clfsr
#
#' @method cal_clfsr mixture_normal_per_scale
#
#' @export cal_clfsr.mixture_normal_per_scale
#
#
#
#' @importFrom ashr set_data
#' @importFrom ashr get_fitted_g
#
#' @export
#' @keywords internal
#

cal_clfsr.mixture_normal_per_scale <- function( G_prior ,
                                                Bhat,
                                                Shat,
                                                indx_lst,...  )
{
  t_col_post <- function(t ){

    t_m_post <- function(s){
      m <- G_prior [[ s]]

      data <-  ashr::set_data(Bhat[t,indx_lst[[s]] ],
                              Shat[t, indx_lst[[s]] ]
      )
      return(calc_lfsr( m ,data))
    }
    return(unlist(lapply( c((length(indx_lst)  -1): 1,length(indx_lst)   ), #important to maintain the ordering of the wavethresh package !!!!
                          t_m_post )))
  }

  out <- lapply( 1:(dim(Bhat)[1]),
                 FUN= t_col_post)


  out <- t(Reduce("cbind", out))
  class(out) <- "clfsr_wc"
  return(out)
}



#' @title Compute conditional local false sign rate
#' @description Compute conditional local false sign rate
#'
#' @param clfsr_wc an object of the class "clfsr_wc" generated by \code{\link{cal_clfsr}}
#'
#' @param alpha  vector of length P containning
#'
#' @return  J local false sign rate
#'
#' @export
#' @keywords internal

cal_lfsr  <- function (clfsr_wc, alpha){
  out <-  c(alpha %*%clfsr_wc)
  return( out)
}

#' @title Compute assignment probabilities from log Bayes factors
#'
#' @description Add description here.
#'
#' @param lBF vector of log Bayes factors
#' @export
#' @keywords internal
cal_zeta <- function(lBF)
{
  out <- exp(lBF - max(lBF ) ) /sum( exp(lBF - max(lBF ) ))
  return(out)
}





# @title Regress column l of Y on column j of X
#
# @description Regress column l of Y on column j of X
#
# @param Y  functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelets that assume that J is of the form J^2. If J is not a power of 2, susiF internally remaps the data into a grid of length 2^J
#
# @param X matrix of size n by p in
#
# @param v1 vector of 1 of length n
#
# @param l Y column index
#
# @param j X column index
#
# @param lowc_wc wavelet coefficient with low count to be discarded
#
#
# @return vector of 2 containing the regression coefficient and standard error
#
# @export
#
fit_lm <- function( l,j,Y,X,v1, lowc_wc =NULL ,...)  ## Speed Gain
{

  if( !is.null(lowc_wc)){

    if (l %in%lowc_wc){
      return(c(0,
               1
      )
      )

    }else{


      return(  fast_lm(x=X[,j] ,
                       y= Y[,l]
      )
      )
    }

  }else{
    return(  fast_lm(X[,j] ,
                     Y[,l]
    )
    )
  }





}

# @title Fit ash of coefficient from scale s
#
# @description Add description here.
#
# @param Bhat  matrix pxJ regression coefficient, Bhat[j,t] corresponds to regression coefficient of Y[,t] on X[,j]
#
# @param Shat matrix pxJ standard error, Shat[j,t] corresponds to standard error of the regression coefficient of Y[,t] on X[,j]
#
# @param s scale of interest
#
# @param indx_lst list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class mixture_normal_per_scale
#
# @param lowc_wc wavelet coefficient with low count to be discarded
#
# @return an ash object
#
#' @importFrom ashr ash
#
# @export
#
fit_ash_level <- function (Bhat, Shat, s, indx_lst, lowc_wc,...)
{
  if( !is.null(lowc_wc)){


    t_ind <-indx_lst[[s]]
    t_ind <-  t_ind[which(t_ind %!in% lowc_wc)]


    if( length(t_ind)==0){ #create a ash object with full weight on null comp
      out <- ashr::ash(rnorm(100, sd=0.1),
                       rep(1,100),
                       mixcompdist = "normal" ,
                       outputlevel=0)
      out$fitted_g$pi  <- c(1, rep(0, (length(out$fitted_g$pi )-1) ))
    }else{
      out <- ashr::ash(as.vector(Bhat[,t_ind]),
                       as.vector(Shat[,t_ind]),
                       mixcompdist = "normal" ,
                       outputlevel=0)
    }

  }else{
    out <- ashr::ash(as.vector(Bhat[,indx_lst[[s]]]),
                     as.vector(Shat[,indx_lst[[s]]]),
                     mixcompdist = "normal" ,
                     outputlevel=0)
  }

  return(out)
}
#' @importFrom ashr ash calc_loglik
#' @importFrom ashr ash
#' @importFrom stats dnorm
#' @importFrom ashr calc_loglik
fit_hmm <- function (x,sd,
                     halfK=100,
                     mult=3,
                     smooth=FALSE,
                     thresh=0.00001,
                     prefilter=TRUE,
                     thresh_prefilter=1e-30,
                     maxiter=3,
                     max_zscore=20,
                     thresh_sd=1e-30,
                     epsilon=1e-2
){


  # deal with case where very close to zero sds
  if( length(which(sd< thresh_sd))>0){
    sd[ which(sd< thresh_sd)] <- thresh_sd
  }
  if (sum(is.na(sd))>0){
    x [ which( is.na(sd))]<- 0
    sd[ which( is.na(x))]<- 1
  }
  if(sum(!is.finite(sd))>0){
    x [which(!is.finite(sd))]=0
    sd[which(!is.finite(sd))]=1

  }


  if( length(which(abs(x/sd)> max_zscore))>0){ #avoid underflow  a z-score of 20=> pv< e-90


    sd[which(abs(x/sd)> max_zscore)] <- abs(x[which(abs(x/sd)> max_zscore)])/ max_zscore

  }

  K = 2*halfK-1
  sd=  sd
  X <-x

  pos <- seq(  0, 1 ,   length.out=halfK)

  #define the mean states
  mu <- (pos^(1/mult))*1.5*max(abs(X)) # put 0 state at
  #the firstplace
  mu <- c(mu, -mu[-1] )


  min_delta <- abs(mu[2]-mu[1])
  if( prefilter){
    tt <- apply(
      do.call( rbind, lapply(1:length(x), function( i){
        tt <-dnorm(x[i],mean = mu, sd=sd[i])
        return( tt/ sum(tt))
      } )),

      2,
      mean,na.rm=TRUE)
    temp_idx <- which(tt > thresh_prefilter)
    if( 1 %!in% temp_idx){
      temp_idx <- c(1,temp_idx)
    }
    mu <- mu[temp_idx]
    K <- length(mu)
  }

  P <- diag(0.5,K)+ matrix(epsilon, ncol=K, nrow=K) #this ensure that the HMM can only "transit via null state"
  P[1,-1] <- 0.5
  P[-1,1] <- 0.5



  pi = rep( 1/length(mu), length(mu)) #same initial guess

  emit = function(k,x,t){
    dnorm(x,mean=mu[k],sd=sd[t]   )
  }


  alpha_hat = matrix(  nrow = length(X),ncol=K)
  alpha_tilde = matrix(  nrow = length(X),ncol=K)
  G_t <- rep(NA, length(X))
  for(k in 1:K){
    alpha_hat[1, ] = pi* emit(1:K,x=X[1],t=1)
    alpha_tilde[1, ] = pi* emit(1:K,x=X[1],t=1)
  }



  # Forward algorithm
  for(t in 1:(length(X)-1)){
    m = alpha_hat[t,] %*% P

    alpha_tilde[t+1, ] = m *emit(1:K,x=X[t+1], t= t+1 )
    G_t[t+1] <- sum( alpha_tilde[t+1,])
    alpha_hat[t+1,] <-  alpha_tilde[t+1,]/ ( G_t[t+1])
  }

  beta_hat = matrix(nrow =  length(X),ncol=K)

  beta_tilde = matrix(nrow =  length(X),ncol=K)
  C_t <- rep(NA, length(X))
  # Initialize beta
  for(k in 1:K){
    beta_hat[ length(X),k] = 1
    beta_tilde [ length(X),k] = 1
  }

  # Backwards algorithm
  for(t in ( length(X)-1):1){
    emissio_p <- emit(1:K,X[t+1],t=t+1)
    beta_tilde [t, ] = apply( sweep( P,2, beta_hat[t+1,]*emissio_p ,"*" ),1,sum)
    C_t[t] <- max(beta_tilde[t,])
    beta_hat[t,] <-  beta_tilde [t, ] /C_t[t]
  }


  ab = alpha_hat*beta_hat
  prob = ab/rowSums(ab)

  #image(prob)#plot(apply(prob[,-1],1, sum), type='l')
  #plot(x)
  #lines(1-prob[,1])



  ab = alpha_hat*beta_hat
  prob = ab/rowSums(ab)


  xi <- array(0, dim = c(K, K))
  for (t in 1:(length(X) - 1)) {
    xi_t <- outer(alpha_hat[t, ], beta_hat[t+1, ] * emit(1:K, X[t+1], t+1)) * P
    xi_t <- xi_t / sum(xi_t)
    xi <- xi + xi_t
  }

  transition=matrix( 0, ncol=ncol(P), nrow=length(X))
  row_sums <- rowSums(xi)
  row_sums[row_sums == 0] <- 1  # prevent division by zero
  P <- xi / row_sums


  idx_comp <- which( apply(prob, 2, mean) >thresh )
  if ( !(1%in% idx_comp) ){ #ensure 0 is in the model

    idx_comp<- c(1, idx_comp)
  }

  ash_obj <- list()
  x_post <- 0*x


  for ( i in 2:length(idx_comp)){
    mu_ash <-mu[idx_comp[i] ]
    weight <- prob[,idx_comp[i]]

    ash_obj[[i]]  <- ash(x,sd,
                         weight=weight,
                         mode=mu_ash,
                         mixcompdist = "normal"
    )
    x_post <-  x_post +weight*ash_obj[[i]]$result$PosteriorMean

  }

  # plot( X)
  prob <-  prob[ ,idx_comp]
  iter=1

  K= length(idx_comp)
  P= P[idx_comp, idx_comp]

  while( iter <maxiter){


    alpha_hat = matrix(nrow = length(X),ncol=K)
    alpha_tilde = matrix(nrow = length(X),ncol=K)
    G_t <- rep(NA, length(X))

    data0 <-  set_data(X[1],sd[1])

    alpha_hat[1, ] <- prob[1, ]
    alpha_hat[1, ] <- prob[1, ]





    # Forward algorithm
    for(t in 1:(length(X)-1)){
      m = alpha_hat[t,] %*% P
      data0 <-  set_data(X[t],sd[t])



      alpha_tilde[t+1, ] = m  *c(dnorm(X[t+1], mean=0, sd=sd[t+1]),
                                 sapply( 2:K, function( k) exp(ashr::calc_loglik(ash_obj[[k]],
                                                                                 data0)
                                 )
                                 )
      )
      G_t[t+1] <- sum( alpha_tilde[t+1,])
      alpha_hat[t+1,] <-  alpha_tilde[t+1,]/ ( G_t[t+1])
    }

    beta_hat = matrix(nrow =  length(X),ncol=K)

    beta_tilde = matrix(nrow =  length(X),ncol=K)
    C_t <- rep(NA, length(X))
    # Initialize beta
    for(k in 1:K){
      beta_hat[ length(X),k] = 1
      beta_tilde [ length(X),k] = 1
    }

    # Backwards algorithm
    for(t in ( length(X)-1):1){

      data0 <-  set_data(X[t+1],sd[t+1])
      emissio_p <- c(dnorm(X[t+1], mean=0, sd=sd[t+1]),
                     sapply( 2:K, function( k) exp(ashr::calc_loglik(ash_obj[[k]],
                                                                     data0)
                     )
                     )
      )



      beta_tilde [t, ] = apply( sweep( P,2, beta_hat[t+1,]*emissio_p ,"*" ),1,sum)

      C_t[t] <- max(beta_tilde[t,])
      beta_hat[t,] <-  beta_tilde [t, ] /C_t[t]
    }




    ab = alpha_hat*beta_hat
    prob = ab/rowSums(ab)


    # image(prob)#plot(apply(prob[,-1],1, sum), type='l')
    #plot(x)
    #lines(1-prob[,1])

    ash_obj <- list()
    x_post <- 0*x

    for ( k in 2:K){
      # mu_ash <- sum(prob[,k]*X)/(sum(prob[,k])) #M step for the mean
      mu_ash <-mu[k ]



      weight <- prob[,k]

      ash_obj[[k]]  <- ash(x,sd,
                           weight=weight,
                           mode=mu_ash,
                           mixcompdist = "normal"
      )
      x_post <-  x_post +weight*ash_obj[[k]]$result$PosteriorMean

    }

    #Baum_Welch
    #Baum_Welch <-  function(X,sd,mu,P, prob, alpha, beta ){

    ab = alpha_hat*beta_hat
    prob = ab/rowSums(ab)


    xi <- array(0, dim = c(K, K))
    for (t in 1:(length(X) - 1)) {
      xi_t <- outer(alpha_hat[t, ], beta_hat[t+1, ] * emit(1:K, X[t+1], t+1)) * P
      xi_t <- xi_t / sum(xi_t)
      xi <- xi + xi_t
    }

    transition=matrix( 0, ncol=ncol(P), nrow=length(X))
    row_sums <- rowSums(xi)
    row_sums[row_sums == 0] <- 1  # prevent division by zero
    P <- xi / row_sums
    iter =iter +1
    #lines( x_post, col=iter)


    #print( sum(log(G_t[-1])))
  }


  lfsr_est <-prob[,1]
  for ( k in 2: K){
    lfsr_est <- lfsr_est + prob[,k]*ash_obj[[k]]$result$lfsr
  }


  out <- list( prob =prob,
               x_post = x_post,
               lfsr =  lfsr_est,
               mu= mu)




}

#'@title Compute refined estimate using HMM regression
#'
#' @description    Compute refined estimate using HMM regression
#'
#' @param obj  a susiF object
#' @param Y  matrix of responses
#' @param X matrix containing the covariates
#' @param verbose logical
#' @param fit_indval logical if set to true compute fitted value (default value TRUE)
#'
#' @param \dots Other arguments.
#' @export


HMM_regression<- function (obj,
                           Y,
                           X,
                           verbose=TRUE,
                           fit_indval=TRUE,...)
  UseMethod("HMM_regression")


#' @rdname HMM_regression
#'
#' @method HMM_regression susiF
#'
#' @importFrom stats lm
#'
#' @export HMM_regression.susiF
#'
#' @export
#'
HMM_regression.susiF <- function( obj,
                                  Y ,
                                  X ,
                                  verbose=TRUE,
                                  fit_indval=TRUE ,...
){

  if(verbose){
    print( "Fine mapping done, refining effect estimates using HMM regression")
  }
  idx <- do.call( c, lapply( 1:length(obj$cs),
                             function(l){
                               tp_id <-  which.max(obj$alpha[[l]])
                             }
  )
  )




  N <- nrow(X)
  sub_X <- data.frame (X[, idx])
  if(length(idx)> length(unique(idx))){

    sub_X[,duplicated(idx)]<- 0* sub_X[,duplicated(idx)]

  }

  fitted_trend  <- list()
  fitted_lfsr   <- list()


  tt <- lapply( 1: ncol(Y),
                function(j){

                  summary(lm(Y[,j]~ .,data=sub_X))$coefficients[ -1,c(1,2,4 )]


                }


  )

  fitted_trend <- list()
  fitted_lfsr   <- list()


  for ( l in 1: length(idx)){
    fitted_lfsr [[l]] <- rep(1 , ncol(Y))
    fitted_trend[[l]] <- rep(0 , ncol(Y))
  }


  if (length(idx) ==1){
    est  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 1]))
    sds  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 2]))
    pvs  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 3]))

    tsds <- sds#pval2se(est,pvs) # t -likelihood correction usefull to contrl lfsr in small sample size
    tsds[ which( tsds==0)]   <- sds[ which( tsds==0)]
    if (sum(is.na(tsds))>0){
      est [ which( is.na(tsds))]<- 0
      tsds[ which( is.na(tsds))]<- 1
    }


    s =  fit_hmm(x=est ,sd=tsds ,halfK=20 )
    fitted_lfsr [[1]] <- s$lfsr
    fitted_trend[[1]] <- s$x_post
  }else{
    for (  lp in 1: length(idx))
    {

      idx_cs <-  which( colnames(sub_X) %in% rownames(tt[[1]])[lp] )
      est  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ lp  ,1]))

      sds  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][lp,2]))
      pvs  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][lp,3]))
      tsds <- sds#pval2se(est,pvs) # t -likelihood correction usefull to contrl lfsr in small sample size
      tsds[ which( tsds==0)]<- sds[ which( tsds==0)]
      if (sum(is.na(tsds))>0){
        est [ which( is.na(tsds))]<- 0
        tsds[ which( is.na(tsds))]<- 1
      }
      s =  fit_hmm(x=est ,sd=tsds ,halfK=20 )
      fitted_lfsr [[idx_cs]] <- s$lfsr
      fitted_trend[[idx_cs]] <- s$x_post

    }

  }


  fitted_trend <- lapply(1:length(idx), function(l)
    fitted_trend[[l]]/obj$csd_X[idx[l]]
  )


  obj$fitted_func <- fitted_trend
  obj$lfsr_func   <- fitted_lfsr

 # if( fit_indval ){

  #   mean_Y <- attr(Y, "scaled:center")
  #   obj$ind_fitted_func <- matrix(mean_Y,
  #                                byrow=TRUE,
  #                                 nrow=nrow(Y),
  #                                 ncol=ncol(Y))+Reduce("+",
  #                                                      lapply(1:length(obj$alpha),
  #                                                            function(l)
                                                                #                                                              matrix( X[,idx[[l]]] , ncol=1)%*%
  #                                                              t(obj$fitted_func[[l]] )*(attr(X, "scaled:scale")[idx[[l]]])
  #                                                     )
#                                )

  #}

  return(obj)
}


#' @title Compute Log-Bayes Factor
#'
#' @description Compute Log-Bayes Factor
#'
#' @param G_prior Mixture normal prior.
#'
#' @param Bhat p x J matrix of regression coefficients;
#' \code{Bhat[j,t]} gives regression coefficient of \code{Y[,t]} on
#' \code{X[,j]}.
#'
#' @param Shat p x J matrix of standard errors; \code{Shat[j,t]}
#'   gives standard error of the regression coefficient of
#'   \code{Y[,t]} on {X[,j]}.
#'
#' @param indx_lst List generated by \code{\link{gen_wavelet_indx}}
#'   for the given level of resolution, used only with class
#'   \dQuote{mixture_normal_per_scale}
#'
#' @param lowc_wc wavelet coefficient with low count to be discarded
#' @param df degree of freedom, if set to NULL use normal distribution
#' @param indx_lst list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class mixture_normal_per_scale
#' @param \dots Other arguments.
#' @return  The log-Bayes factor for each covariate.
#'
#' @export
#' @keywords internal
log_BF <- function (G_prior, Bhat, Shat,lowc_wc,indx_lst, df=NULL,...)
  UseMethod("log_BF")

#' @rdname log_BF
#'
#' @importFrom stats dnorm
#'
#' @method log_BF mixture_normal
#'
#' @export log_BF.mixture_normal
#'
#' @export
#' @keywords internal
log_BF.mixture_normal <- function (G_prior,
                                   Bhat,
                                   Shat,
                                   lowc_wc,
                                   indx_lst,
                                   df=NULL, ...) {




  ## Deal with overfitted cases
  Shat[ Shat<=0 ] <- 1e-32
  if (is.null(df)){

    t_col_post <- function (t,lowc_wc) {

      m    <- G_prior[[1]]
      if(!is.null(lowc_wc)){
        tt   <- rep(0,length(Shat[t,-lowc_wc] ))
      }else{
        tt   <- rep(0,length(Shat[t,]))
      }

      pi_k <- m$fitted_g$pi
      sd_k <- m$fitted_g$sd

      # Speed Gain: could potential skip the one that are exactly zero.
      # Speed Gain: could potential skip the one that are exactly zero.


      if (!is.null(lowc_wc)) {
        tt <- Reduce("+", lapply(1:length(m$fitted_g$pi), function(k) {
          pi_k[k] * dnorm(Bhat[t, -lowc_wc], sd = sqrt(sd_k[k]^2 + Shat[t, -lowc_wc]^2))
        }))
        out <- sum(log(tt) - dnorm(Bhat[t, -lowc_wc], sd = Shat[t, -lowc_wc], log = TRUE))
      } else {
        tt <- Reduce("+", lapply(1:length(m$fitted_g$pi), function(k) {
          pi_k[k] * dnorm(Bhat[t, ], sd = sqrt(sd_k[k]^2 + Shat[t, ]^2))
        }))
        out <- sum(log(tt) - dnorm(Bhat[t, ], sd = Shat[t, ], log = TRUE))
      }


      # tt <- ifelse(tt==Inf,max(10000, 100*max(tt[-which(tt==Inf)])),tt)

      return(out)
    }
  }else{


    t_col_post <- function (t,lowc_wc) {

      m    <- G_prior[[1]]
      if(!is.null(lowc_wc)){
        tt   <- rep(0,length(Shat[t,-lowc_wc] ))
      }else{
        tt   <- rep(0,length(Shat[t,]))
      }

      pi_k <- m$fitted_g$pi
      sd_k <- m$fitted_g$sd

      # Speed Gain: could potential skip the one that are exactly zero.
      # Speed Gain: could potential skip the one that are exactly zero.

      if(!is.null(lowc_wc)){
        for (k in 1:length(m$fitted_g$pi))
        {
          tt <- tt + pi_k[k] *   LaplacesDemon::dstp(Bhat[t,-lowc_wc],tau = 1/(sd_k[k]^2 + Shat[t,-lowc_wc]^2), nu=df)
        }
        out <- sum(log(tt) - LaplacesDemon::dstp(Bhat[t,-lowc_wc],tau = 1/Shat[t,-lowc_wc]^2,nu=df,log = TRUE))
      }else{
        for (k in 1:length(m$fitted_g$pi))
        {
          tt <- tt + pi_k[k] *LaplacesDemon::dstp(Bhat[t,],tau = 1/(sd_k[k]^2 + Shat[t, ]^2), nu=df)
        }
        out <- sum(log(tt) - LaplacesDemon::dstp(Bhat[t,],tau = 1/Shat[t, ]^2,nu=df,log = TRUE))
      }

      # tt <- ifelse(tt==Inf,max(10000, 100*max(tt[-which(tt==Inf)])),tt)

      return(out)
    }
  }


  out <- lapply(1:nrow(Bhat),function(k) t_col_post(k, lowc_wc))
  lBF <- do.call(c,out)

  # Avoid extreme overflow problem when little noise is present.
  if (prod(is.finite(lBF)) == 0) {
    lBF <- ifelse(lBF==Inf,max(10000, 100*max(lBF[-which(lBF==Inf)])),lBF)
    lBF <- ifelse(lBF== -Inf,max(-10000, -100*max(lBF[-which(lBF== -Inf)])),
                  lBF)
  }

  return(lBF)
}

#' @rdname log_BF
#'
#' @method log_BF mixture_normal_per_scale
#'
#' @export log_BF.mixture_normal_per_scale
#'
#' @importFrom ashr set_data
#'
#'
#'
#' @export
#' @keywords internal
log_BF.mixture_normal_per_scale <- function (G_prior,
                                             Bhat,
                                             Shat,
                                             lowc_wc,
                                             indx_lst,
                                             df=NULL,
                                             ...) {


  ## Deal with overfitted cases
  Shat[ Shat<=0 ] <- 1e-32
  if (is.null(df)){
    t_col_post <- function (t) {
      t_s_post <- function (s) {

        if( !is.null(lowc_wc)){


          t_ind <-indx_lst[[s]]
          t_ind <-  t_ind[which(t_ind %!in% lowc_wc)]

          if( length(t_ind)==0){ #create a ash object with full weight on null comp
            return(0)
          }else{
            m <- G_prior[[s]] # Speed Gain: could potential skip the one that are exactly zero.
            data <- ashr::set_data(Bhat[t,t_ind],Shat[t,t_ind])
            return(ashr::calc_logLR(ashr::get_fitted_g(m),data))
          }
        }

        else{
          m <- G_prior[[s]] # Speed Gain: could potential skip the one that are exactly zero.
          data <- ashr::set_data(Bhat[t,indx_lst[[s]]],Shat[t,indx_lst[[s]]])
          return(ashr::calc_logLR(ashr::get_fitted_g(m),data))
        }

      }

      # NOTE: Maybe replace unlist(lapply(...)) with sapply(...).
      return(sum(unlist(lapply(1:(log2(ncol(Bhat))+1), # Important to maintain the ordering of the wavethresh package !!!!
                               t_s_post))))
    }
  }
  else{
    t_col_post <- function (t) {
      t_s_post <- function (s) {
        m <- G_prior[[s]]


        pi_k <- m$fitted_g$pi
        sd_k <- m$fitted_g$sd
        if( !is.null(lowc_wc)){


          t_ind <-indx_lst[[s]]
          t_ind <-  t_ind[which(t_ind %!in% lowc_wc)]
          tt   <- rep(0,length(Shat[t, t_ind]))

          if( length(t_ind)==0){ #create a ash object with full weight on null comp
            return(0)
          }else{
            for (k in 1:length(m$fitted_g$pi))
            {

              tt <- tt + pi_k[k] *LaplacesDemon::dstp(Bhat[t,t_ind],tau = 1/(sd_k[k]^2 + Shat[t,t_ind]^2), nu=df)
            }
            out <- sum(log(tt) - LaplacesDemon::dstp(Bhat[t,t_ind],tau = 1/Shat[t,t_ind]^2,nu=df,log = TRUE))

            return(out)
          }
        }

        else{
          # Speed Gain: could potential skip the one that are exactly zero.
          t_ind <-indx_lst[[s]]
          t_ind <-  t_ind[which(t_ind %!in% lowc_wc)]
          tt   <- rep(0,length(Shat[t, t_ind]))

          for (k in 1:length(m$fitted_g$pi))
          {
            tt <- tt + pi_k[k] *LaplacesDemon::dstp(Bhat[t,t_ind],tau = 1/(sd_k[k]^2 + Shat[t,t_ind]^2), nu=df)
          }
          out <- sum(log(tt) - LaplacesDemon::dstp(Bhat[t,t_ind],tau = 1/Shat[t,t_ind]^2,nu=df,log = TRUE))


          return(out )
        }

      }

      # NOTE: Maybe replace unlist(lapply(...)) with sapply(...).
      return(sum(unlist(lapply(1:(log2(ncol(Bhat))+1), # Important to maintain the ordering of the wavethresh package !!!!
                               t_s_post))))
    }
  }



  out <- lapply(1:nrow(Bhat),FUN = t_col_post)
  lBF <- do.call(c,out)

  if( prod(is.finite(lBF) )==0) # Avoid extreme overflow problem when little noise is present
  {
    lBF <-  ifelse(lBF==Inf,max(10000, 100*max(lBF[-which(lBF==Inf)])),lBF)
    lBF <-  ifelse(lBF== -Inf,max(-10000, -100*max(lBF[-which(lBF== -Inf)])),lBF)
  }
  return(lBF)
}

#' @title Compute posterior mean under mixture normal prior
#'
#' @description Add description here.
#'
#' @param G_prior mixture normal prior
#'
#' @param Bhat  matrix pxJ regression coefficient, Bhat[j,t] corresponds to regression coefficient of Y[,t] on X[,j]
#'
#' @param Shat matrix pxJ standard error, Shat[j,t] corresponds to standard error of the regression coefficient of Y[,t] on X[,j]
#'
#' @param lBF log BF
#'
#' @param indx_lst list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class mixture_normal_per_scale
#'
#' @param lowc_wc wavelet coefficient with low count to be discarded
#' @param e threshold value to avoid computing posterior that have low alpha value. Set it to 0 to compute the entire posterior. default value is 0.001
#' @param \dots Other arguments.
#' @return pxJ matrix of posterior mean
#'
#'
#'
#' @export
#' @keywords internal

post_mat_mean <- function (G_prior, Bhat, Shat,lBF,lowc_wc,indx_lst,
                           e,...)
  UseMethod("post_mat_mean")


#' @rdname post_mat_mean
#'
#' @method post_mat_mean mixture_normal
#'
#' @export post_mat_mean.mixture_normal
#'
#'
#'
#' @importFrom ashr set_data
#'
#' @export
#' @keywords internal
#'
post_mat_mean.mixture_normal  <- function( G_prior,
                                           Bhat,
                                           Shat,
                                           lBF=NULL,
                                           lowc_wc,
                                           indx_lst,
                                           e=0.001,...  )
{



  if( !is.null( lBF)){

    alpha  <- cal_zeta(   lBF)
    idx_c <-  which( alpha >e )
  }else{
    idx_c=NULL
  }

  if ( length(idx_c)==0|| is.null( lBF)){

    t_col_post <- function(t){
      m <- G_prior [[1]]
      data <-  ashr::set_data(Bhat[t, ] ,Shat[t,] )
      return(ashr::postmean(ashr::get_fitted_g(m),data))
    }



    out <- lapply( 1:(dim(Bhat)[1]),
                   FUN= t_col_post)

    out <- t(Reduce("cbind", out))


    if( !is.null(lowc_wc)){
      out[, lowc_wc] <-0
    }
  }else{
    t_out <- 0*Bhat



    t_col_post <- function(t){
      m <- G_prior [[1]]
      data <-  ashr::set_data(Bhat[t, ] ,Shat[t,] )
      return(ashr::postmean(ashr::get_fitted_g(m),data))
    }



    out <- lapply(idx_c,
                  FUN= t_col_post)

    out <- t(Reduce("cbind", out))



    t_out[idx_c,] <- out
    out <- t_out
    if( !is.null(lowc_wc)){
      out[, lowc_wc] <-0
    }
  }



  return(out)
}



#' @rdname post_mat_mean
#'
#' @method post_mat_mean mixture_normal_per_scale
#'
#' @export post_mat_mean.mixture_normal_per_scale
#'
#' @export
#' @keywords internal
#'
#' @importFrom ashr set_data
post_mat_mean.mixture_normal_per_scale <- function( G_prior,
                                                    Bhat,
                                                    Shat,
                                                    lBF=NULL,
                                                    lowc_wc,
                                                    indx_lst,
                                                    e=0.001,
                                                    ...  )
{



  if( !is.null( lBF)){

    alpha  <- cal_zeta(   lBF)
    idx_c <-  which( alpha >e )
  }else{
    idx_c=NULL
  }


  if ( length(idx_c)==0|| is.null( lBF)){

    t_col_post <- function(t  ){

      t_m_post <- function(s ){
        m <- G_prior [[ s]]

        data <-  ashr::set_data(Bhat[t,indx_lst[[s]] ],
                                Shat[t, indx_lst[[s]] ]
        )
        return(ashr::postmean(ashr::get_fitted_g(m),data))
      }
      return(unlist(lapply( c((length(indx_lst)  -1): 1,length(indx_lst)   ), #important to maintain the ordering of the wavethresh package !!!!
                            t_m_post  )
      )
      )
    }



    out <- lapply( 1:(dim(Bhat)[1]),
                   FUN= t_col_post)

    out <- t(Reduce("cbind", out))


    if( !is.null(lowc_wc)){
      out[, lowc_wc] <-0
    }
  }else{

    t_out <- 0*Bhat

    t_col_post <- function(t  ){

      t_m_post <- function(s ){
        m <- G_prior [[ s]]

        data <-  ashr::set_data(Bhat[t,indx_lst[[s]] ],
                                Shat[t, indx_lst[[s]] ]
        )
        return(ashr::postmean(ashr::get_fitted_g(m),data))
      }
      return(unlist(lapply( c((length(indx_lst)  -1): 1,length(indx_lst)   ), #important to maintain the ordering of the wavethresh package !!!!
                            t_m_post  )
      )
      )
    }



    out <- lapply( idx_c,
                   FUN= t_col_post)

    out <- t(Reduce("cbind", out))

    t_out[idx_c,] <- out
    out <- t_out
    if( !is.null(lowc_wc)){
      out[, lowc_wc] <-0
    }

  }


  return(out)

}




#'@title Compute posterior standard deviation under mixture normal prior
#'
#' @description Compute posterior standard deviation under mixture normal prior
#'
#' @param G_prior mixture normal prior
#'
#' @param Bhat  matrix pxJ regression coefficient, Bhat[j,t] corresponds to regression coefficient of Y[,t] on X[,j]
#'
#' @param Shat matrix pxJ standard error, Shat[j,t] corresponds to standard error of the regression coefficient of Y[,t] on X[,j]
#'
#' @param lBF log BF
#'
#' @param indx_lst list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class mixture_normal_per_scale
#'
#' @param lowc_wc wavelet coefficient with low count to be discarded
#' @param e threshold value to avoid computing posterior that have low alpha value
#' @param \dots Other arguments.
#' @return pxJ matrix of posterior standard deviation
#'
#' @export
#' @keywords internal


post_mat_sd <- function (G_prior, Bhat, Shat,lBF,lowc_wc,indx_lst,
                         e,...)
  UseMethod("post_mat_sd")

#' @rdname post_mat_sd
#'
#' @method post_mat_sd mixture_normal
#'
#' @export post_mat_sd.mixture_normal
#'
#'
#'
#' @importFrom ashr set_data
#'
#' @export
#' @keywords internal
#'
post_mat_sd.mixture_normal  <- function( G_prior,
                                         Bhat,
                                         Shat,
                                         lBF=NULL,
                                         lowc_wc,
                                         indx_lst,
                                         e=0.001,...  )
{

  if( !is.null( lBF)){

    alpha  <- cal_zeta(   lBF)
    idx_c <-  which( alpha >e )
  }else{
    idx_c=NULL
  }

  if ( length(idx_c)==0|| is.null( lBF)){
    t_col_post <- function(t){
      m <- G_prior [[1]]
      data <-  ashr::set_data(Bhat[t,  ] ,Shat[t, ] )
      return(ashr::postsd(ashr::get_fitted_g(m),data))
    }



    out <- lapply(1:(dim(Bhat)[1] ),t_col_post )


    out <- t(Reduce("cbind", out))


    if( !is.null(lowc_wc)){
      out[, lowc_wc] <-1
    }
  }else{


    t_out <- 0*Shat+1
    t_col_post <- function(t){
      m <- G_prior [[1]]
      data <-  ashr::set_data(Bhat[t,  ] ,Shat[t, ] )
      return(ashr::postsd(ashr::get_fitted_g(m),data))
    }



    out <- lapply(idx_c,t_col_post )


    out <- t(Reduce("cbind", out))

    t_out[idx_c,] <- out
    out <- t_out
    if( !is.null(lowc_wc)){
      out[, lowc_wc] <-1
    }

  }


  return(out)
}

#' @rdname post_mat_sd
#'
#' @method post_mat_sd mixture_normal_per_scale
#'
#' @export post_mat_sd.mixture_normal_per_scale
#'
#'
#'
#' @importFrom ashr set_data
#'
#' @export
#' @keywords internal
#'
post_mat_sd.mixture_normal_per_scale <-  function( G_prior,
                                                   Bhat,
                                                   Shat,
                                                   lBF=NULL,
                                                   lowc_wc,
                                                   indx_lst,
                                                   e=0.001,...  )
{


  if( !is.null( lBF)){

    alpha  <- cal_zeta(   lBF)
    idx_c <-  which( alpha >e )
  }else{
    idx_c=NULL
  }

  if ( length(idx_c)==0|| is.null( lBF)){

    t_col_post <- function(t  ){

      t_sd_post <- function(s ){
        m <- G_prior [[ s]]

        data <-  ashr::set_data(Bhat[t,indx_lst[[s]] ],
                                Shat[t, indx_lst[[s]] ]
        )
        return(ashr::postsd(ashr::get_fitted_g(m),data))
      }
      return(unlist(lapply( c((length(indx_lst)  -1): 1,length(indx_lst)   ), #important to maintain the ordering of the wavethresh package !!!!
                            t_sd_post  )
      )
      )
    }



    out <- lapply(1:(dim(Bhat)[1] ),t_col_post )


    out <- t(Reduce("cbind", out))


    if( !is.null(lowc_wc)){
      out[, lowc_wc] <-1
    }
  }else{

    t_out <- 0*Shat+1

    t_col_post <- function(t  ){

      t_sd_post <- function(s ){
        m <- G_prior [[ s]]

        data <-  ashr::set_data(Bhat[t,indx_lst[[s]] ],
                                Shat[t, indx_lst[[s]] ]
        )
        return(ashr::postsd(ashr::get_fitted_g(m),data))
      }
      return(unlist(lapply( c((length(indx_lst)  -1): 1,length(indx_lst)   ), #important to maintain the ordering of the wavethresh package !!!!
                            t_sd_post  )
      )
      )
    }



    out <- lapply(idx_c,t_col_post )


    out <- t(Reduce("cbind", out))

    t_out[idx_c,] <- out
    out <- t_out
    if( !is.null(lowc_wc)){
      out[, lowc_wc] <-1
    }


  }

  return(out)
}




#'@title Compute refined estimate using translation invariant wavelet transform
#'
#' @description  Compute refined estimate using translation invariant wavelet transform
#'
#' @param obj  a susiF object
#'
#' @param Y  matrix of responses
#'
#' @param X matrix containing the covariates
#' @param verbose logical
#' @param filter.number see wd description in wavethresh package description
#' @param family  see wd description in wavethresh package description
#' @param alpha required confidence level
#' @param \dots Other arguments.
#' @export


TI_regression <- function (obj,Y,X, verbose ,
                           filter.number , family   ,
                           alpha,
                           ...)
  UseMethod("TI_regression")


#' @rdname TI_regression
#'
#' @method TI_regression susiF
#'
#' @export TI_regression.susiF
#'
#'
#' @importFrom ashr ash
#' @importFrom wavethresh wd
#' @importFrom wavethresh convert
#' @importFrom wavethresh nlevelsWT
#' @importFrom wavethresh av.basis
#'
#' @export
#'



TI_regression.susiF <- function( obj,Y,X, verbose=TRUE,
                                 filter.number = 1, family = "DaubExPhase" ,
                                 alpha=0.01,
                                 ... ){

  if(verbose){
    print( "Fine mapping done, refining effect estimates using cylce spinning wavelet transform")
  }

  dummy_station_wd <- wavethresh::wd(Y[1,], type="station",
                                     filter.number = filter.number ,
                                     family = family)


  Y_f <- do.call(rbind, lapply(1:nrow(Y),
                               function( i) wavethresh::wd(Y[i,],
                                                           type="station",
                                                           filter.number = filter.number ,
                                                           family = family
                               )$D))
  Y_c <- do.call(rbind, lapply(1:nrow(Y),
                               function( i)  wavethresh::wd(Y[i,],
                                                            type="station",
                                                            filter.number = filter.number ,
                                                            family = family)$C))

  refined_est <- list(wd=list(),
                      wdC=list(),
                      wd2=list(),
                      fitted_func=list(),
                      fitted_var=list(),
                      idx_lead_cov = list()
  )



  for ( l in 1: length(obj$cs)){
    refined_est$wd[[l]]  <- rep( 0, ncol(Y_f))
    refined_est$wdC[[l]] <- rep( 0, ncol(Y_c))
    refined_est$wd2[[l]] <- rep( 0, ncol(Y_f))
    refined_est$wdC2[[l]] <- rep( 0, ncol(Y_c))
    refined_est$idx_lead_cov[[l]]  <- which.max(obj$alpha[[l]])
  }


  if(  length(obj$cs)==1){


      res <- cal_Bhat_Shat(Y_f, matrix(X[,refined_est$idx_lead_cov[[1]]],
                                       ncol=1))
      wd   <- rep( 0 ,length(res$Bhat))
      wd2  <- rep(0, length(res$Shat))
      temp <-   lapply(1:nrow(dummy_station_wd$fl.dbase$first.last.d),
                       function(s){
                         level <- s
                         n <- 2^wavethresh::nlevelsWT(dummy_station_wd)
                         first.last.d <-dummy_station_wd$fl.dbase$first.last.d
                         first.level  <- first.last.d[level, 1]
                         last.level   <- first.last.d[level, 2]
                         offset.level <- first.last.d[level, 3]
                         first.level  <- first.last.d[level, 1]
                         idx   <- (offset.level + 1 - first.level):(offset.level +n - first.level)
                         t_ash <- ashr::ash(c( res$Bhat[idx]), (c(res$Shat[idx])),
                                            nullweight=300,mixcompdist = "normal")

                         wd [idx] <- t_ash$result$PosteriorMean
                         wd2[idx] <- t_ash$result$PosteriorSD^2

                         out  <- list( wd,
                                       wd2)
                       }
      )


      wd  <- Reduce("+", lapply(1:length(temp), function(s) temp[[s]][[1]]))
      wd2 <-  Reduce("+", lapply(1:length(temp), function(s) temp[[s]][[2]]))





      refined_est$wd[[l]] <- wd
      refined_est$wd2[[l]]<-  wd2


      res <- cal_Bhat_Shat(Y_c, matrix(X[,refined_est$idx_lead_cov[[1]]],
                                       ncol=1))
      t_ash <-ashr::ash(c( res$Bhat),c(res$Shat), nullweight=3,mixcompdist = "normal")
      refined_est$wdC[[l]] <- t_ash$result$PosteriorMean




  }else{



      for (k in 1:5){

        for ( l in 1: length(obj$cs) ){
          par_res<-  Y_f -Reduce("+",
                                 lapply( (1: length(refined_est$idx_lead_cov))[-l],
                                         function(j)
                                           X[,refined_est$idx_lead_cov[[j]]] %*%t( refined_est$wd[[j]] )
                                 )
          )

          res <- cal_Bhat_Shat(par_res, matrix(X[,refined_est$idx_lead_cov[[l]]], ncol=1))
          wd <- rep( 0 ,length(res$Bhat))
          wd2 <- rep(0, length(res$Shat))
          temp <-   lapply(1:nrow(dummy_station_wd$fl.dbase$first.last.d),
                           function(s){
                             level <- s
                             n <- 2^wavethresh::nlevelsWT(dummy_station_wd)
                             first.last.d <-dummy_station_wd$fl.dbase$first.last.d
                             first.level <- first.last.d[level, 1]
                             last.level <- first.last.d[level, 2]
                             offset.level <- first.last.d[level, 3]
                             first.level <- first.last.d[level, 1]
                             idx <- (offset.level + 1 - first.level):(offset.level +n - first.level)
                             t_ash <- ashr::ash(c( res$Bhat[idx]), (c(res$Shat[idx])),
                                                nullweight=300,mixcompdist = "normal")

                             wd [idx] <- t_ash$result$PosteriorMean
                             wd2[idx] <- t_ash$result$PosteriorSD^2

                             out  <- list( wd,
                                           wd2)
                           }
          )


          wd <- Reduce("+", lapply(1:length(temp), function(s) temp[[s]][[1]]))
          wd2 <-  Reduce("+", lapply(1:length(temp), function(s) temp[[s]][[2]]))



          refined_est$wd[[l]] <- wd
          refined_est$wd2[[l]]<-  wd2


          par_resc<-  Y_c -Reduce("+",
                                  lapply( (1: length(refined_est$idx_lead_cov))[-l],
                                          function(j)
                                            X[,refined_est$idx_lead_cov[[j]]] %*%t( refined_est$wdC[[j]] )
                                  )
          )

          res <- cal_Bhat_Shat(par_resc, matrix(X[,refined_est$idx_lead_cov[[l]]], ncol=1))
          t_ash <- ashr::ash(c( res$Bhat),c(res$Shat), nullweight=300,mixcompdist = "normal")
          refined_est$wdC[[l]] <- t_ash$result$PosteriorMean
          refined_est$wdC2[[l]] <- t_ash$result$PosteriorSD^2


        }

      }
    }


  coeff= qnorm(1-( alpha)/2)
  obj$fitted_var =list()
  for( l in 1:length(obj$cs)){

    dummy_station_wd$C <- refined_est$wdC[[l]]
    dummy_station_wd$D <- refined_est$wd[[l]]
    mywst <- wavethresh::convert(dummy_station_wd  )
    nlevels <-wavethresh::nlevelsWT(mywst)
    refined_est$fitted_func[[l]]=  wavethresh::av.basis(mywst, level = (dummy_station_wd$nlevels-1), ix1 = 0,
                                                        ix2 = 1, filter = mywst$filter) *1/(obj$csd_X[ which.max(obj$alpha[[l]])] )
    mv.wd = wd.var(rep(0, ncol(Y)),   type = "station")
    mv.wd$D <-  (refined_est$wd2[[l]])
    mv.wd$C <-  ( refined_est$wdC2[[l]])

    obj$fitted_var[[l]]   <-  AvBasis.var(convert.var(mv.wd))*(1/(obj$csd_X[ which.max(obj$alpha[[l]])] )^2)
    obj$fitted_func[[l]] <-  refined_est$fitted_func[[l]]
    up                         <-  obj$fitted_func[[l]]+ coeff* sqrt(obj$fitted_var[[l]]) #*sqrt(obj$N-1)
    low                        <-  obj$fitted_func[[l]]- coeff*sqrt(obj$fitted_var[[l]]) #*sqrt(obj$N-1)
    obj$cred_band[[l]]   <- rbind(up, low)
    rownames(obj$cred_band[[l]]) <- c("up","low")

    names(obj$fitted_func)[l] <- paste("Effect_estimate_",l, sep = "")
    names(obj$fitted_var)[l]  <- paste("Variance_estimate_effect_",l, sep = "")
    names(obj$cred_band)[l]   <- paste("credible_band_effect_",l, sep = "")
  }

  rm( refined_est)

  return(obj)
}




#univariateHMM regression
# Y is  a matrix of observed function
# X is a 1 column matrix
univariate_HMM_regression <- function( Y,X,
                                       filter.number = 1 ,
                                       family = "DaubExPhase",
                                       alpha=0.05){


  X     <- colScale(X)
  csd_X <-   attr(X, "scaled:scale")

  tt <- lapply( 1: ncol(Y),
                function(j){

                  summary(lm(Y[,j]~ X[,1]))$coefficients[ ,c(1,2,4 )]


                }


  )


  est  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 1   ]))

  sds  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 2]))
  pvs  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 3]))
  tsds <- sds#pval2se(est,pvs) # t -likelihood correction usefull to contrl lfsr in small sample size
  tsds[ which( tsds==0)]<- sds[ which( tsds==0)]
  if (sum(is.na(tsds))>0){
    est [ which( is.na(tsds))]<- 0
    tsds[ which( is.na(tsds))]<- 1
  }
  s =  fit_hmm(x=est ,sd=tsds ,halfK=20 )
  out = list( effect_estimate=s$x_post*1/(csd_X),
              lfsr=s$lfsr)
  return(out)

}


#univariate TI regression
# Y is  a matrix of observed function
# X is a 1 column matrix
univariate_TI_regression <- function( Y,X,
                                      filter.number = 1 ,
                                      family = "DaubExPhase",
                                      alpha=0.05){
  X     <- colScale(X)
  csd_X <-   attr(X, "scaled:scale")

  dummy_station_wd <- wavethresh::wd(Y[1,], type="station",
                                     filter.number = filter.number ,
                                     family = family)


  Y_f <- do.call(rbind, lapply(1:nrow(Y),
                               function( i) wavethresh::wd(Y[i,],
                                                           type="station",
                                                           filter.number = filter.number ,
                                                           family = family
                               )$D))
  Y_c <- do.call(rbind, lapply(1:nrow(Y),
                               function( i)  wavethresh::wd(Y[i,],
                                                            type="station",
                                                            filter.number = filter.number ,
                                                            family = family)$C))

  refined_est <- list(wd=rep( 0, ncol(Y_f)),
                      wdC= rep( 0, ncol(Y_c)),
                      wd2=rep( 0, ncol(Y_f)),
                      fitted_func=list(),
                      fitted_var=list(),
                      idx_lead_cov = 1
  )



  res <- cal_Bhat_Shat(Y_f, matrix(X[,1], ncol=1))
  wd <- rep( 0 ,length(res$Bhat))
  wd2 <- rep(0, length(res$Shat))
  temp <-   lapply(1:nrow(dummy_station_wd$fl.dbase$first.last.d),
                   function(s){
                     level <- s
                     n <- 2^wavethresh::nlevelsWT(dummy_station_wd)
                     first.last.d <-dummy_station_wd$fl.dbase$first.last.d
                     first.level <- first.last.d[level, 1]
                     last.level <- first.last.d[level, 2]
                     offset.level <- first.last.d[level, 3]
                     first.level <- first.last.d[level, 1]
                     idx <- (offset.level + 1 - first.level):(offset.level +n - first.level)
                     t_ash <- ashr::ash(c( res$Bhat[idx]), (c(res$Shat[idx])),
                                        nullweight=300,mixcompdist = "normal")

                     wd [idx] <- t_ash$result$PosteriorMean
                     wd2[idx] <- t_ash$result$PosteriorSD^2

                     out  <- list( wd,
                                   wd2)
                   }
  )


  wd <- Reduce("+", lapply(1:length(temp), function(s) temp[[s]][[1]]))
  wd2 <-  Reduce("+", lapply(1:length(temp), function(s) temp[[s]][[2]]))


  refined_est$wd  <- wd
  refined_est$wd2 <-  wd2
  res <- cal_Bhat_Shat(Y_c, matrix(X[,1], ncol=1))
  t_ash <- ashr::ash(c( res$Bhat),c(res$Shat), nullweight=3,mixcompdist = "normal")
  refined_est$wdC  <- t_ash$result$PosteriorMean



  coeff= qnorm(1-( alpha)/2)



  dummy_station_wd$C <- refined_est$wdC
  dummy_station_wd$D <- refined_est$wd
  mywst <- wavethresh::convert(dummy_station_wd  )
  nlevels <-wavethresh::nlevelsWT(mywst)
  refined_est$fitted_func =  wavethresh::av.basis(mywst, level = (dummy_station_wd$nlevels-1), ix1 = 0,
                                                  ix2 = 1, filter = mywst$filter) *1/(csd_X)
  mv.wd = wd.var(rep(0, ncol(Y)),   type = "station")
  mv.wd$D <-  (refined_est$wd2 )

  fitted_var   <-  AvBasis.var(convert.var(mv.wd)) *(1/(csd_X)^2)
  fitted_func  <-  refined_est$fitted_func
  up                         <-   fitted_func + coeff* sqrt( fitted_var ) #*sqrt(obj$N-1)
  low                        <-   fitted_func - coeff*sqrt( fitted_var ) #*sqrt(obj$N-1)
  cred_band   <- rbind(up, low)


  rownames( cred_band ) <- c("up","low")
  out = list( effect_estimate=fitted_func,
              cred_band=cred_band,
              fitted_var= fitted_var)
  return(out)

}



univariate_smash_regression=  function( Y,X, alpha=0.05){


  coeff= qnorm(1-( alpha)/2)
  X     <- colScale(X)
  csd_X <-   attr(X, "scaled:scale")

  tt <- lapply( 1: ncol(Y),
                function(j){

                  summary(lm(Y[,j]~-1+X[,1]))$coefficients[ ,c(1,2,4 )]


                }


  )
  est  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 1   ]))

  sds  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 2]))
  pvs  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 3]))
  tsds <-  sds#pval2se(est,pvs) # t -likelihood correction usefull to contrl lfsr in small sample size
  tsds[ which( tsds==0)]<- sds[ which( tsds==0)]
  if (sum(is.na(tsds))>0){
    est [ which( is.na(tsds))]<- 0
    tsds[ which( is.na(tsds))]<- 1
  }

  s =  smashr::smash.gaus(x=est ,
                          sigma = tsds,#  sqrt(tsds) ,#mean(tsds),,
                          ashparam = list(optmethod="mixVBEM" ),
                          post.var = TRUE  )


 # s =  smashr::smash.gaus(x=est ,
 # #                        sigma = sd (  X %*%s$mu.est  - Y)/sqrt(nrow(Y))  ,#(tsds) ,#mean(tsds),,
 #                         ashparam = list(optmethod="mixVBEM" ),
 #                         post.var = TRUE  )





  fitted_trend   <- s$mu.est
  fitted_var <- s$mu.est.var
  fitted_func <- fitted_trend


  up                         <-   fitted_func + coeff* sqrt(fitted_var )
  low                        <-   fitted_func - coeff*sqrt(fitted_var )
  cred_band    <- rbind(up, low)

  out = list( effect_estimate=fitted_func,
              cred_band=cred_band,
              fitted_var= fitted_var)
  return(out)


}


#univariate TI regression
# Y is  a matrix of observed function
# X is a 1 column matrix
univariate_TI_regression_IS <- function( Y,X,
                                         filter.number = 1 ,
                                         family = "DaubExPhase",
                                         alpha=0.05){
  X     <- colScale(X)
  csd_X <-   attr(X, "scaled:scale")

  dummy_station_wd <- wavethresh::wd(Y[1,], type="station",
                                     filter.number = filter.number ,
                                     family = family)


  Y_f <- do.call(rbind, lapply(1:nrow(Y),
                               function( i) wavethresh::wd(Y[i,],
                                                           type="station",
                                                           filter.number = filter.number ,
                                                           family = family
                               )$D))
  Y_c <- do.call(rbind, lapply(1:nrow(Y),
                               function( i)  wavethresh::wd(Y[i,],
                                                            type="station",
                                                            filter.number = filter.number ,
                                                            family = family)$C))

  refined_est <- list(wd=rep( 0, ncol(Y_f)),
                      wdC= rep( 0, ncol(Y_c)),
                      wd2=rep( 0, ncol(Y_f)),
                      fitted_func=list(),
                      fitted_var=list(),
                      idx_lead_cov = 1
  )



  res <- cal_Bhat_Shat(Y_f, matrix(X[,1], ncol=1))
  wd <- rep( 0 ,length(res$Bhat))
  wd2 <- rep(0, length(res$Shat))





  t_ash <- ashr::ash(c( res$Bhat ), (c(res$Shat )),
                     nullweight=3,mixcompdist = "normal")

  wd <- t_ash$result$PosteriorMean
  wd2 <-t_ash$result$PosteriorSD^2


  refined_est$wd  <- wd
  refined_est$wd2 <-  wd2
  res <- cal_Bhat_Shat(Y_c, matrix(X[,1], ncol=1))
  t_ash <- ashr::ash(c( res$Bhat),c(res$Shat), nullweight=3,mixcompdist = "normal")
  refined_est$wdC  <- t_ash$result$PosteriorMean



  coeff= qnorm(1-( alpha)/2)



  dummy_station_wd$C <- refined_est$wdC
  dummy_station_wd$D <- refined_est$wd
  mywst <- wavethresh::convert(dummy_station_wd  )
  nlevels <-wavethresh::nlevelsWT(mywst)
  refined_est$fitted_func =  wavethresh::av.basis(mywst, level = (dummy_station_wd$nlevels-1), ix1 = 0,
                                                  ix2 = 1, filter = mywst$filter) *1/(csd_X)
  mv.wd = wd.var(rep(0, ncol(Y)),   type = "station")
  mv.wd$D <-  (refined_est$wd2 )

  fitted_var   <-  AvBasis.var(convert.var(mv.wd)) *(1/(csd_X)^2)
  fitted_func  <-  refined_est$fitted_func
  up                         <-   fitted_func + coeff* sqrt( fitted_var ) #*sqrt(obj$N-1)
  low                        <-   fitted_func - coeff*sqrt( fitted_var ) #*sqrt(obj$N-1)
  cred_band   <- rbind(up, low)


  rownames( cred_band ) <- c("up","low")
  out = list( effect_estimate=fitted_func,
              cred_band=cred_band,
              fitted_var= fitted_var)
  return(out)

}



#'@title Wrapper for univariate functional  regression used on susiF
#'
#' @description  Compute refined estimate using translation invariant wavelet transform
#'
#' @param Y  matrix of responses
#'
#' @param X a one column matrix containing the covariate of interest
#' @param verbose logical
#' @param filter.number see wd description in wavethresh package description
#' @param family  see wd description in wavethresh package description
#' @param alpha required confidence level
#' @param \dots Other arguments.
#' @export
#' @examples
#' library(susiF)
#'


univariate_functional_regression <- function(Y,X,
                                             method=c("TI", "HMM","smash"),
                                             filter.number = 1 ,
                                             family = "DaubExPhase",
                                             alpha=0.05){
  method          <- match.arg(method)
  if( method=="TI"){

    out  =  univariate_TI_regression( Y=Y,
                                      X=X,
                                      filter.number =filter.number  ,
                                      family = family,
                                      alpha=  alpha)
  }
  if( method=="HMM"){
    out= univariate_HMM_regression( Y=Y,
                                    X=X,
                                    filter.number =filter.number  ,
                                    family = family,
                                    alpha=  alpha)

  }
  if(method=="smash"){
    out= univariate_smash_regression ( Y=Y,
                                       X=X,
                                       alpha=  alpha
    )
  }

  return(out)
}


#'@title Compute refined estimate using translation invariant wavelet transform
#'
#' @description e Compute refined estimate using translation invariant wavelet transform
#'
#' @param obj  a susiF object
#'
#' @param Y  matrix of responses
#'
#' @param X matrix containing the covariates
#' @param verbose logical
#' @param filter.number see wd description in wavethresh package description
#' @param family  see wd description in wavethresh package description
#' @param alpha required confidence level
#' @param \dots Other arguments.
#' @export


smash_regression <- function (obj,Y,X, verbose ,
                              filter.number , family  ,
                              alpha  ,...)
  UseMethod("smash_regression")


#' @rdname smash_regression
#'
#' @method smash_regression susiF
#'
#' @export smash_regression.susiF
#'
#'
#' @importFrom ashr ash
#' @importFrom wavethresh wd
#' @importFrom wavethresh convert
#' @importFrom wavethresh nlevelsWT
#' @importFrom wavethresh av.basis
#'
#' @export
#'
smash_regression.susiF <- function(  obj,Y,X, verbose=TRUE,
                                     filter.number = 10, family = "DaubLeAsymm" ,
                                     alpha=0.95,...
){
  if(verbose){
    print( "Fine mapping done, refining effect estimates using cylce spinning wavelet transform")
  }
  idx <- do.call( c, lapply( 1:length(obj$cs),
                             function(l){
                               tp_id <-  which.max(obj$alpha[[l]])
                             }
  )
  )




  N <- nrow(X)
  sub_X <- data.frame (X[, idx])
  if(length(idx)> length(unique(idx))){

    sub_X[,duplicated(idx)]<- 0* sub_X[,duplicated(idx)]

  }

  fitted_trend  <- list()
  fitted_lfsr   <- list()


  tt <- lapply( 1: ncol(Y),
                function(j){
                  t_fit= summary(lm(Y[,j]~ .,data=sub_X))
                  cbind(t_fit$coefficients[ -1,c(1,2,4 )],
                        rep( t_fit$sigma, (nrow(t_fit$coefficients)-1))
                        )


                }


  )



  fitted_trend <- list()
  fitted_var   <- list()
  #some problem with low variance counts
  sigma = (do.call(c, lapply( 1: length(tt) ,function (j) mean (tt[[j]][  ,2]))))
  if( sum(is.na(sigma ))>0){
    sigma [which(is.na(sigma ))]= mean (sigma , na.rm=TRUE)
  }
 # sigma2=mean(sigma2)
  for ( l in 1: length(idx)){
    fitted_lfsr [[l]] <- rep(1 , ncol(Y))
    fitted_trend[[l]] <- rep(0 , ncol(Y))
  }


  if (length(idx) ==1){
    est  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 1]))
    sds  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 2]))
    pvs  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ 3]))

    tsds <- sds#pval2se(est,pvs) # t -likelihood correction usefull to contrl lfsr in small sample size
    tsds[ which( tsds==0)]   <- sds[ which( tsds==0)]
    if (sum(is.na(tsds))>0){
      est [ which( is.na(tsds))]<- 0
      tsds[ which( is.na(tsds))]<- 1
    }

     s =  smashr::smash.gaus(x=est ,
                             sigma =    (tsds),
                             ashparam = list(optmethod="mixVBEM" ),
                            post.var = TRUE  )
    #browser()
    # s =  smash_lw(noisy_signal=est ,
    #               noise_level  =   tsds)
    fitted_trend[[1]] <- s$mu.est
    fitted_var  [[1]] <- s$mu.est.var
  }else{

    for (  lp in 1: length(idx))
    {

      idx_cs <-  which( colnames(sub_X) %in% rownames(tt[[1]])[lp] )
      est  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][ lp  ,1]))

      sds  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][lp,2]))
      pvs  <- do.call(c, lapply( 1: length(tt) ,function (j) tt[[j]][lp,3]))
      tsds <- sds#pval2se(est,pvs) # t -likelihood correction usefull to contrl lfsr in small sample size
      tsds[ which( tsds==0)]<- sds[ which( tsds==0)]
      if (sum(is.na(tsds))>0){
        est [ which( is.na(tsds))]<- 0
        tsds[ which( is.na(tsds))]<- 1
      }


       s =  smashr::smash.gaus(x=est ,
                               sigma =     (tsds),
                               ashparam =list(optmethod="mixVBEM"),
                             post.var = TRUE  )
      # s =  smash_lw(noisy_signal=est ,
      #                noise_level  =   tsds)
      fitted_trend[[idx_cs]]  <- s$mu.est
      fitted_var[[idx_cs]]  <- s$mu.est.var
    }

  }


  fitted_trend <- lapply(1:length(idx), function(l)
    fitted_trend[[l]]/obj$csd_X[idx[l]]
  )
  fitted_var <- lapply(1:length(idx), function(l)
    fitted_var[[l]]/(obj$csd_X[idx[l]]^2)
  )


  obj$fitted_func <- fitted_trend
  obj$fitted_var  <-  fitted_var

  coeff= qnorm(1-(1-alpha)/2)
  for( l in 1:length(obj$cs)){



    up                         <-  obj$fitted_func[[l]]+ coeff* sqrt(fitted_var[[l]])
    low                        <-  obj$fitted_func[[l]]- coeff*sqrt(fitted_var[[l]])
    obj$cred_band[[l]]   <- rbind(up, low)
    names(obj$cred_band[[l]]) <- c("up","low")
    names(obj$cred_band)[l] <- paste("credible_band_effect_",l, sep = "")

  }


  return(obj)
}
smash_lw  <- function(noisy_signal, noise_level = 1, n.shifts = 50,
                      family = "DaubExPhase", filter.number = 10 ) {
  library(wavethresh)
  library(ashr)

  x <- noisy_signal
  n <- length(x)
  if( length(noise_level)==1){

    sds <- rep(noise_level, n)
  }else{
    sds=noise_level
  }
  if ((log2(n) %% 1) != 0) {
    next_pow2 <- 2^ceiling(log2(n))
    extra <- next_pow2 - n
    reflect_part <- rev(x[1:extra])
    x_padded <- c(x, reflect_part)
    sds_padded <- c(sds, rev(sds[1:extra]))
  } else {
    x_padded <- x
    sds_padded <- sds
  }

  x_reflect <- c(x_padded, rev(x_padded))
  s_reflect <- c(sds_padded, rev(sds_padded))
  pos_interest <- 1:n
  pos_interest_padded <- (length(x_padded) + 1):(length(x_padded) + n)

  n_r <- length(x_reflect)
  k <- floor(n / n.shifts)

  W <- GenW(n = length(x_reflect), filter.number = filter.number, family = family)
  wavelet_var <- apply(W^2, 1, function(row) sum(row * (s_reflect^2)))

  est <- list()
  est_var <- list()

  for (i in 1:n.shifts) {
    shifted_x <- c(x_reflect[(i * k + 1):n_r], x_reflect[1:(i * k)])
    shifted_var <- wavelet_var
    wd_shifted <- wavethresh::wd(shifted_x, filter.number = filter.number, family = family)

    idx_wave <- gen_wavelet_indx(log2(length(shifted_x)))

    temp <- lapply(1:(length(idx_wave) - 1), function(s) {
      t_ash <- ashr::ash(wd_shifted$D[idx_wave[[s]]],
                         sqrt(shifted_var[idx_wave[[s]]]),
                         nullweight = 300)
      list(wd = t_ash$result$PosteriorMean,
           wd2 = t_ash$result$PosteriorSD^2)
    })

    d <- rep(0, length(wd_shifted$D))
    d2 <- rep(0, length(wd_shifted$D)  +1 )
    for (s in 1:(length(idx_wave) - 1)) {
      d[idx_wave[[s]]] <-  temp[[s]]$wd
      d2[idx_wave[[s]]] <- temp[[s]]$wd2
    }

    wd_shifted$D <- d
    tt <- wavethresh::wr(wd_shifted)

    #mv.wd = wd.var(rep(0,length(noisy_signal)),   type = "station")
    # mv.wd$D=d2
    #ab_var=   AvBasis.var(convert.var( mv.wd))
    s_var <- apply(W^2, 1, function(row) sum(row * d2))

    est[[i]] <- tt
    est_var[[i]] <- s_var
  }

  recover_x <- function(shifted_x, i, k) {
    n <- length(shifted_x)
    shift_back <- (n - i * k) %% n
    c(shifted_x[(shift_back + 1):n], shifted_x[1:shift_back])
  }

  est_f <- lapply(1:length(est), function(i) recover_x(est[[i]], i, k = k))
  est_v_f <- lapply(1:length(est_var), function(i) recover_x(est_var[[i]], i, k = k))

  est_f <- do.call(rbind, est_f)
  est_v_f <- do.call(rbind, est_v_f)

  est_final <- list()
  est_var_final <- list()

  for (i in 1:nrow(est_f)) {
    est_final[[i]] <- 0.5 * (est_f[i, pos_interest] + rev(est_f[i, pos_interest_padded]))
    est_var_final[[i]] <- 0.5 * (est_v_f[i, pos_interest] + rev(est_v_f[i, pos_interest_padded]))
  }

  est_mean <- colMeans(do.call(rbind, est_final))
  est_var <- colMeans(do.call(rbind, est_var_final))

  return(list(mu.est = est_mean, mu.est.var = est_var))
}

smash_2lw= function( noisy_signal, noise_level=1, n.shifts=50 ){

  if( n.shifts>length(noisy_signal)){
    n.shifts=length(noisy_signal)-1
  }
  family =  "DaubExPhase"
  filter.number=1
  if( length(noise_level)==1){
    noise_level=rep(noise_level, length(noisy_signal))
  }


  x=noisy_signal
  n= length(x)
  if( length(noise_level)==1){
    sds= rep( noise_level, length(x))
  }else{
    sds=noise_level
  }


  if ((log2(n) %% 1) != 0) {
    next_pow2 <- 2^ceiling(log2(n))
    extra <- next_pow2 - n
    reflect_part <- rev(x[1:extra])
    x_padded <- c(  x, reflect_part)

    sds_padded <- c( sds, rev(sds[1:extra]))

  } else {
    x_padded <- x
    sds_padded <- sds
  }



  if ((log2(n) %% 1) != 0) {
    next_pow2 <- 2^ceiling(log2(n))
    extra <- next_pow2 - n
    reflect_part <- rev(x[1:extra])
    x_padded <- c(x, reflect_part)
    sds_padded <- c(sds, rev(sds[1:extra]))
  } else {
    x_padded <- x
    sds_padded <- sds
  }

  x_reflect <- c(x_padded, rev(x_padded))
  s_reflect <- c(sds_padded, rev(sds_padded))
  pos_interest <- 1:n
  pos_interest_padded <- (length(x_padded) + 1):(length(x_padded) + n)

  n_r <- length(x_reflect)
  k <- floor(n / n.shifts)

  W <-  wavethresh::GenW(n = length(x_reflect), filter.number = filter.number, family = family)

  # Compute wavelet coefficient variances
  wavelet_var <- apply(W^2, 1, function(row) sum(row * (s_reflect^2)))

  est <- list()
  est_var <- list()


  for (i in 1:n.shifts) {
    shifted_x <- c(x_reflect[(i * k + 1):n_r], x_reflect[1:(i * k)])
    shifted_var <- wavelet_var
    wd_shifted <- wavethresh::wd(shifted_x, filter.number = filter.number, family = family)

    idx_wave <- gen_wavelet_indx(log2(length(shifted_x)))

    temp <- lapply(1:(length(idx_wave) - 1), function(s) {
      t_ash <- ashr::ash(c(wd_shifted$D[idx_wave[[s]]]), sqrt(shifted_var[idx_wave[[s]]]) ,mixcompdist = "normal")
      out <- list(wd = t_ash$result$PosteriorMean, wd2 = t_ash$result$PosteriorSD^2)
    })

    d <- rep(0, length(wd_shifted$D))
    d2 <- rep(0, length(wd_shifted$D)+1)
    for (s in 1:(length(idx_wave) - 1)) {
      d[idx_wave[[s]]] <- temp[[s]]$wd
      d2[idx_wave[[s]]] <- temp[[s]]$wd2
    }

    wd_shifted$D <- d
    tt <- wavethresh::wr(wd_shifted)  # Wavelet filter

    # Corrected posterior variance calculation
    s_var <- apply(W^2, 1, function(row) sum(row * d2))

    est[[i]] <- tt
    est_var[[i]] <- s_var
  }



  recover_x <- function(shifted_x, i, k ) {
    n <- length(shifted_x)
    shift_back <- (n  - i*k  ) %% n  # Ensure non-negative index
    recovered_x <- c(shifted_x[(shift_back + 1):n], shifted_x[1:shift_back])
    return(recovered_x)
  }

  est_f=list()
  est_v_f=list()
  for ( i in 1: length(est)){
    est_f[[i]]=recover_x(shifted_x=est[[i]],i,k=k)
    est_v_f[[i]]=recover_x(shifted_x=est_var[[i]],i,k=k)
  }

  est_f=do.call(rbind, est_f)
  est_v_f=do.call(rbind, est_v_f)

  est_final=list()
  est_var_final=list()

  for( i in 1:nrow(est_f)){
    est_final[[i]] = 0.5* ( est_f[ i ,pos_interest]+ rev(est_f[i ,pos_interest_padded]))
    est_var_final[[i]] = 0.5* ( est_v_f[ i ,pos_interest]+ rev(est_v_f[i ,pos_interest_padded]))

  }



  return( list(mu.est=  apply(do.call(rbind, est_final),2,mean) , # colMeans(do.call(rbind, est_final)),
               mu.est.var=   apply(do.call(rbind, est_var_final),2,mean) )# colMeans(do.call(rbind, est_var_final)))
  )


}
