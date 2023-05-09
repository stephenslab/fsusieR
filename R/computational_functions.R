################################## susiF computational FUNCTIONS ############################


#' @title Compute correlation between credible sets
#
#' @param susiF.obj a susiF object
#' @param  X matrix of covariates
#'
#' @importFrom stats median
#' @importFrom stats cor
#'
#' @export
#' @keywords internal
cal_cor_cs <- function(susiF.obj,X){

  if(length(susiF.obj$cs)==1)
  {return(susiF.obj)}else{

    mat_cor <- matrix(NA, ncol = length(susiF.obj$cs),
                      nrow= length(susiF.obj$cs))

    for ( l in 1:length(susiF.obj$cs)){

      for ( k in 1:length(susiF.obj$cs)){

        temp <- max(do.call( c,
                             sapply( 1:length(susiF.obj$cs[[l]]),
                                     function(l1)
                                       sapply( 1:length(susiF.obj$cs[[k]]),
                                               function( k1)  cor(X[,c(susiF.obj$cs[[l]][l1],
                                                                       susiF.obj$cs[[k]][k1])])[2,1],
                                               simplify=FALSE)
                             )
        )
        )
        mat_cor[l,k] <-temp
        mat_cor[k,l] <-temp

      }
    }
  }
  susiF.obj$cs_cor <- mat_cor
  return(susiF.obj)
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

cal_Bhat_Shat   <- function(Y, X ,v1 ,resid_var=1, lowc_wc=NULL,
                            ind_analysis, cor_small2=FALSE, ...  )
{



    if(missing(ind_analysis)){


      d <- colSums(X^2)
      Bhat <- (t(X)%*%Y )/d

      Shat  <- do.call( cbind,
                        lapply( 1:ncol(Bhat),
                                function(i) matrixStats::colSds(Y [,i] -sweep( X,2, Bhat[,i], "*"))
                        )
      )
      Shat <- Shat/sqrt(nrow(Y))

    }else{
      if( is.list(ind_analysis) ){ #usefull for running multiple univariate regression with different problematic ind


        Bhat <-  do.call(cbind,lapply(1:length(ind_analysis),
                                      function(l){
                                        d   <- colSums(X[ind_analysis[[l]], ]^2)
                                        out <- (t(X[ind_analysis[[l]], ])%*%Y[ind_analysis[[l]], l])/d
                                        return(out)
                                      }


        ) )


        Shat  <-   matrix(mapply(function(l,j)
                                           sqrt(fast_var(Y[ind_analysis[[l]],l] - X[ind_analysis[[l]], j]  *  Bhat[j,l]) /(length(ind_analysis[[l]])-1)),
                                 l=rep(1:dim(Y)[2],each= ncol(X)),
                                 j=rep(1:dim(X)[2], ncol(Y))
                                 ),
                           ncol=dim(Y)[2]
                          )


      }else{
        d <- colSums(X[ind_analysis , ]^2)
        Bhat <- (t(X[ind_analysis , ])%*%Y[ind_analysis , ])/d

        Shat  <- do.call( cbind,
                          lapply( 1:ncol(Bhat),
                                  function(i) matrixStats::colSds(Y [ind_analysis,i] -sweep( X,2, Bhat[ind_analysis,i], "*"))
                          )
        )
        Shat <- Shat/sqrt(nrow(Y[ind_analysis,]))

      }

    }



    # sd_res <- sqrt(resid_var)
  #
   if( cor_small2){
    Z <- Bhat/ Shat# (sqrt(n*Shat))
  # if(missing(ind_analysis)){
      n <- nrow(Y)-1
      p <-   2 * pt(abs(Z ), df=n ,
                    lower.tail = FALSE)
    Bhat  <-  sign(Z)*stats::qnorm(p / 2, sd=c(Se) ,lower.tail=FALSE)
        # for ( i in 1: nrow(Shat)){
        #   for ( j in 1:ncol(Shat)){
        #     Shat[i,j]<-  pval2se(bhat=Bhat[i,j], p=p[i,j])
        #  }
      # }

  #}else{
  #
  #
  #   if ( is.list(ind_analysis)){
  #
  #    }else{

  #      }
  #
  #      }
       }
    if( !is.null(lowc_wc)){
      Bhat[,lowc_wc] <- 0
      Shat[,lowc_wc] <- 1
    }
    out  <- list( Bhat = Bhat,
                  Shat = Shat)

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

# @title Compute conditional local false sign rate
#
# @description Compute conditional local false sign rate
# @param G_prior mixture normal prior
#
# @param Bhat  matrix pxJ regression coefficient, Bhat[j,t] corresponds to regression coefficient of Y[,t] on X[,j]
#
# @param Shat matrix pxJ standard error, Shat[j,t] corresponds to standard error of the regression coefficient of Y[,t] on X[,j]
#
# @param indx_lst list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class mixture_normal_per_scale
#
#
#
# @return pxJ matrix of local false sign rate
#
# @export


cal_clfsr <- function (G_prior, Bhat, Shat,...)
  UseMethod("cal_clfsr")


# @rdname cal_clfsr
#
# @method cal_clfsr mixture_normal
#
# @export cal_clfsr.mixture_normal
#
#
#
# @importFrom ashr set_data
# @importFrom ashr get_fitted_g
# @export
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


# @rdname cal_clfsr
#
# @method cal_clfsr mixture_normal_per_scale
#
# @export cal_clfsr.mixture_normal_per_scale
#
#
#
# @importFrom ashr set_data
# @importFrom ashr get_fitted_g
#
# @export
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

  return(out)
}



#@title Compute conditional local false sign rate
# @description Compute conditional local false sign rate
#
# @param clfsr_wc an object of the class "clfsr_wc" generated by \code{\link{cal_clfsr}}
#
# @param alpha  vector of length P containning
#
# @return  J local false sign rate
#
# @export


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
# @importFrom stats var
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
# @importFrom ashr ash
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
#' @param \dots Other arguments.
#' @return  The log-Bayes factor for each covariate.
#'
#' @export
#' @keywords internal
log_BF <- function (G_prior, Bhat, Shat,lowc_wc, df=NULL,...)
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
log_BF.mixture_normal <- function (G_prior, Bhat, Shat,lowc_wc,df=NULL, ...) {


  if (is.null(df)){

    print(paste("df  is null"))
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
          tt <- tt + pi_k[k] * dnorm(Bhat[t,-lowc_wc],sd = sqrt(sd_k[k]^2 + Shat[t,-lowc_wc]^2))
        }
        out <- sum(log(tt) - dnorm(Bhat[t,-lowc_wc],sd = Shat[t,-lowc_wc],log = TRUE))
      }else{
        for (k in 1:length(m$fitted_g$pi))
        {
          tt <- tt + pi_k[k] * dnorm(Bhat[t,],sd = sqrt(sd_k[k]^2 + Shat[t,]^2))
        }
        out <- sum(log(tt) - dnorm(Bhat[t,],sd = Shat[t,],log = TRUE))
      }

      # tt <- ifelse(tt==Inf,max(10000, 100*max(tt[-which(tt==Inf)])),tt)

      return(out)
    }
  }else{


    print(paste("df  is not null"))
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
            t_ind <-indx_lst[[s]]

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
#' @param indx_lst list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class mixture_normal_per_scale
#'
#' @param lowc_wc wavelet coefficient with low count to be discarded
#' @param \dots Other arguments.
#' @return pxJ matrix of posterior mean
#'
#'
#'
#' @export
#' @keywords internal

post_mat_mean <- function (G_prior, Bhat, Shat,lowc_wc,...)
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
post_mat_mean.mixture_normal  <- function( G_prior ,
                                           Bhat,
                                           Shat,
                                           lowc_wc,
                                           indx_lst,...  )
{

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
                                                    lowc_wc,
                                                    indx_lst,...  )
{



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
#' @param indx_lst list generated by   gen_wavelet_indx   for the given level of resolution, used only with class mixture_normal_per_scale
#'
#' @param lowc_wc wavelet coefficient with low count to be discarded
#' @param \dots Other arguments.
#' @return pxJ matrix of posterior standard deviation
#'
#' @export
#' @keywords internal


post_mat_sd <- function (G_prior, Bhat, Shat,lowc_wc,...)
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
post_mat_sd.mixture_normal  <- function( G_prior ,
                                         Bhat,
                                         Shat,
                                         lowc_wc,
                                         indx_lst,...  )
{

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
post_mat_sd.mixture_normal_per_scale <-  function( G_prior ,
                                                   Bhat,
                                                   Shat,
                                                   lowc_wc,
                                                   indx_lst,...  )
{


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
  return(out)
}

