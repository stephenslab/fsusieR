





HMM_regression.susiF <- function( susiF.obj,Y,X , maxit=5   ){


  idx <- do.call( c, lapply( 1:length(susiF.obj$cs),
                             function(l){
                               tp_id <-  which.max( susiF.obj$pip[susiF.obj$cs[[l]]])
                               susiF.obj$cs[[l]][tp_id]
                             }
  )
  )

  temp_Y <- Y
  fitted_trend <- list()
  fitted_lfdr   <- list()
  for ( k in 1:maxit){
    for (l in 1:length(idx)){
      res <- cal_Bhat_Shat(temp_Y,X )

      s = fit_hmm(x=res$Bhat[idx[l],],sd=res$Shat[idx[l],],halfK=100 )
      fitted_lfdr [[l]] <- s$prob[,1]
      fitted_trend[[l]] <- s$x_post
      if( j ==length(idx)){
        idx_var <- (1:length(idx)) [- (1)]
      }else{
        idx_var <- (1:length(idx))[- (l+1)]
      }


      temp_Y <- Y - Reduce("+", lapply( idx_var, function( j){
        X[,idx[l] ]%*%t(fitted_trend[[l]])
      }

      )
      )

    }

  }

  fitted_trend <- lapply(1:length(idx), function(l)
    fitted_trend[[l]]/susiF.obj$csd_X[idx[l]]
  )


 susiF.obj$fitted_func <- fitted_trend
 susiF.obj$lfsr_func   <- fitted_lfdr

 mean_Y          <- attr(Y, "scaled:center")
 susiF.obj$ind_fitted_func <- matrix(mean_Y,
                                     byrow=TRUE,
                                     nrow=nrow(Y),
                                     ncol=ncol(Y))+Reduce("+",
                                                          lapply(1:length(susiF.obj$alpha),
                                                                 function(l)
                                                                   matrix( X[,idx[[l]]] , ncol=1)%*%  t(susiF.obj$fitted_func[[l]] )
                                                          )
                                     )

  return(susiF.obj)
}
