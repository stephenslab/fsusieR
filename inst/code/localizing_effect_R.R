iter_trend_filtering <- function(susiF.obj, Y,X,L_points=20){

  idx <- do.call( c, lapply( 1:length(susiF.obj$cs),
                             function(l){
                               tp_id <-  which.max( susiF.obj$pip[susiF.obj$cs[[l]]])
                               susiF.obj$cs[[l]][tp_id]
                             }
  )
  )

  temp_Y <- Y
  fitted_trend <- list()
  for ( k in 1:3){
    for (j in 1:length(idx)){
      res <- cal_Bhat_Shat(temp_Y,X )
      #plot(res$Bhat[idx[j],])
      s =susieR:: susie_trendfilter(res$Bhat[idx[j],],
                                    order=0,
                                    L=L_points)
      #plot(predict(s))
      fitted_trend[[j]] <- predict(s)
      if( j ==length(idx)){
        idx_var <- (1:length(idx)) [- (1)]
      }else{
        idx_var <- (1:length(idx))[- (j+1)]
      }


      temp_Y <- Y - Reduce("+", lapply( idx_var, function( j){
        X[,idx[j] ]%*%t(fitted_trend[[j]])
      }

      )
      )

    }

  }

  fitted_trend <- lapply(1:length(idx), function(l)
    fitted_trend[[l]]/susiF.obj$csd_X[idx[l]]
  )



}
