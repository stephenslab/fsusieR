DWT2_pois <-  function (data )
{



  output <- do.call(rbind, lapply( 1:nrow(data), function (i) hft(data[i,])))

  return(output)
}

cal_Bhat_Shat2 <- function(Y, X ,v1 ,resid_var=1, lowc_wc=NULL,ind_analysis,  ...  )
{


  if(missing(ind_analysis)){

    out <- t(mapply(function(l,j)  fast_lm2(x=cbind(v1,X[,j] ),
                                           y= Y[,l]
    )
    ,
    l=rep(1:dim(Y)[2],each= ncol(X)),
    j=rep(1:dim(X)[2], ncol(Y))
    )
    )




  }else{
    if( is.list(ind_analysis) ){
      out <- t(mapply(function(l,j)  fast_lm2(x=cbind(v1[ind_analysis[[l]]],X[ind_analysis[[l]],j]  )  ,
                                             y= Y[ind_analysis[[l]],l]
      )
      ,
      l=rep(1:dim(Y)[2],each= ncol(X)),
      j=rep(1:dim(X)[2], ncol(Y))
      )
      )
    }else
      out <- t(mapply(function(l,j)  fast_lm2(x=cbind(v1[ind_analysis ],X[ind_analysis ,j] ),
                                             y= Y[ind_analysis ,l]
      )
      ,
      l=rep(1:dim(Y)[2],each= ncol(X)),
      j=rep(1:dim(X)[2], ncol(Y))
      )
      )
  }



  # sd_res <- sqrt(resid_var)
  B0hat   <-  matrix( unlist(out[,1]), nrow=ncol(X))
  Bhat   <-  matrix( unlist(out[,2]), nrow=ncol(X))
  Shat   <-  matrix( unlist(out[,3]), nrow=ncol(X))#matrix( sd_res , nrow=ncol(X), ncol=ncol(Y))/sqrt(attr(X,"d"))
  if( !is.null(lowc_wc)){
    B0hat [,lowc_wc] <- 0
    Bhat  [,lowc_wc] <- 0
    Shat  [,lowc_wc] <- 1
  }
  out  <- list( B0hat = B0hat,
                Bhat = Bhat,
                Shat = Shat)

  return(out)
}


fast_lm2 <- function(x,y)
{

  be <- solve(crossprod(x),crossprod(x,y))
  sd <-  sqrt(fast_var(y - x %*% be)/(length(x)-1))

   return(c(be ,sd))
  #return( summary(lm(y~x))$coefficients[2,1:2 ])
}



