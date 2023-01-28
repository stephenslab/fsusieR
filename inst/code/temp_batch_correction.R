

#' @title Adjusting functions for covariates
#
#' @description Simple function to adjust the observed curves for covariates (e.g., age, population structure)
#' @param Y  observed curves.
#' @param X matrix containing the covariates
#' @param Batch vector of character specifying the Batch
#' @param pos Original position of the column of Y
#' @param thresh_lowcount numeric, used to check the wavelet coefficients have
#'  problematic distribution (very low dispersion even after standardization).
#'  Basically check if the median of the absolute value of the distribution of
#'   a wavelet coefficient is below this threshold, if yes, the algorithm discard
#'   this wavelet coefficient (setting its estimate effect to 0 and estimate sd to 1).
#'   Set to 0 by default. It can be useful when analyzing sparse data from sequence
#'    based assay or small samples.
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings are printed to the
#' console.
#' @export
#


adjust_FM_covariate <- function(Y,X=NULL,Batch,pos=NULL, thresh_lowcount=0, verbose){



  if (is.null(pos))
  {
    pos <- 1:dim(Y)[2]
  }
  if(!is.wholenumber(log2(dim(Y)[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check whether dim(Y) not equal to 2^J or if the data are unevenly spaced
  {

    inter_pol.obj <- interpol_mat(Y, pos)
    Y             <- inter_pol.obj$Y
    bp            <- inter_pol.obj$bp

    start_pos <- min( pos)
    end_pos <-max(pos)
    outing_grid   <- start_pos + (end_pos-start_pos)/(length(pos))*inter_pol.obj$grid
    if(verbose)
    {
      message( "Response matrix dimensions not equal to nx 2^J \n or unevenly spaced data \n interpolation procedure used")
    }
  }  else{

    outing_grid   <- pos
  }


  W <- DWT2(Y)
  Y_f      <-  cbind( W$D,W$C)

  lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)


  if(missing(Batch)){

    W <- DWT2(Y)
    Y_f      <-  cbind( W$D,W$C)


    v1 <- rep(1, dim(X)[2])
    lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)


    Est <- lapply( 1:ncol(Y_f), function(j){
      fit <-  lm(Y_f[,j]~ X)
      out <-  list(coef=summary(fit)$coefficients[,1] ,
                   sd  =summary(fit)$coefficients[,2])
      return(out)
    }




    )
    Coef <- t(do.call(rbind,lapply(1:length(Est), function(j) Est[[j]]$coef)))
    Sd   <- t(do.call(rbind,lapply(1:length(Est), function(j) Est[[j]]$sd)))



    indx_lst <- gen_wavelet_indx(log2(dim(Y_f)[2]))
    if (is.null(nrow(Coef))){

      coef_p <- ashr::ash(Coef,Sd )$result$PosteriorMean
      temp <- wavethresh::wd(rep(0, length(Coef)[2]))

      temp$D                     <-  coef_p[ -indx_lst[[length(indx_lst)]]]
      temp$C[length(temp$C)]     <-  coef_p[ indx_lst[[length(indx_lst)]]]
      fitted_coef   <-  wavethresh::wr(temp)
    }else{
      fitted_coef <-list()
      for ( i in 1:nrow(Coef ))
      {
        coef_p <- ashr::ash(Coef[i,],Sd[i,])$result$PosteriorMean
        temp <- wavethresh::wd(rep(0, dim(Coef)[2]))

        temp$D                     <-  coef_p[  -indx_lst[[length(indx_lst)]]]
        temp$C[length(temp$C)]     <-   coef_p[  indx_lst[[length(indx_lst)]]]
        fitted_coef[[i]]   <-  wavethresh::wr(temp)
      }
      fitted_coef <- do.call(cbind,   fitted_coef )
    }


    Y_adjusted  = Y-X%*%t(fitted_coef[,-1])#removing estimated baseline
    out <- list( Y_adjusted  = Y_adjusted,

                 fitted_coef = fitted_coef[,-1],#removing estimated baseline
                 pos         = outing_grid
    )
  }else{

    Batch <- as.factor(Batch)
    contrasts(Batch) <- contr.sum(levels(Batch))
    Batch <- model.matrix(~Batch)[,-1,drop=FALSE]
    if(!is.null(X)){

      X <- cbind(Batch,X)
    }else{
      X <- Batch
    }

    Est <- lapply( 1:ncol(Y_f), function(j){
      fit <-  lm(Y_f[,j]~ X)
      out <-  list(coef=summary(fit)$coefficients[,1] ,
                   sd  =summary(fit)$coefficients[,2])
      return(out)
    }




    )
    Coef <- t(do.call(rbind,lapply(1:length(Est), function(j) Est[[j]]$coef)))
    Sd   <- t(do.call(rbind,lapply(1:length(Est), function(j) Est[[j]]$sd)))


    indx_lst <- gen_wavelet_indx(log2(dim(Y_f)[2]))
    if (is.null(nrow(Coef))){

      coef_p <- do.call(c, lapply(1:length(indx_lst),function(j)
                                   ashr::ash(Coef[ indx_lst[[j]]],
                                             Sd  [ indx_lst[[j]]]
                                          )$result$PosteriorMean
                                  )
                          )
      temp <- wavethresh::wd(rep(0, length(Coef)[2]))

      temp$D                     <-  coef_p[ -indx_lst[[length(indx_lst)]]]
      temp$C[length(temp$C)]     <-  coef_p[ indx_lst[[length(indx_lst)]]]
      fitted_coef   <-  wavethresh::wr(temp)
    }else{
      fitted_coef <-list()
      for ( i in 1:nrow(Coef ))
      {
        #coef_p <- ashr::ash(Coef[i,],Sd[i,])$result$PosteriorMean
        coef_p <- do.call(c, lapply(1:length(indx_lst),function(j)
                                    ashr::ash(Coef[i,indx_lst[[j]]],
                                              Sd  [i,indx_lst[[j]]]
                                              )$result$PosteriorMean
                                     )
                          )

        temp <- wavethresh::wd(rep(0, dim(Coef)[2]))

        temp$D                     <-  coef_p[  -indx_lst[[length(indx_lst)]]]
        temp$C[length(temp$C)]     <-   coef_p[  indx_lst[[length(indx_lst)]]]
        fitted_coef[[i]]   <-  wavethresh::wr(temp)
      }
      fitted_coef <- do.call(cbind,   fitted_coef )
    }


    Y_adjusted  = Y-cbind(rep(1,nrow(X)),X)%*% t(fitted_coef)
  }

  out <- list( Y_adjusted  = Y_adjusted,

               fitted_coef = fitted_coef,
               pos         = outing_grid
  )


  return(out)
}


