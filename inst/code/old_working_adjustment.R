temp <- as.data.frame(cbind(PCA_t, meta_t))
if( table(temp$Batch)[4] < 2)
{
  temp$Batch[ which( temp$Batch==4)] <- 1
}
ff <- ~-1+PCA1+Age+factor(Batch)  #ensure that each Batch effect is estimated individually
utils::str(m <- model.frame(ff, temp))
mat <- model.matrix(ff, m)

Cov <-   mat
adjust_data <- adjust_FM_covariate(Y=Y_t,X=Cov )



adjust_FM_covariate <- function (Y, X, pos = NULL,Batch=NULL ,thresh_lowcount = 0, verbose)
{
  if (is.null(pos)) {
    pos <- 1:dim(Y)[2]
  }
  if (!is.wholenumber(log2(dim(Y)[2])) | !(sum(duplicated(diff(pos))) ==
                                           (length(pos) - 2))) {
    inter_pol.obj <- interpol_mat(Y, pos)
    Y <- inter_pol.obj$Y
    bp <- inter_pol.obj$bp
    start_pos <- min(pos)
    end_pos <- max(pos)
    outing_grid <- start_pos + (end_pos - start_pos)/(length(pos)) *
      inter_pol.obj$grid
    if (verbose) {
      message("Response matrix dimensions not equal to nx 2^J \n or unevenly spaced data \n interpolation procedure used")
    }
  }
  else {
    outing_grid <- pos
  }

  W   <- DWT2(Y)
  Y_f <- cbind(W$D, W$C)
  v1  <- rep(1, dim(X)[2])
  lowc_wc  <- which_lowcount(Y_f, thresh_lowcount)
  coef_Y_f <- cal_Bhat_Shat(Y_f, X, v1)
  indx_lst <- gen_wavelet_indx(log2(dim(Y_f)[2]))
  if (is.null(nrow(coef_Y_f$Bhat))) {
    coef_p <- ashr::ash(coef_Y_f$Bhat, coef_Y_f$Shat)$result$PosteriorMean
    temp <- wavethresh::wd(rep(0, length(coef_Y_f$Bhat)[2]))
    temp$D <- coef_p[-indx_lst[[length(indx_lst)]]]
    temp$C[length(temp$C)] <- coef_p[indx_lst[[length(indx_lst)]]]
    fitted_coef <- wavethresh::wr(temp)
  }
  else {
    fitted_coef <- list()
    for (i in 1:nrow(coef_Y_f$Bhat)) {
      coef_p <- ashr::ash(coef_Y_f$Bhat[i, ], coef_Y_f$Shat[i,])$result$PosteriorMean
      temp <- wavethresh::wd(rep(0, dim(coef_Y_f$Bhat)[2]))
      temp$D <- coef_p[-indx_lst[[length(indx_lst)]]]
      temp$C[length(temp$C)] <- coef_p[indx_lst[[length(indx_lst)]]]
      fitted_coef[[i]] <- wavethresh::wr(temp)
    }
    fitted_coef <- do.call(cbind, fitted_coef)
  }
  Y_adjusted = Y - X %*% t(fitted_coef)
  out <- list(Y_adjusted = Y_adjusted, fitted_coef = fitted_coef,
              pos = outing_grid)
  return(out)
}
