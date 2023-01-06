################################## WAVELET UTILITY FUNCTIONS ############################

# @title Lifting scheme for non decimated and unevenly spaced data on matrix
#
# @description  Interpolation procedure from Kovac and Silvermann 2000 for a matrix of functions
#
#
# @details We suppose that we observe n functions/curves at T different time points. We remap these functions into a grid of length equal to a power of two using the lifting scheme of Kovac and Silvermann 2000.
#
# @param Y matrix of curves of size nXT, where the curves are stored by row
#
# @param pos vector of positive values, each the vector component corresponds to the measurement time. (the values of pos should increase along with the index)
#
# @importFrom  wavethresh makegrid
#
# @export
#
interpol_mat <- function(Y, pos)
{

  bp    <- (pos- min(pos))/(max(pos)-min(pos))
  Y_new <- t(apply(Y, 1, interpolKS, bp=bp))
  grid  <- wavethresh::makegrid(t=bp, y = 1:dim(Y)[2], gridn = 2^(floor(log(length(pos)-1,2)) + 1)   )$gridy

  out <- list(Y    = Y_new,
              grid = grid,
              bp   = bp)
  return(out)
}



# @title Interpolation procedure from Kovac and Silverman 2000
#
# @description  Interpolation procedure from Kovac and Silveramnn 2000
#
# @param y vector of observed curves
#
# @param bp scaled observation time points on 0,1
#
# @importFrom wavethresh makegrid
#
# @export
interpolKS <-  function (y, bp)
{
  wavethresh::makegrid(t=bp,y=y, gridn = 2^(floor(log(length(bp)-1,2)) + 1)   )$gridy
}


# @title  wavelet transform on a matrix of functions
#
# @description function adapted from grove R package from Ma and Soriano. Perform wavelet transform of each row of a matrix
#
# @param data matrix of size NxJ where J is a power of two
#
# @param filter.number This selects the smoothness of the wavelet that you want to use in the decomposition. By default, this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments. Description from Guy Nason (wavethresh package)
#
# @param family specifies the family of wavelets that you want to use. Two popular options are "DaubExPhase" and "DaubLeAsymm" but see the help for filter.select for more possibilities.. Description from Guy Nason (wavethresh package)
#
# @return A list with the following components
#
#\item{C}{ vector length n containing the wavelet C coefficient}
#\item{D}{ matrix of size nx 2^(J -1) where each row contains the wavelet D coefficients, ordered in the same way as in the wavethresh package}
#\item{family}{ used for the wavelet transform}
#\item{filter.number}{ }
#
# @importFrom wavethresh accessC
# @importFrom wavethresh wd
# @export
DWT2 <- function (data, filter.number = 10, family = "DaubLeAsymm")
{

  J <- ncol(data)
  n <- nrow(data)
  D <- matrix(NA, nrow = n, ncol = J - 1)
  C <- rep(NA, n)
  for (i in 1:n) { ## Speed Gain
    temp <- wavethresh::wd(data[i, ], filter.number = filter.number,
               family = family)
    D[i, ] <- temp$D
    C[i] <- wavethresh::accessC(temp, level = 0)
  }
  output <- list(C = C, D = D, J = log2(J), filter.number = filter.number,
                 family = family)
  class(output) <- "DWT"
  return(output)
}


#' @title Generate list of index wavelet coefficients
#
#' @description Generate list of index corresponding to where the wavelet coefficient of scale s are stored in wd$D.
#First element of the list corresponds to the index for the wavelet coefficient at scale 0, the second component corresponds to the index for the wavelet coefficient at scale 1 (and so on). The last componnet of the list correspond to where the C coefficients is stored in init_prior for class mixture_normal_per_scale
#
#' @param lev_res integer, corresponding to log2 of the signal length. WARNING the ordering change for different values of s.
#
#' @export
#
#' @examples
#' library(wavethresh)
#' tem_func <- rnorm( 2^8)
#' twav <- wd(tem_func)
#' indx_lst <- gen_wavelet_indx(8)
#' plot(accessD(twav,level=6), ( twav$D[unlist(indx_lst[(6+1)])]) )
#'  #should a straightline
gen_wavelet_indx <- function(lev_res)
{
  indx_lst <- list()
  indx_lst[[1]] <- 2^lev_res -1 #coefficient
  for ( s in 1:(lev_res-1))
  {

    indx  <- 2^(lev_res)- (2^((s+1))-1) :(2^s)


    indx_lst[[s+1]] <- indx
  }
  indx_lst[[length(indx_lst)+1]] <- 2^lev_res# C coefficient
  out <- indx_lst
  return(out)

}




#' @title Adjusting functions for covariates
#
#' @description Simple function to adjust the observed curves for covariate (e.g. age, population structure)
#' @param Y  observed curves.
#' @param X matrix of containning the covaraite
#' @param pos Original position of the column of Y
#' @param thresh_lowcount numeric, use to check the wavelet coefficients have
#'  problematic distribution (very low dispersion even after standardization).
#'  Basically check if the median of the absolute value of the distribution of
#'   a wavelet coefficient is below this threshold, if yes the algorithm discard
#'   this wavelet coefficient (setting its estimate effect to 0 and estimate sd to 1).
#'   Set to 0 by default. Can be useful when analyzing sparse data from sequence
#'    based assay or small samples.
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings, are printed to the
#' console.
#' @export
#


adjust_FM_covariate <- function(Y,X,pos=NULL, thresh_lowcount=0, verbose){


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


  v1 <- rep(1, dim(X)[2])
  lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)

  fitted_Y_f <- do.call( cbind,lapply(1:dim(Y_f)[2],
                                      function(j) if(j%in%lowc_wc){
                                        return(rep(0,nrow(Y)))
                                      }else{
                                        return(fitted(lm(Y[,j]~X)))
                                      }


  )

  )


  Y_fitted <- 0*Y
  indx_lst <- gen_wavelet_indx(lev_res = log2(ncol(fitted_Y_f)))
  temp <- wavethresh::wd(rep(0, dim(fitted_Y_f)[2]))
  for ( i in 1:nrow(Y))
  {
    temp$D                     <-  fitted_Y_f[i,-indx_lst[[length(indx_lst)]]]
    temp$C[length(temp$C)]     <-   fitted_Y_f[i,indx_lst[[length(indx_lst)]]]
    Y_fitted [i,] <-  wavethresh::wr(temp)
  }


  Y_adjusted = Y - Y_fitted

  out <- list( Y_adjusted  = Y_adjusted,
               Y_fitted    = Y_fitted,
               pos         = outing_grid
  )

  return(out)
}


