################################## WAVELET UTILITY FUNCTIONS ############################

#' @title Lifting scheme for non decimated and unevenly spaced data on matrix
#'
#' @description  Interpolation procedure from Kovac and Silvermann 2000 for a matrix of functions
#'
#'
#' @details We suppose that we observe n functions/curves at T different time points. We remap these functions into a grid of length equal to a power of two using the lifting scheme of Kovac and Silvermann 2000.
#'
#' @param Y matrix of curves of size nXT, where the curves are stored by row
#'
#' @param pos vector of positive values, each the vector component corresponds to the measurement time. (the values of pos should increase along with the index)
#'
#' @importFrom  wavethresh makegrid
#'
#' @export
#'
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



#' @title Interpolation procedure from Kovac and Silverman 2000
#'
#' @description  Interpolation procedure from Kovac and Silveramnn 2000
#'
#' @param y vector of observed curves
#'
#' @param bp scaled observation time points on 0,1
#'
#' @importFrom wavethresh makegrid
#'
#' @export
interpolKS <-  function (y, bp)
{
  wavethresh::makegrid(t=bp,y=y, gridn = 2^(floor(log(length(bp)-1,2)) + 1)   )$gridy
}


#' @title  wavelet transform on a matrix of functions
#'
#' @description function adapted from grove R package from Ma and Soriano. Perform wavelet transform of each row of a matrix
#'
#' @param data matrix of size NxJ where J is a power of two
#'
#' @param filter.number This selects the smoothness of the wavelet that you want to use in the decomposition. By default, this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments. Description from Guy Nason (wavethresh package)
#'
#' @param family specifies the family of wavelets that you want to use. Two popular options are "DaubExPhase" and "DaubLeAsymm" but see the help for filter.select for more possibilities.. Description from Guy Nason (wavethresh package)
#'
#' @return A list with the following components
#'
#'\item{C}{ vector length n containing the wavelet C coefficient}
#'\item{D}{ matrix of size nx 2^(J -1) where each row contains the wavelet D coefficients, ordered in the same way as in the wavethresh package}
#'\item{family}{ used for the wavelet transform}
#'\item{filter.number}{ }
#'
#' @importFrom wavethresh accessC
#' @importFrom wavethresh wd
#' @export
DWT2 <- function (data, filter.number = 10, family = "DaubLeAsymm")
{

  J <- ncol(data)
  n <- nrow(data)
  D <- matrix(NA, nrow = n, ncol = J - 1)
  C <- rep(NA, n)
  for (i in 1:n) { ## Speed Gain
    temp <- wd(data[i, ], filter.number = filter.number,
               family = family)
    D[i, ] <- temp$D
    C[i] <- accessC(temp, level = 0)
  }
  output <- list(C = C, D = D, J = log2(J), filter.number = filter.number,
                 family = family)
  class(output) <- "DWT"
  return(output)
}


#' @title Generate list of index wavelet coefficients
#'
#' @description Generate list of index corresponding to where the wavelet coefficient of scale s are stored in wd$D.
#'First element of the list corresponds to the index for the wavelet coefficient at scale 0, the second component corresponds to the index for the wavelet coefficient at scale 1 (and so on). The last componnet of the list correspond to where the C coefficients is stored in init_prior for class mixture_normal_per_scale
#'
#' @param s integer, corresponding to log2 of the signal length. WARNING the ordering change for different values of s.
#'
#' @export
#'
#' @examples
#' tem_func <- rnorm( 2^8)
#' twav <- wd(tem_func)
#' indx_lst <- gen_wavelet_indx(8)
#' plot(accessD(twav,level=6), ( twav$D[unlist(indx_lst[(6+1)])]) )
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


