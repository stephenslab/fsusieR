
#'@title Find wavelet with critically low variance
#' @description Find wavelet with critically low variance
#' @param Y_f a matrix of wavelet coefficients (col= wc, row =ind)
#' @param thresh_lowcount mad
#' @importFrom stats median
#' @export
#' @keywords internal
which_lowcount <- function( Y_f, thresh_lowcount ){
  tt <- which( apply( abs(Y_f),2,median) <= thresh_lowcount )

  if(length(tt)==0)
  {
    return(NULL)
  }else{
    return(tt)
  }
}

# Testing if x is a whole number.
is.wholenumber <- function (x, tol = .Machine$double.eps^0.5)
  abs(x - round(x)) < tol

# Not in operator

'%!in%' <- function(x,y)!('%in%'(x,y))

# Based on Rfast implementation.
#' @importFrom Rfast cova
fast_lm <- function(x,y)
{

  be <- solve(crossprod(x),crossprod(x,y))
  sd <-  sqrt(Rfast::cova(y - x %*% be) /(length(x)-1))


  return(c(be,sd))
}






#Circular permutation on vector
# Code adapted from https://mzuer.github.io
#' @importFrom utils head
#' @importFrom utils tail
shifter <- function(x, n = 1) {
  # if (n == 0) x else c(tail(x, -n), head(x, n))
  if (n == 0) x else c(tail(x, n), head(x, -n))
}

#shifter(c(1:10), n=-1)
# [1]  1  2  3  4  5  6  7  8  9 10
#shifter(c(1:10), n=1)
# [1] 10  1  2  3  4  5  6  7  8  9
#shifter(c(1:10), n=2)
# [1]  9 10  1  2  3  4  5  6  7  8
#
#' @importFrom stats qqnorm
Quantile_transform  <- function(x)
{

  x.rank = rank(x, ties.method="random")
  #x.rank = rank(x, ties.method="average")
  return(qqnorm(x.rank,plot.it = F)$x)
}


#' @title  Scaling function from r-blogger
#'
#' @description from https://www.r-bloggers.com/2016/02/a-faster-scale-function/
#' @importFrom matrixStats colSds
#
#' @param x a matrix
#' @param center logical if true center column
#' @param scale logical if true scale column
#' @param add_attr logical  https://www.r-bloggers.com/2016/02/a-faster-scale-function/
#' @param rows logical  https://www.r-bloggers.com/2016/02/a-faster-scale-function/
#' @param  cols logical  https://www.r-bloggers.com/2016/02/a-faster-scale-function/
#' @export


colScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {

  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }

  ################
  # Get the column means
  ################
  cm =  colMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = matrixStats::colSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }

  # Handle zero variance columns
  zero_var_cols = which(csd == 0)
  for (col in zero_var_cols) {

    csd[col] <- 1  # Set the scaled:scale attribute for this column to 1
  }

  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = t((t(x) - cm) / csd)

  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
    n <- nrow(x)
    d = n * cm^2 + (n - 1) * csd^2
    d = (d - n * cm^2) / csd^2
    attr(x, "d") <- d
  }

  return(x)
}



gen_EM_out <- function(tpi_k , lBF){

  out <- list(tpi_k = tpi_k,lBF = lBF)
  class(out) <- c("EM_pi","list")
  return(out)

}

#' @title  Compute purity of a list of CS
#'
#' @description Compute purity
#' @param l_cs the list of Credible sets
#' @param X the regression matrix
#' @export
#' @keywords internal
cal_purity <- function(l_cs,X){
  tt <- list()
  for (k in 1:length(l_cs)){
    if(length(unlist(l_cs[[k]]))==1 ){
      tt[[k]] <- 1
    }else{
      x <-abs( cor(X[,unlist(l_cs[[k]]   ) ]))


      tt[[k]] <-  min( x[col(x) != row(x)])
    }
  }
  return( tt )
}



#' @title Extract coordinates of the regions affected by the different CS
#'
#' @description Extract coordinates of the regions affected by the different CS
#'
#' @details return a matrix with 3 columns. Each lines corresponds to a regions
#' in which the estimated credible bands are "crossing zero"/i.e. the effects are likely not to be 0 in this region.
#' the second column corresponds to the start of the region  and the third to the end of the affected region
#'
#' @param obj at fitted obj object
#' @param lfsr_thresh threshold for affected region when using HMM postprocessing
#' @importFrom stats complete.cases
#'
#' @export
#'
affected_reg <- function( obj, lfsr_thresh=0.05){
  outing_grid <- obj$outing_grid

  reg <-  list()
  h <- 1
  if(!is.null (obj$cred_band)){

    for (   l in 1:length(obj$cs)){

      pos_up <-  which(obj$cred_band[[l]][1,]<0)
      pos_low <- which(obj$cred_band[[l]][2,]>0)


      reg_up <- split( pos_up,cumsum(c(1,diff( pos_up)!=1)))

      reg_low <- split( pos_low,cumsum(c(1,diff( pos_low)!=1)))
      if( length(reg_up[[1]]) >0){
        for( k in 1:length(reg_up)){
          reg[[h]] <- c(l, outing_grid[reg_up[[k]][1]], outing_grid[reg_up[[k]][length(reg_up[[k]])]])

          h <- h+1
        }
      }

      if( length(reg_low[[1]]) >0){
        for( k in 1:length(reg_low )){
          reg[[h]] <- c(l, outing_grid[reg_low [[k]][1]], outing_grid[reg_low [[k]][length(reg_low [[k]])]])

          h <- h+1
        }
      }



    }
  }

  if(!is.null(obj$lfsr_func)){

    for (   l in 1:length(obj$cs)){

      pos_up <-  which(obj$lfsr_func[[l]] < lfsr_thresh)
      pos_low <- which(obj$cred_band[[l]][2,]>0)

      reg_up <- split( pos_up,cumsum(c(1,diff( pos_up)!=1)))
      for( k in 1:length(reg_up)){
        reg[[h]] <- c(l, outing_grid[reg_up[[k]][1]], outing_grid[reg_up[[k]][length(reg_up[[k]])]])

        h <- h+1
      }



    }
  }

  reg <-  do.call(rbind, reg)
  colnames(reg) <- c("CS", "Start","End")
  reg <- as.data.frame(reg)
  reg <- reg[order(reg$CS, reg$Start),]
  reg <- reg[complete.cases(reg),]
  return(reg)
}



#From Lu and Stephens
#p is  a log p to avoid underflow
#
#' @importFrom stats pt
#' @importFrom stats qnorm
effective.effect=function(betahat,se,df){
  lp =log(  2)  + pt(abs(betahat/se  ), df=df ,
                     lower.tail = FALSE,
                     log.p=TRUE)

  sign(betahat)*stats::qnorm( lp- log(  2),
                              sd=se ,
                              lower.tail=FALSE,
                              log.p=TRUE )

}
# from the ashr package
pval2se = function(bhat,p){
  z = qnorm(1-p/2)
s = abs(bhat/z)

return(s)}



update_Shat_pois <- function(Shat, indx_lst, lowc_wc){

  for ( k in 1:length(indx_lst)){

    idx <-  indx_lst[[k]]
    if (!is.null(lowc_wc )){
      if( length(which (lowc_wc %in% indx_lst[[k]] ))>0){
        idx <-  indx_lst[[k]][ -which( lowc_wc%in%indx_lst[[k]])]
      }
    }

    est_sd <-  ifelse(mean(Shat[,idx])>0,
                      mean(Shat[, idx]),
                      1e-10 )

    Shat [,indx_lst[[k]]] <- est_sd
  }
  return(Shat)
}

#' @title Binned reponsed data
#'
#' @param  Y a matrix of observed individual function measure in more than 1024 positions
# (i.e., the result of running fsusie)
#' @param base_pair the corresponding position of each column
#' @param n_bins number of bins for the ouput
#' @export

#' @export
bin_Y <- function(Y, base_pair, n_bins = 1024) {
  total_bp    <- ncol(Y)
  bin_size    <- total_bp / n_bins
  binned_data <- matrix(ncol = n_bins, nrow = nrow(Y))
  start_bin   <- rep(NA, n_bins)

  for (i in 1:n_bins) {
    # Calculate exact start and end positions of bins
    start_col <- min(which(base_pair >= (min(base_pair) + (i - 1) * bin_size)))
    end_col   <- max(which(base_pair < (min(base_pair) + i * bin_size)))

    if (!is.na(start_col) & !is.na(end_col) & start_col <= end_col) {
      binned_data[, i] <- rowSums(Y[, start_col:end_col, drop = FALSE], na.rm = TRUE)
    } else {
      binned_data[, i] <- 0  # Prevent missing values
    }

    start_bin[i] <- base_pair[start_col]  # Ensure correct bin placement
  }

  return(list(binned_data = binned_data,
              pos = start_bin,
              bin_size = bin_size))
}




log1pexp = function (x){
  indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), right = TRUE,
                   include.lowest = TRUE)
  kk <- which(indx == 1)
  if (length(kk)) {
    x[kk] <- exp(x[kk])
  }
  kk <- which(indx == 2)
  if (length(kk)) {
    x[kk] <- log1p(exp(x[kk]))
  }
  kk <- which(indx == 3)
  if (length(kk)) {
    x[kk] <- x[kk] + exp(-x[kk])
  }
  return(x)
}
















reflect_vec <- function (x)
{
  n = length(x)
  J = log2(n)
  if ((J%%1) == 0) {
    x = c(x, x[n:1])
    return(list(x = x, idx = 1:n))
  }
  else {
    n.ext = 2^ceiling(J)
    lnum = round((n.ext - n)/2)
    rnum = n.ext - n - lnum
    if (lnum == 0) {
      x.lmir = NULL
    }
    else {
      x.lmir = x[lnum:1]
    }
    if (rnum == 0) {
      x.rmir = NULL
    }
    else {
      x.rmir = x[n:(n - rnum + 1)]
    }
    x.ini = c(x.lmir, x, x.rmir)
    x.mir = x.ini[n.ext:1]
    x = c(x.ini, x.mir)
    return(list(x = x, idx = (lnum + 1):(lnum + n)))
  }
}




#' @title Perform Haar-Fisz tranform on  count matrix
#
#' @param   count.data a N by T matrix of count data where each row is an observed Poisson processt
#
#
#
#' @return a matrix of size N by T  of the Haar-Fisz transformed data
#
#' @export
#'@export
HFT<- function(count.data){
  lst <- list()
  for ( i in 1:nrow(count.data)){
    lst[[i]] <-haarfisz::hft(count.data[i,])
  }
  out <- do.call( rbind,lst)
  return(out)
}









cs_relation <- function(cs1, cs2) {
  if (setequal(cs1, cs2)) return("identical")
  if (all(cs1 %in% cs2))  return("cs1_in_cs2")
  if (all(cs2 %in% cs1))  return("cs2_in_cs1")
  return("none")
}

