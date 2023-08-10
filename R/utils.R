
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
fast_lm <- function(x,y)
{

    be <- solve(crossprod(x),crossprod(x,y))
    sd <-  sqrt(fast_var(y - x %*% be) /(length(x)-1))


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

fast_var <- function (x)
{
  .Call(stats:::C_cov, x, x, 5, FALSE)
}


#' @importFrom matrixStats colSds
#from https://www.r-bloggers.com/2016/02/a-faster-scale-function/
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
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = t( (t(x) - cm) / csd )
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
    n <- nrow(x)
    d = n*cm^2 + (n-1)*csd^2
    d = (d - n*cm^2)/csd^2
    attr(x, "d") <- d
  }
  return(x)
}


gen_EM_out <- function(tpi_k , lBF){

  out <- list(tpi_k = tpi_k,lBF = lBF)
  class(out) <- c("EM_pi","list")
  return(out)

}

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
#' @param susiF.obj at fitted susiF.obj object
#' @export
affected_reg <- function( susiF.obj){
  outing_grid <- susiF.obj$outing_grid

  reg <-  list()
  h <- 1
  for (   l in 1:length(susiF.obj$cs)){

    pos_up <-  which(susiF.obj$cred_band[[l]][1,]<0)
    pos_low <- which(susiF.obj$cred_band[[l]][2,]>0)


    reg_up <- split( pos_up,cumsum(c(1,diff( pos_up)!=1)))

    reg_low <- split( pos_low,cumsum(c(1,diff( pos_low)!=1)))
    for( k in 1:length(reg_up)){
      reg[[h]] <- c(l, outing_grid[reg_up[[k]][1]], outing_grid[reg_up[[k]][length(reg_up[[k]])]])

      h <- h+1
    }
    for( k in 1:length(reg_low )){
      reg[[h]] <- c(l, outing_grid[reg_low [[k]][1]], outing_grid[reg_low [[k]][length(reg_low [[k]])]])

      h <- h+1
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

effective.effect=function(betahat,se,df){

  p = 2 * pt(abs(betahat/se  ), df=n ,
             lower.tail = FALSE)

  sign(betahat)*stats::qnorm(p / 2, sd=se ,lower.tail=FALSE)


}




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
