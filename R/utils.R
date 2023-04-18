
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
  cm = colMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = colSds(x, center = cm)
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
