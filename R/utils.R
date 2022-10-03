

which_lowcount <- function( Y_f, thresh_lowcount ){
  tt <- which( apply( abs(Y_f),2,median) <thresh_lowcount )

  if(length(tt)==0)
  {
    return(NULL)
  }else{
    return(tt)
  }
}

#testing if x is a wholenumber
#'
is.wholenumber <- function (x, tol = .Machine$double.eps^0.5)
  abs(x - round(x)) < tol

# Not in operator

'%!in%' <- function(x,y)!('%in%'(x,y))

#based on Rfast implementation
#'
fast_lm <- function(x,y)
{
  be <- solve(crossprod(x),crossprod(x,y))
  #resid <-  y - x %*% be
  #out <- list(be = be,
  #            residuals = resid)

  #return(out)
  return(be)
}


#Circular permutation on vector
# Code adapted from https://mzuer.github.io
#'
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
