
#testing if x is a wholenumber
#'@export
is.wholenumber <- function (x, tol = .Machine$double.eps^0.5)
  abs(x - round(x)) < tol

# Not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))

#based on Rfast implementation
#'@export
fast_lm <- function(x,y)
{
  be <- solve(crossprod(x),crossprod(x,y))
  resid <-  y - x %*% be
  out <- list(be = be,
              residuals = resid)
  return(out)
}
