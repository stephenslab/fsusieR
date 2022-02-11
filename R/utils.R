
#testing if x is a wholenumber
#'@export
is.wholenumber <- function (x, tol = .Machine$double.eps^0.5)
  abs(x - round(x)) < tol

#' @title  Not in operator
#'@export
'%!in%' <- function(x,y)!('%in%'(x,y))

#based on Rfast implementation
#'@export
fast_lm <- function(x,y)
{

  be <- solve(crossprod(x), crossprod(x, y))
  resid <-  y- x%*%xbe
  out <- list( be=be,
               residuals=resid)
  return(out)
}
