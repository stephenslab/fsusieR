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
