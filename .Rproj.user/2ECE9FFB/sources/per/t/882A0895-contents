library(susiF.alpha)
library(ashr)
library(microbenchmark)
        x0 <- rnorm(200)
        x11 <- rnorm(199)
        x1 <- c(NA, x11)
   x21 <-      rnorm(100)
x2 <- c(rep(NA, 100), x21)
x31 <- rnorm(1)
x3 <- c(rep(NA, 199),x31)

f1 <- function(x, idx)
{
   x <- x[idx]
   sapply(1:100, function(i) crossprod(x))
}





idx0 <- which(!is.na(x0))
idx1 <- which(!is.na(x1))
idx2 <- which(!is.na(x2))
idx3 <- which(!is.na(x3))

microbenchmark(times=100000L,

  v1 =sapply(1:100, function(i) crossprod(x11)),
  v01= f1(x1 ,idx1 ),
  v2= sapply(1:100, function(i) crossprod(x21)),
  v02= f1(x2 ,idx2 ),
  v3= sapply(1:100, function(i) crossprod(x11)),
  v03= f1(x3 ,idx3 )
)
