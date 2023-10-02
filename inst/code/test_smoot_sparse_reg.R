set.seed(2)
library(wavethresh)
library(ashr)
library(smashr)
source("D:/Document/Serieux/Travail/Package/susiF.alpha/inst/code/fit_hmm.R" )
nois_lev <- 3
f1 <-   DJ.EX( )$blocks
f2 <-   DJ.EX( )$doppler

f2[200:300] <- 3
f2[600: 800] <- 0



f3 <-   DJ.EX( )$heavi

f3[100:150]<- 0
f3[800:900]<- 3

f4 <-  DJ.EX( )$bumps



f1_noisy <-  f1 +rnorm(f1,sd=nois_lev)


tt1 <- fit_hmm(x=f1_noisy,
               sd=rep(nois_lev, length(f1)),
               mult=2,
               halfK= 50,maxiter=10,
               smooth = FALSE)


tt2 <-fit_hmm(x=f1_noisy,
              sd=rep(nois_lev, length(f1)),
              mult=2,
              halfK= 50,maxiter=10 )


plot( f1, type="l" ,lwd=2   )

lines(smash(f1_noisy), col="blue")
lines( (tt1$x_post ) , lwd=2, col="red")
lines( (tt2$x_post   ) , lwd=2, col="green")


f2_noisy <-  f2 +rnorm(f1,sd=nois_lev)


tt1 <- fit_hmm(x=f2_noisy,
               sd=rep(nois_lev, length(f1)),
               mult=2,
               halfK= 20,maxiter=10,
               smooth = FALSE)


tt2 <-fit_hmm(x=f2_noisy,
              sd=rep(nois_lev, length(f1)),
              mult=2,
              halfK= 20,maxiter=10 )




plot( f2, type="l" ,lwd=2   )
lines(smash(f2_noisy), col="blue")
 lines( (tt1$x_post ) , lwd=2, col="red")
 lines( (tt2$x_post   ) , lwd=2, col="green")






f3_noisy <-  f3 +rnorm(f1,sd=nois_lev)


tt1 <- fit_hmm(x=f3_noisy,
               sd=rep(nois_lev, length(f1)),
               mult=2,
               halfK= 50,maxiter=10,
               smooth = FALSE)


tt2 <-fit_hmm(x=f3_noisy,
              sd=rep(nois_lev, length(f1)),
              mult=2,
              halfK= 50,maxiter=10 )


plot( f3, type="l" ,lwd=2   )
 lines(smash(f3_noisy), col="blue")


lines( (tt1$x_post ) , lwd=2, col="red")


lines( (tt2$x_post  ) , lwd=2, col="green")


f4_noisy <-  f4 +rnorm(f1,sd=nois_lev)


tt1 <- fit_hmm(x=f4_noisy,sd=rep(nois_lev, length(f1)),
               halfK= 200,mult=2 ,smooth = FALSE,maxiter=4)
tt2 <-fit_hmm(f4_noisy  ,
              sd=rep(nois_lev, length(f1)),
              mult=2,
              halfK= 200 ,maxiter=4)

plot( f4, type="l" ,lwd=2   )
lines(smash(f4_noisy), col="blue")


lines( (tt1$x_post ) , lwd=2, col="red")

lines( ( tt2$x_post  ) , lwd=2, col="green")


data <- set_data(0,0.2)
calc_loglik(beta.ash,data)
beta.ash$fitted_g
