## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.height = 3,
  fig.align = "center",
  fig.cap = "&nbsp;",
  dpi = 175
  )

## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	message = FALSE,
	warning = FALSE
)
library(fsusieR)
library(susieR)
library(wavethresh)
set.seed(2)
data(N3finemapping)
attach(N3finemapping)

rsnr <- 1 #expected root signal noise ratio

pos1 <- 25   #Position of the causal covariate for effect 1
pos2 <- 75   #Position of the causal covariate for effect 2
lev_res <- 7#length of the molecular phenotype (2^lev_res)
f1 <-  simu_IBSS_per_level(lev_res )$sim_func#first effect
f2 <- simu_IBSS_per_level(lev_res )$sim_func #second effect

f1_cov  <- simu_IBSS_per_level(lev_res )$sim_func #effect cov 1
f2_cov  <- simu_IBSS_per_level(lev_res )$sim_func #effect cov 2
f3_cov  <- simu_IBSS_per_level(lev_res )$sim_func #effect cov 3

 

## -----------------------------------------------------------------------------

Geno <- N3finemapping$X[,1:100]
tt <- svd(N3finemapping$X[,500:700])

Cov <- matrix(rnorm(3*nrow(Geno ),sd=2), ncol=3)
target.data  <-list()
noisy.data  <-list()
for ( i in 1:nrow(Geno))
{
  f1_obs <- f1
  f2_obs <- f2
  noise <- rnorm(length(f1), sd=  (1/  rsnr ) * var(f1))
  target.data [[i]] <-  Geno [i,pos1]*f1_obs +noise +Geno [i,pos2]*f2_obs 
  noisy.data [[ i]] <-   Cov [i,1]*f1_cov  + Cov  [i,2]*f2_cov + Cov  [i,3]*f3_cov

}
technical.noise <- do.call(rbind, noisy.data)
target.data <- do.call(rbind, target.data)
Y <- technical.noise+target.data


## -----------------------------------------------------------------------------

Est_effect <- EBmvFR(Y,X=Cov,adjust=TRUE )
plot(Est_effect$fitted_func[1,],type = "l",
     col="blue",
     lty=2,
     lwd=2)
lines(f1_cov)
legend(x=60,y=-1,
       lwd=c(1,2),
       lty=c(1,2),
       legend=c('effect','estimated'))



Y_corrected  <-Est_effect$Y_adjusted

## -----------------------------------------------------------------------------
  plot(  Y_corrected, Y-Cov %*%Est_effect$fitted_func )

## -----------------------------------------------------------------------------

par(mfrow=c(1,2))
plot( Y , target.data)
plot( Y_corrected , target.data)
par(mfrow=c(1,1))

## ----echo=FALSE---------------------------------------------------------------
out <- susiF(Y=Y_corrected  ,
             X=Geno,
             L=10)

## -----------------------------------------------------------------------------
 
 plot_susiF(out)
 


## ----echo=FALSE---------------------------------------------------------------
out_unadj <- susiF(Y=Y,
             X=Geno,
             L=10 )

## -----------------------------------------------------------------------------
 
 plot( f1, type="l", main="Estimated effect 1", xlab="")
 lines(get_fitted_effect(out,l=2),col='blue' )
 lines(get_fitted_effect(out_unadj,l=2),col='blue' , lty=2)
 
 abline(a=0,b=0)
 legend(x= 20,
        y=-0.5,
        lty= c(1,1,2),
        legend = c("effect 1","fSuSiE ajd", "fSuSiE unajd " ),
        col=c("black","blue","blue" )
 )
 plot( f2, type="l", main="Estimated effect 2", xlab="")
 lines(get_fitted_effect(out,l=1),col='darkgreen' )
 lines(get_fitted_effect(out_unadj,l=1),col='darkgreen' , lty=2)
 
 abline(a=0,b=0)
 legend(x= 10,
        y=2.5,
        lty= c(1,1,2),
        legend = c("effect 1","fSuSiE ajd", "fSuSiE unajd " ),
        col=c("black","darkgreen","darkgreen" )
 )
  

