## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.height = 3,
  fig.align = "center",
  fig.cap = "&nbsp;",
  ddpi = 175
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
set.seed(1)
data(N3finemapping)
attach(N3finemapping)

rsnr <- 0.5 #expected root signal noise ratio

pos1 <- 25    #Position of the causal covariate for effect 1
pos2 <- 75    #Position of the causal covariate for effect 2
lev_res <- 7#length of the molecular phenotype (2^lev_res)
f1 <-  simu_IBSS_per_level(lev_res )$sim_func#first effect
f2 <- simu_IBSS_per_level(lev_res )$sim_func #second effect

plot( f1, type ="l", ylab="effect", col="blue")
abline(a=0,b=0)
lines(f2, type="l", col="green")
legend(x=100,
       y=3,
       lty = rep(1,3),
      legend= c("effect 1", "effect 2" ),
       col=c( "blue","green"))


## -----------------------------------------------------------------------------

noisy.data  <- list()
X <- N3finemapping$X[,1:100]
for ( i in 1:nrow(X))
{
  f1_obs <- f1
  f2_obs <- f2
  noise <- rnorm(length(f1), sd=  (1/  rsnr ) * var(f1))
  noisy.data [[i]] <-  X[i,pos1]*f1_obs +X[i,pos2]*f2_obs + noise

}
noisy.data <- do.call(rbind, noisy.data)

Y <- noisy.data



## -----------------------------------------------------------------------------
 out <- susiF(Y,X,L=10) 


## -----------------------------------------------------------------------------
out$cs

## -----------------------------------------------------------------------------
 
 plot_susiF(out )
 


## -----------------------------------------------------------------------------
 plot( f1, type="l", main="Estimated effect 1", xlab="")
 lines(get_fitted_effect(out,l=2),col='blue' )
 abline(a=0,b=0)
 legend(x= 35,
        y=3,
        lty= rep(1,2),
        legend = c("effect 1"," fSuSiE est "),
        col=c("black","blue" )
 )
 plot( f2, type="l", main="Estimated effect 2", xlab="")
 lines(get_fitted_effect(out,l=1),col='green' )
 abline(a=0,b=0)
 legend(x= 20,
        y=-1.5,
        lty= rep(1,2),
        legend = c("effect 2"," fSuSiE est "),
        col=c("black","green" )
 )
  

## -----------------------------------------------------------------------------
 
out$pip[1:10]

plot_susiF(out, pip_only=TRUE)

## -----------------------------------------------------------------------------
 affected_reg( out)
  

## -----------------------------------------------------------------------------
true_sig  <- matrix( X[ ,pos1], ncol=1)%*%t(f1_obs) +matrix( X[ ,pos2], ncol=1)%*%t(f2_obs)             


plot(out$ind_fitted_func,  true_sig,
     xlab= "predicted value", ylab="true value") 
abline(a=0,b=1)

## -----------------------------------------------------------------------------

plot(out$ind_fitted_func,  Y,
     xlab= "predicted value", ylab="Noisy observation value") 
abline(a=0,b=1)

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
set.seed(2)
data(N3finemapping)
attach(N3finemapping)

rsnr <- 0.5 #expected root signal noise ratio

pos1 <- 25    #Position of the causal covariate for effect 1
pos2 <- 75    #Position of the causal covariate for effect 2
lev_res <- 9#length of the molecular phenotype (2^lev_res)
f1 <-  simu_IBSS_per_level(lev_res )$sim_func  #first effect
f2 <- simu_IBSS_per_level(lev_res )$sim_func  #second effect


## Suppose we observe these functions on 200 points
sampling_pos <- sample( 1:length(f1), size=200)
sampling_pos <-sampling_pos  [order(sampling_pos)]

plot(f2, type="l", main="underlying function and the different sampling position")
 points(sampling_pos, f2[sampling_pos], col="red",pch=20)
 points(sampling_pos, rep(0, length(sampling_pos)), col="red",pch=3)
 abline(0,0)
 

## ----message=FALSE, warning=FALSE---------------------------------------------
  f1_obs <- f1[sampling_pos]
  f2_obs <- f2[sampling_pos]
noisy.data  <- list()
X <- N3finemapping$X[,1:100]
for ( i in 1:nrow(X))
{

  noise <- rnorm(length(f1_obs ), sd=  2)
  noisy.data [[i]] <-  X[i,pos1]*f1_obs +X[i,pos2]*f2_obs + noise

}
noisy.data <- do.call(rbind, noisy.data)

Y <- noisy.data


## ----message=FALSE, warning=FALSE---------------------------------------------
 out <- susiF(Y,X,L=3 ,pos= sampling_pos)

## -----------------------------------------------------------------------------
plot(   x=1:512,           y=f2,type="l")
points( x=sampling_pos,    y=f2_obs, col="red",pch=20)
lines(  x=out$outing_grid, y=out$fitted_func[[1]]    ,col='darkblue', lwd=2)
lines(  x=out$outing_grid, y= out$cred_band[[1]][1,] ,col='darkblue',lty=2 )
lines(  x=out$outing_grid, y= out$cred_band[[1]][2,] ,col='darkblue' ,lty=2 )

## -----------------------------------------------------------------------------
 out1 <- susiF(Y,X,L=3 , prior = 'mixture_normal_per_scale',verbose=FALSE)
 out1$runtime
 out2 <- susiF(Y,X,L=3 , prior = 'mixture_normal',verbose=FALSE)
 out2$runtime

