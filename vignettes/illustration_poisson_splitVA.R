## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(fsusieR)
set.seed(1)
n =200

intensities= exp( 3*sin( 1:n/20)+0.1 + rnorm(n, sd=0.5) )
x = rpois(n,intensities)
 
plot(x,  pch=19)
lines(exp(3*sin( 1:n/20)), lwd=3, col="lightgreen")


## -----------------------------------------------------------------------------
 
 tt= ebpm_normal(x)
tt$fitted_g
plot(log(intensities), type="l",lwd=3, col="lightgreen")
points(tt$posterior$mean_log, pch=19) 

## -----------------------------------------------------------------------------
 tt2= pois_mean_GP(x, prior_mean = 3*sin( 1:n/20) ,
                   prior_var = 1)
 
 plot(tt2$posterior$posteriorMean_latent, pch=19, col="blue")
points(tt$posterior$mean_log, pch=19)
 lines(log(intensities), ,lwd=3, col="lightgreen")
 
legend("bottomright", legend = c(  "Split VA Poisson Gaussian", "Poisson Gaussian with prior mean"), col = c(  "black" , "blue"), pch = 19)
 
 mean((3*sin( 1:n/20)- tt$posterior$mean_log)^2)
 mean((3*sin( 1:n/20)- tt2$posterior$posteriorMean_latent)^2)

## -----------------------------------------------------------------------------
 
 plot(tt2$posterior$posteriorMean_latent, pch=19, col="blue")
points(tt$posterior$mean_log, pch=19)
points(log(x+1),col="red3" ,pch=19)

 lines(log(intensities),lwd=3, col="lightgreen")
 
legend("bottomright", legend = c("log1p", "Split VA Poisson Gaussian", "Poisson Gaussian with prior mean"), col = c("red3", "black" , "blue"), pch = 19)
 mean((3*sin( 1:n/20)- tt$posterior$mean_log)^2)
 mean((3*sin( 1:n/20)- tt2$posterior$posteriorMean_latent)^2)

## -----------------------------------------------------------------------------
s=rep(1, length(x))
init= ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var=NULL),fix_g = c(T,F))
m =init$posterior$mean_log
sigma2 = var(  m - smashr::smash.gaus( m))   

## -----------------------------------------------------------------------------
for ( k in 1:20){
    qb=smashr::smash.gaus(m,sqrt(sigma2), post.var = TRUE)
    Eb = qb$mu.est
    Eb2 = Eb^2+qb$mu.est.var
    
    opt = vga_pois_solver(init_val=m,
                          x=x,s=s,
                          beta= Eb ,
                          sigma2=sigma2,
                          tol=1e-05,
                          maxiter = 100)
    
    
    
    m = opt$m
    v = opt$v
    
    sigma2_new = mean(m^2+v+ Eb2+ 2* Eb -2*m*  Eb  ) 
    
    
}

## -----------------------------------------------------------------------------
plot(log(intensities),lwd=3, col="lightgreen" , type="l")
points(m, pch=19, col="blue") 
points(log(x+1), pch=19 ) 
legend("bottomright", legend = c("log1p of the count", "smooth Split VA Poisson Gaussian" ), col = c( "black" , "blue"), pch = 19)

