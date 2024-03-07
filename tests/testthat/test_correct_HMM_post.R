rm(list=ls())
devtools::load_all(".")

library(fsusieR)
library(ashr)
library(wavethresh)

#Example using curves simulated under the Mixture normal per scale prior
sd_noise <- 3 #expected root signal noise ratio
N <- 100     #Number of individuals
P <- 10     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 5   #Position of the causal covariate for effect 2
lev_res <-7#length of the molecular phenotype (2^lev_res)
f1 <-  rep(0, 2^lev_res)
f1[ 20:25] <-2
f1[ 50:55] <-1
f1[ 20:25] <-2

f1 <-  0*simu_IBSS_per_level(lev_res )$sim_func
f1[60:length(f1)] <-0
f1[ 70] <- -2
#first effect)
f2 <-  0.1*DJ.EX(128)$bumps
plot( f2, type ="l", ylab="effect", col="blue")
lines(f1, col="red")
beta0       <- 0
beta1       <- 1
beta2       <- 1
noisy.data  <- list()
#f2 <-f1
G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
X <-G
for ( i in 1:N)
{
  noise <- rnorm(length(f1), sd=  sd_noise)
  noisy.data [[i]] <-   beta1*G[i,pos1]*f1 + beta1*G[i,pos2]*f2  + noise
  
}
noisy.data <- do.call(rbind, noisy.data)
Y <- noisy.data

out <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale',
               
             post_processing = "HMM" )
 
out1 <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale'   )
lines(out$fitted_func[[1]])


susiF.obj <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale',
             
             post_processing = "HMM" )
 

devtools::load_all(".")

Xb    = colScale(X)
Yb      <- colScale(Y, scale=FALSE)
out2  <- HMM_regression( obj=out1,
                        Y=Yb,
                        X=X 
)

#seems related to the variance scale of X


plot (out  $fitted_func[[1]])
lines(susiF.obj$fitted_func[[1]])

lines(out2$fitted_func[[1]])
 

plot(out  $fitted_func[[1]],out2$fitted_func[[1]] )
abline(a=0,b=1)

plot (out  $lfsr_func[[1]])
lines(susiF.obj$lfsr_func[[1]])

lines(out2$lfsr_func[[1]])

plot(out  $lfsr_func[[1]],out2$lfsr_func[[1]] )
abline(a=0,b=1)
