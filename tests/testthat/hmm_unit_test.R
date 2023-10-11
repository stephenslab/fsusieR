##unit test HMM

low variance


sharp transition


rm(list=ls())
devtools::load_all(".")
source("D:/Document/Serieux/Travail/Package/susiF.alpha/inst/code/fit_hmm.R", echo=TRUE)

library(susiF.alpha)
library(ashr)
library(wavethresh)
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
sd_noise <- 0.1 #expected root signal noise ratio
N <- 100    #Number of individuals
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
f1[ 70] <- -1
#first effect)
f2 <-  0.1*DJ.EX(128)$blocks

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
Y[,10:20] <- 0*Y[,10:20]
out <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale', filter.number =8  )


X <- colScale(X)
# centering input
Y <- colScale(Y, scale=FALSE)
susiF.obj <- out
L_points=20

idx <- do.call( c, lapply( 1:length(susiF.obj$cs),
                           function(l){
                             tp_id <-  which.max( susiF.obj$pip[susiF.obj$cs[[l]]])
                             susiF.obj$cs[[l]][tp_id]
                           }
)
)

temp_Y <- Y
res <- cal_Bhat_Shat(temp_Y,X )



s =fit_hmm(x=res$Bhat[idx[1],],sd=res$Shat[idx[1],],halfK=50 )
