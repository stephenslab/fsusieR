rm(list=ls())
library(susiF.alpha)
lst <- list()
est_lst<- list()
for ( o in 1:100){

  #Example using curves simulated under the Mixture normal per scale prior
  rsnr <- runif(n=1, min=0.1,max=2) #expected root signal noise ratio
  N <- 100    #Number of individuals
  P <- 3     #Number of covariates/SNP
  pos1 <- 1   #Position of the causal covariate for effect 1
  pos2 <- 2   #Position of the causal covariate for effect 2
  lev_res <- 7#length of the molecular phenotype (2^lev_res)
  f1 <-  simu_IBSS_per_level(lev_res )$sim_func#first effect
  f2 <- simu_IBSS_per_level(lev_res )$sim_func #second effect


  G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
  beta0       <- 0
  beta1       <- 1
  beta2       <- 1
  noisy.data  <- list()

  for ( i in 1:N)
  {
    f1_obs <- f1
    f2_obs <- f2
    noise <- rnorm(length(f1), sd=  (1/  rsnr ) * var(f1))
    noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs + beta2*G[i,pos2]*f2_obs + noise

  }
  noisy.data <- do.call(rbind, noisy.data)



  Y <- noisy.data
  X <- G


  res <- EBmvFR (Y=noisy.data,X=X)
  res2 <- EBmvFR2 (Y=noisy.data,X=X)

  est_lst [[o]] <- sqrt(res$sigma2)
  lst [[o]]     <- (1/  rsnr ) * var(f1)

}


plot( do.call(c,est_lst), do.call(c, lst) )

plot( do.call(c,est_lst), do.call(c, lst), xlim=c(0,8),ylim=c(0,8))
abline(a=0,b=1)
