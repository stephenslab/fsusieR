library(ashr)
library(wavethresh)
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
rsnr <- 0.2 #expected root signal noise ratio
N <- 100    #Number of individuals
P <- 100     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 5   #Position of the causal covariate for effect 2
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
       col=c("black","blue","yellow"))
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
  noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs + beta2*G[i,pos2]*f2_obs+ beta1*G[i,15]*f1_obs + noise
  
}
noisy.data <- do.call(rbind, noisy.data)




plot( noisy.data[1,], type = "l", col=(G[1, pos1]*3+1),
      main ="Observed curves \n   colored by the causal effect" , ylim= c(-40,40), xlab="")
for ( i in 2:N)
{
  lines( noisy.data[i,], type = "l", col=(G[i, pos1]*3+1))
  
}
legend(x=0.3,
       y=-10,
       lty = rep(1,3),
       legend= c("0", "1","2"),
       col=c("black","blue","yellow"))



Y <- noisy.data
X <- G
#Running fSuSiE

out <- susiF(Y,X,L=3 , prior = 'mixture_normal_per_scale')
out2 <- susiF(Y,X,L=3 , prior = 'mixture_normal_per_scale',
             post_processing = "HMM")

test_that("some_function does not crash", {
  expect_error(plot_susiF_effect(out,effect = c(1,3),cred_band = FALSE), NA)
  expect_error(plot_susiF_effect(out,effect = c(1,3)), NA)
  expect_error(plot_susiF_effect(out), NA)
  expect_error(plot_susiF_pip(out), NA)
  expect_error(plot_susiF(out), NA)
  expect_error(plot(out), NA)
  expect_error(plot_susiF_effect(out2,effect = c(1,3),lfsr_curve  = FALSE), NA)
  expect_error(plot_susiF_effect(out2,effect = c(1,3)), NA)
  expect_error(plot_susiF_effect(out2), NA)
  expect_error(plot_susiF_pip(out2), NA)
  expect_error(plot_susiF(out2), NA)
  expect_error(plot(out2), NA)
  
  
}) 
