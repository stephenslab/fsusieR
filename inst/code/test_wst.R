library (susiF.alpha)
library(ashr)
library(wavethresh)
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
rsnr <- 2 #expected root signal noise ratio
N <- 100    #Number of individuals
P <- 10     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 5   #Position of the causal covariate for effect 2
lev_res <- 10#length of the molecular phenotype (2^lev_res)
f1 <-  0*simu_IBSS_per_level(lev_res )$sim_func#first effect
f2 <- 0*DJ.EX(2^7)$bumps #second effect

f1 [50 ] <-1
f2 [100 ] <-1

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
  noise <- rnorm(length(f1), sd=  1)
  noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs + beta2*G[i,pos2]*f2_obs + noise

}
noisy.data <- do.call(rbind, noisy.data)




plot( noisy.data[1,], type = "l", col=(G[1, pos1]*3+1),
      main="Observed curves \n colored by the causal effect", ylim= c(-40,40), xlab="")
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

out <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale',TI=FALSE)
out2 <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale'  )



plot_susiF(out2, cred.band = TRUE)



#### Findout which regions are affected by the different CS

#This corresponds to the regions where the credible bands are "crossing zero"/i.e. the effects are likely not to be 0 in this region.

#You can also access the information directly in the output of susiF  as follow
par(mfrow=c(1,1))


plot( f1, type="l", main="Estimated effect 1", xlab="", ylim=c(-0.3,1.3), xlim=c(40,200))
lines(unlist(out$fitted_func[[2]]),col='blue' )
lines(unlist(out2$fitted_func[[2]]),col='green' )
legend(x= 35,
       y=3,
       lty= rep(1,2),
       legend = c("effect 1"," fSuSiE est "),
       col=c("black","blue" )
)
plot( f1, type="l", main="Estimated effect 2", xlab="" , ylim=c(-0.3,1.3), xlim=c(40,200), lwd=2)
lines(unlist(out$fitted_func[[2]]),col='blue',lwd=1.3 )
lines(unlist(out2$fitted_func[[2]]),col='darkgreen',lwd=1.3  )
abline(a=0,b=0)
legend(x= 20,
       y=-1.5,
       lty= rep(1,2),
       legend = c("effect 2"," fSuSiE est "),
       col=c("black","green" )
)


par(mfrow=c(1,1))
