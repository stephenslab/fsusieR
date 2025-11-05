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
library(fsusieR)
library(susieR)
library(wavethresh)
set.seed(1)
data(N3finemapping)
attach(N3finemapping)

rsnr <- 0.5 #expected root signal noise ratio

pos1 <- 25    #Position of the causal covariate for effect 1
pos2 <- 25    #Position of the causal covariate for effect 2
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


## -----------------------------------------------------------------------------

noisy.data  <- list()
noisy.data2 <- list()
X <- N3finemapping$X[,1:100]
for ( i in 1:nrow(X))
{
  f1_obs <- f1
  f2_obs <- f2
  noise1 <- rnorm(length(f1), sd=  (1/  rsnr ) * var(f1))
  noise2 <- rnorm(length(f1), sd=  (1/  rsnr ) * var(f2))
  noisy.data [[i]] <-  X[i,pos1]*f1_obs + noise1
  noisy.data2 [[i]] <-   X[i,pos2]*f2_obs + noise2
}
noisy.data <- do.call(rbind, noisy.data)
noisy.data2 <- do.call(rbind, noisy.data2)

Y1 <- noisy.data
Y2 <- noisy.data



## -----------------------------------------------------------------------------
 out1 <- susiF(Y1,X,L=10 )
 out2 <- susiF(Y2,X,L=10 )

## -----------------------------------------------------------------------------
library(coloc)



bf1 <-  exp(out1$lBF[[1]])
names(bf1) <- paste("SNP", rep(1:100))
bf2 <-  exp(out2$lBF[[1]])
names(bf2) <- paste("SNP", rep(1:100))
 
coloc.bf_bf(bf1=bf1, bf2=bf2)

