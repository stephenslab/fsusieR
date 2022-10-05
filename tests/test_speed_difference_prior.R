library(susiF.alpha)
library(ashr)
library(microbenchmark)
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
sim  <- simu_test_function(rsnr=1,pos2= 2 ,is.plot = FALSE)
Y <- sim$noisy.data
X <- sim$G

microbenchmark(
  v1= susiF(Y,X,L=10 , prior = 'mixture_normal_per_scale'),
  v2= susiF(Y,X,L=10 , prior = 'mixture_normal')
)



