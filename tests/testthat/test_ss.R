library(testthat)
library(ashr)
library(wavethresh)
library(mixsqp)
N =100
P=10
set.seed(23)
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1 <- 5
noisy.data  <- list()
set.seed(2)
f1 <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay =1.5)

plot(f1$sim_func, type="l", ylab="y")
N=500
P=10
set.seed(23)
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1 <- 5
noisy.data  <- list()
rsnr=10
for ( i in 1:N)
{
  f1_obs <- f1$sim_func
  noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs +  rnorm(length(f1$sim_func), sd=  (1/  rsnr ) *sd(f1$sim_func))

}
noisy.data <- do.call(rbind, noisy.data)


Y <- noisy.data
X <- G
W <- DWT2(Y)
update_D <- W
Y_f <- cbind( W$D,W$C) #Using a column like phenotype
update_Y <-Y_f
v1 <- rep(1, dim(X)[2])
tt <-  cal_Bhat_Shat(Y_f,X,v1)
indx_lst <- gen_wavelet_indx(9)
Bhat <- tt$Bhat
Shat <- tt$Shat

XtX <- t(X)%*%X
Xty <- t(X)%*%Y_f
yty <- t(Y_f)%*%Y_f
var_y <- apply(Y_f, 2,var)
R <- cor(X)
data_suff <- make_data_suff_stat (Bhat , Shat , R , N  , var_y , XtX , Xty , yty )
class(data_suff)

test_that("Class of data_suff is", {

  data_suff <- make_data_suff_stat (Bhat , Shat , R , N  , var_y , XtX , Xty , yty )

  expect_equal(  class(data_suff)[1],
  "suff_stat"
  )
})
indx_lst <- gen_wavelet_indx(9)

G_prior <- init_prior(data_suff, prior="mixture_normal_per_scale",  indx_lst )
test_that("Class of the prior is", {
  G_prior <- init_prior(data_suff, prior="mixture_normal_per_scale",  indx_lst )
  expect_equal(class( G_prior)[1],
  "mixture_normal_per_scale"
  )
})



susiF_ss.obj   <-  init_susiF_ss_obj(L=2, G_prior, data_suff )
class(susiF_ss.obj)

test_that("Class of the prior is", {

  susiF_ss.obj   <-  init_susiF_ss_obj(L=2, G_prior, data_suff )

  expect_equal(class(susiF_ss.obj)[1],
               "susiF_ss"
  )
})
test_that("Class of the pis should be", {


  expect_equal(class(get_pi(susiF_ss.obj,l=1))[1],
               "pi_mixture_normal_per_scale"
  )
})


get_pi(susiF_ss.obj,l=1)

susiF_ss  (Bhat, Shat, R, N , var_y, XtX, Xty, yty, L = 2,
                     pos = NULL,
                     prior = "mixture_normal_per_scale",
                     verbose = TRUE,
                     plot_out = TRUE,
                     maxit = 100,
                     tol = 1e-3,
                     cov_lev = 0.95

)
