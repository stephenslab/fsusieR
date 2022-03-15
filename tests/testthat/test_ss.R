library(testthat)
library(susiF.alpha)
library(ashr)
library(wavethresh)
library(mixsqp)
N =1000
P=10
set.seed(23)
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1 <- 5
noisy.data  <- list()
set.seed(2)
f1 <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay = .5)

plot(f1$sim_func, type="l", ylab="y")
N=200
P=10
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1        <- 5
noisy.data  <- list()
rsnr        =  0.02
for ( i in 1:N)
{
  f1_obs           <- f1$sim_func
  noisy.data [[i]] <- beta1*G[i,pos1]*f1_obs +  rnorm(length(f1$sim_func), sd= 1)
}
noisy.data <- do.call(rbind, noisy.data)


Y <- noisy.data
X <- G

Y <- apply(Y, 2,scale)

v1 <- rep(1,nrow(Y) )


X <- apply(X, 2,scale)

W <- DWT2(Y)
update_D <- W
Y_f <- cbind( W$D,W$C) #Using a column like phenotype

tt <-  cal_Bhat_Shat(Y ,X,v1)
Bhatb <- tt$Bhat
Shatb <- tt$Shat


update_Y <- Y_f
v1       <- rep(1, dim(X)[2])
tt       <- cal_Bhat_Shat(Y_f,X,v1)
indx_lst <- gen_wavelet_indx(9)
Bhat     <- tt$Bhat
Shat     <- tt$Shat

XtX   <- t(X)%*%X
Xty   <- t(X)%*%Y_f
yty   <- t(Y_f)%*%Y_f
var_y <- apply(Y_f, 2,var)
R     <- cor(X)
data_suff  <- make_data_suff_stat (Bhat , Shat , R , N  , var_y , XtX , Xty , yty,wav_trans = TRUE )
data_suffb <- make_data_suff_stat (Bhatb , Shatb , R , N  , var_y , XtX , Xty , yty)

class(data_suff)


plot( data_suff$Bhat, data_suffb$Bhat)
abline(a=0,b=1)
plot( data_suff$Shat, data_suffb$Shat)
abline(a=0,b=1)

plot( data_suff$Bhat/data_suff$Shat, data_suffb$Bhat/data_suffb$Shat)
abline(a=0,b=1)


plot( data_suff$Bhat/data_suff$Shat, data_suffb$Bhat/data_suffb$Shat, xlim=c(-10,10), ylim=c(-10,10))
abline(a=0,b=1)

test_that("Class of data_suff is", {

  data_suff <- make_data_suff_stat (Bhat , Shat , R , N  , var_y , XtX , Xty , yty )

  expect_equal(  class(data_suff)[1],
  "suff_stat"
  )
})

test_that("The transform coefficient should be equal to", {

  data_suff  <- make_data_suff_stat (Bhat , Shat , R , N  , var_y , XtX , Xty , yty,wav_trans = TRUE )
  data_suffb <- make_data_suff_stat (Bhatb , Shatb , R , N  , var_y , XtX , Xty , yty)
  expect_equal( data_suff$Bhat, data_suffb$Bhat)
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


test_that("Correct updating of the susiF_ss object",{

  EM_pi  <- EM_pi(G_prior, data_suff$Bhat, data_suff$Shat, indx_lst = indx_lst)
  susiF_ss.obj <-   update_susiF_obj(susiF_ss.obj,l=1, EM_pi, data_suff$Bhat, data_suff$Shat, indx_lst = indx_lst)
  expect_equal(class(susiF_ss.obj), "susiF_ss")

})



test_that("Correct updating of the susiF_ss object",{
  update_data <- get_partial_residual(susiF_ss.obj,update_data,l=1)
  tt <- cal_Bhat_Shat(susiF_ss.obj,update_data,partial=TRUE)

  expect_equal(Bhat, tt$Bhat)

})

get_ER2(susiF_ss.obj,data_suff )


estimate_residual_variance(susiF_ss.obj,data = data_suff)


loglik_SFR_post(susiF_ss.obj,data = data_suff)
Eloglik(susiF_ss.obj,data = data_suff)


cal_KL_l(susiF_ss.obj, 1, data)

lol  <- update_KL.susiF_ss  (susiF_ss.obj, data  )
lol$KL

get_pi(susiF_ss.obj,l=1)
data <-data_suff

EM_pi  <- EM_pi(G_prior, data_suff$Bhat, data_suff$Shat, indx_lst = indx_lst)
l=1
susiF_ss.obj         <-   update_pi(susiF_ss.obj   = susiF_ss.obj   ,
                                    l = l ,
                                    tpi =  EM_pi$tpi_k)
susiF_ss.obj$G_prior <-   update_prior(get_G_prior(susiF_ss.obj  ) , EM_pi$tpi_k  )

susiF_ss.obj$fitted_wc[[l]]   <- post_mat_mean(get_G_prior(susiF_ss.obj) , Bhat=data$Bhat, Shat=data$Shat,indx_lst= indx_lst )
susiF_ss.obj$fitted_wc2[[l]]  <- post_mat_sd  (get_G_prior(susiF_ss.obj) , Bhat=data$Bhat, Shat=data$Shat, indx_lst= indx_lst)^2


new_alpha <- cal_zeta(  EM_pi$lBF)
susiF_ss.obj <- update_alpha(susiF_ss.obj, l, new_alpha)
susiF_ss.obj <- update_lBF(susiF_ss.obj, l, EM_pi$lBF)




get_post_F(susiF_ss.obj,1)

get_post_F2(susiF_ss.obj,1)


cal_expected_residual (susiF_ss.obj , data)


update_data <- data_suff



update_data <- cal_expected_residual(susiF_ss.obj,data)

update_data <- get_partial_residual(susiF_ss.obj,update_data,l=1)
tt <- cal_Bhat_Shat(susiF_ss.obj,update_data,partial=TRUE)
plot( Bhat, tt$Bhat)


plot(c(Shat), c(tt$Shat))
hist( Shat
hist( tt$Shat )





Bhat <- tt$Bhat
Shat <- tt$Shat #UPDATE. could be nicer
tpi <-  get_pi(susiF.obj,l)
G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

EM_out  <- EM_pi(G_prior  = G_prior,
                 Bhat     =  Bhat,
                 Shat     =  Shat,
                 indx_lst =  indx_lst
)





susiF_ss.obj <-  update_susiF_obj(susiF_ss.obj = susiF_ss.obj ,
                                  l         = l,
                                  EM_pi     = EM_out,
                                  Bhat      = Bhat,
                                  Shat      = Shat,
                                  indx_lst  = indx_lst
)
