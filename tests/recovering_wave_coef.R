rm(list=ls())
library(testthat)
library(ashr)
library(wavethresh)
library(Rfast)
library(mixsqp)
library(susiF.alpha)
set.seed(2)
f1 <- simu_IBSS_per_level(lev_res=4, alpha=1, prop_decay =1.5)

plot(f1$sim_func, type="l", ylab="y")
N=100
P=10
set.seed(23)
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1 <- 5
noisy.data  <- list()
rsnr=1
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

tt2 <-  cal_Bhat_Shat(Y,X,v1)
Bhat2 <- tt2$Bhat
Shat2 <- tt2$Shat


W1 <-t(GenW(n=  ncol(Shat2)  , filter.number = 10, family = "DaubLeAsymm"))
# Shat stack C coefficient last
W1[, c(1,ncol(W1 ))]<- W1[, c(ncol(W1),1)]

W <- DWT2(Bhat2 )
WBhat2 <-   cbind( W$D,W$C)
image(WBhat2)
image(Bhat)
plot( Bhat,WBhat2 )

hist(WBhat2-Bhat)


W1 <- (GenW(n=  ncol(Shat2)  , filter.number = 10, family = "DaubLeAsymm"))
W1[, c(1,ncol(W1 ))]<- W1[, c(ncol(W1),1)]

WShat2 <- 0*Shat2
rWShat2 <- 0*Shat2
for ( i in 1: nrow(Shat2)){
  WShat2[i, ]<- diag(t(W1)%*% diag( Shat2[i, ]  )%*% (W1))
  rWShat2[i, ]<- diag(  t(W1)%*% diag(    (Shat[i, ]  ))%*%  (W1))
}


image(WShat2)
image(rWShat2)
image(Shat )
image(rWShat2-Shat)

image(Shat)
image(Shat2)
plot( Shat,WShat2)
abline(a=0,b=1)




G_prior<- init_prior(Y=Y_f,
                     X=X,
                     prior="mixture_normal",
                     v1=v1,
                     indx_lst = indx_lst)

lBF <- log_BF (G_prior, Bhat, Shat , indx_lst)
lBF
tt <-  cal_Bhat_Shat(Y_f,X,v1)
Bhat <- tt$Bhat
Shat <- tt$Shat
test_that("Class of the prior is", {

  expect_equal(class(
    init_prior(Y=Y_f,
               X=X,
               prior="mixture_normal_per_scale",
               v1=v1,
               indx_lst = indx_lst)
  ),
  "mixture_normal_per_scale"
  )
})

### Test validity normal mixture  -----
Bhat <- tt$Bhat
Shat <- tt$Shat
G_prior<- init_prior(Y=Y_f,
                     X=X,
                     prior="mixture_normal",
                     v1=v1,
                     indx_lst = indx_lst)


plot( Bhat,post_mat_mean(G_prior,Bhat,Shat))



plot( Shat,  (post_mat_sd(G_prior,Bhat,Shat, indx_lst) ))
get_pi_G_prior(G_prior)
get_sd_G_prior(G_prior)
L <- L_mixsq(G_prior, Bhat, Shat)
L
zeta <- cal_zeta(lBF)
tpi <- m_step(L, zeta , indx_lst)
tpi
