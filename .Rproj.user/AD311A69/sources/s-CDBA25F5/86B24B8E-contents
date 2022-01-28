library(testthat)
library(ashr)
library(wavethresh)
library(Rfast)
set.seed(2)
f1 <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay =1.5)

plot(f1$sim_func, type="l", ylab="y")
N=10
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
indx_lst <- gen_wavelet_indx(9)
lol <-  cal_Bhat_Shat(Y_f,X,v1)
test_that("Class of the prior is", {
  expect_equal(class(
                      init_prior(Y=Y_f,
                                 X=X,
                                 prior="mixture_normal",
                                 v1=v1,
                                 indx_lst = indx_lst)
                      ),
               "mixture_normal"
               )
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
tt <- cal_Bhat_Shat(Y_f,X,v1)


### Test validity normal mixture  -----
Bhat <- tt$Bhat
Shat <- tt$Shat
G <- init_prior(Y=Y_f,
                X=X,
                prior="mixture_normal",
                v1=v1,
                indx_lst = indx_lst)
G_prior <- G
lBF <- log_BF (G_prior, tt$Bhat, tt$Shat )
lBF
test_that("Max lBF should be in postion",
          {
  expect_equal(which.max(lBF),
              pos1
              )
          }
         )
plot( Bhat,post_mat_mean(G,Bhat,Shat))



plot( Shat,  (post_mat_sd(G,Bhat,Shat, indx_lst) ))
get_pi_G_prior(G)
get_sd_G_prior(G)
L <- L_mixsq(G, Bhat, Shat)
L
zeta <- cal_zeta(lBF)
tpi <- m_step(L, zeta , indx_lst)
tpi
G_update <- update_prior (G_prior, tpi)
test_that("Updated mixture proportion should be equal to provided input",
          {
            expect_equal(identical(get_pi_G_prior(G_update) ,tpi),
                         TRUE
            )
          }
)
### Test validity normal mixture per scale -----
G <- init_prior(Y=Y_f,
                X=X,
                prior="mixture_normal_per_scale",
                v1=v1,
                indx_lst = indx_lst)
lBF <- log_BF (G, tt$Bhat, tt$Shat , indx_lst)
lBF
test_that("Max lBF should be in postion",
          {
            expect_equal(which.max(lBF),
                         pos1
            )
          }
)
plot( Bhat,  post_mat_mean(G,Bhat,Shat, indx_lst) )
plot( Shat,  (post_mat_sd(G,Bhat,Shat, indx_lst) ))
get_pi_G_prior(G)
get_sd_G_prior(G)
L <- L_mixsq(G, Bhat, Shat, indx_lst)
lBF
zeta <- cal_zeta(lBF)
tpi <- m_step(L, zeta , indx_lst)
class(tpi)
G_prior <- G
G_update <- update_prior (G_prior, tpi)
test_that("Updated mixture proportion should be equal to provided input",
          {
            expect_equal(identical(get_pi_G_prior(G_update) ,tpi),
                         TRUE
            )
          }
)

