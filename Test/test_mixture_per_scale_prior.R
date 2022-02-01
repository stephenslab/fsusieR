library(testthat)
library(ashr)
library(wavethresh)
library(Rfast)
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
susiF_obj <- init_susiF_obj(L=1, G_prior,Y,X)


test_that("Susif object pi are expected to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L=1, G_prior,Y,X)

            expect_equal(get_pi(susiF_obj,1), get_pi_G_prior(G_prior)
            )
          }
)


test_that("Susif internal prior to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L=1, G_prior,Y,X)

            expect_equal(get_G_prior (susiF_obj ),  G_prior)

          }
)

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
plot( Bhat,  post_mat_mean(G,Bhat,Shat, indx_lst) )
plot( Shat,  (post_mat_sd(G,Bhat,Shat, indx_lst) ))
test_that("Class of the proportions  is", {

  expect_equal(class(get_pi_G_prior(G))
  ,
  "pi_mixture_normal_per_scale"
  )
})
test_that("Class of the standard deviations  is", {

  expect_equal(class(get_sd_G_prior(G))
               ,
               "sd_mixture_normal_per_scale"
  )
})
get_pi_G_prior(G)
get_sd_G_prior(G)


L <- L_mixsq(G, Bhat, Shat, indx_lst)


test_that("The likelihood computed by L_mixsqp should be of class", {
  L <- L_mixsq(G, Bhat, Shat, indx_lst)

  expect_equal(class(L), "lik_mixture_normal_per_scale"
  )
})




zeta <- cal_zeta(lBF)
test_that("The highest assignation should be equal to", {
  zeta <- cal_zeta(lBF)
  expect_equal(which.max(zeta), pos1
  )
})

tpi <- m_step(L, zeta , indx_lst)
class(tpi)

test_that("The output of the m_step function should of the class", {
  tpi <- m_step(L, zeta , indx_lst)

  expect_equal( class(tpi),"pi_mixture_normal_per_scale"
  )
})

test_that("The output of the m_step for the pi_0 should equal", {
  tpi <- m_step(L, zeta , indx_lst)
  expect_equal( get_pi0(tpi = tpi), c(0,0.5, rep(1, 8)),
                tolerance = 0.01) #allow 1% error in the proportion estimation
})


G_prior <- G
G_update <- update_prior (G_prior, tpi)
test_that("Updated mixture proportion should be equal to provided input",
          {
            expect_equal(identical(get_pi_G_prior(G_update) ,tpi),
                         TRUE
            )
          }
)
outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
test_that("The outputs of the EM_pi function should be  ",
          {
           outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            expect_equal(class(outEM$tpi_k) ,"pi_mixture_normal_per_scale")
            expect_equal(class(outEM$lBF) ,"numeric")
            expect_equal(length(outEM$lBF) ,dim(Bhat)[1])
            expect_equal( get_pi0(tpi = outEM$tpi_k), c(0,0.5, rep(1, 8)),
                          tolerance = 0.01)
          }
)

test_that("The mixture proportion of the updated susiF object should be equal to   ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            susiF_obj <- update_pi_susiF( susiF_obj, 1,  outEM$tpi_k)

            expect_equal( get_pi (susiF_obj , 1),outEM$tpi_k )
          }
)


test_that("The alpha value of  the update susiF object should be equal to   ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            new_alpha <- cal_zeta(outEM$lBF)
            susiF_obj <- update_alpha.susiF(susiF_obj, l, new_alpha)
            expect_equal( get_alpha (susiF_obj , 1), new_alpha )
          }
)

test_that("The update susiF object should have its argument equal to    ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )

            expect_equal( susiF_obj$fitted_wc[[1]],post_mat_mean( G_prior , Bhat, Shat, indx_lst ))
            expect_equal( susiF_obj$fitted_wc2[[1]],post_mat_sd  ( G_prior , Bhat, Shat , indx_lst))
            expect_equal( get_alpha (susiF_obj , 1), cal_zeta(outEM$lBF))
            expect_equal( get_G_prior(susiF_obj) ,G_prior)

          }
)
