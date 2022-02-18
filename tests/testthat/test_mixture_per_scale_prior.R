library(testthat)
library(ashr)
library(wavethresh)
library(mixsqp)
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
G_prior <- init_prior(Y=Y_f,
                      X=X,
                      prior="mixture_normal_per_scale",
                      v1=v1,
                      indx_lst = indx_lst)

lBF <- log_BF (G_prior, tt$Bhat, tt$Shat , indx_lst)
lBF
test_that("Max lBF should be in postion",
          {
            expect_equal(which.max(lBF),
                         pos1
            )
          }
)
susiF_obj <- init_susiF_obj(L=2, G_prior,Y,X)
test_that("susiF object pi are expected to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L=2, G_prior,Y,X)

            expect_equal(get_pi(susiF_obj,1), get_pi_G_prior(G_prior)
            )
            expect_equal(get_pi(susiF_obj,2), get_pi_G_prior(G_prior)
            )
          }
)

susiF_obj <- init_susiF_obj(L=1, G_prior,Y,X)


test_that("susiF object pi are expected to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L=1, G_prior,Y,X)

            expect_equal(get_pi(susiF_obj,1), get_pi_G_prior(G_prior)
            )
          }
)


test_that("susiF internal prior to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L=2, G_prior,Y,X)

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
  )[1],
  "mixture_normal_per_scale"
  )
})
plot( Bhat,  post_mat_mean(G_prior,Bhat,Shat, indx_lst) )
plot( Shat,  (post_mat_sd(G_prior,Bhat,Shat, indx_lst) ))
test_that("Class of the proportions  is", {

  expect_equal(class(get_pi_G_prior(G_prior))[1]
               ,
               "pi_mixture_normal_per_scale"
  )
})
test_that("Class of the standard deviations  is", {

  expect_equal(class(get_sd_G_prior(G_prior))[1]
               ,
               "sd_mixture_normal_per_scale"
  )
})



L <- L_mixsq(G_prior, Bhat, Shat, indx_lst)


test_that("The likelihood computed by L_mixsq  should be of class", {
  L <- L_mixsq(G_prior, Bhat, Shat, indx_lst)

  expect_equal(class(L)[1], "lik_mixture_normal_per_scale"
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

  expect_equal( class(tpi)[1],"pi_mixture_normal_per_scale"
  )
})

test_that("The output of the m_step for the pi_0 should equal", {
  tpi <- m_step(L, zeta , indx_lst)
  expect_equal( get_pi0(tpi = tpi), c(0,0.5, rep(1, 8)),
                tolerance = 0.01) #allow 1% error in the proportion estimation
})


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
            expect_equal(class(outEM$tpi_k)[1] ,"pi_mixture_normal_per_scale")
            expect_equal(class(outEM$lBF)[1] ,"numeric")
            expect_equal(length(outEM$lBF)[1] ,dim(Bhat)[1])
            expect_equal( get_pi0(tpi = outEM$tpi_k), c(0,0.5, rep(1, 8)),
                          tolerance = 0.01)
          }
)

test_that("The mixture proportion of the updated susiF object should be equal to   ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            susiF_obj <- update_pi( susiF_obj, 1,  outEM$tpi_k)

            expect_equal( get_pi (susiF_obj , 1),outEM$tpi_k )
          }
)


test_that("The alpha value of  the update susiF object should be equal to   ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            new_alpha <- cal_zeta(outEM$lBF)
            susiF_obj <- update_alpha(susiF_obj, 1, new_alpha)
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
            expect_equal( susiF_obj$fitted_wc2[[1]],post_mat_sd  ( G_prior , Bhat, Shat , indx_lst)^2)
            expect_equal( get_alpha (susiF_obj , 1), cal_zeta(outEM$lBF))
            expect_equal( get_G_prior(susiF_obj) ,G_prior)
            expect_equal(   susiF_obj$lBF[[1]]  , outEM$lBF )
          }
)

test_that("The partial residual should be    ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )

            update_T <- cal_partial_resid(
              susiF.obj = susiF_obj,
              l         = 1,
              X         = X,
              D         = W$D,
              C         = W$C,
              indx_lst  = indx_lst
            )

            L=2
            l=1
            id_L <- (1:L)[ - ( (l%%L)+1) ]
            update_D  <-  W$D - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF_obj$alpha[[l]], rep.int(N,P))) %*% (susiF_obj$fitted_wc[[l]][,-indx_lst[[length(indx_lst)]]])   ) )
            update_C  <-  W$C - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF_obj$alpha[[l]], rep.int(N,P))) %*% susiF_obj$fitted_wc[[l]][,indx_lst[[length(indx_lst)]]] ) )
            manual_update <- cbind(  update_D, update_C)

            expect_equal(  update_T ,manual_update)

          }
)


test_that("The output update should be equal to    ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )
            tcs <- list()
            tpip <- list()
            for ( l in 1:susiF_obj$L)
            {
              temp        <- susiF_obj$alpha[[l]]
              temp_cumsum <- cumsum( temp[order(temp, decreasing =TRUE)])
              max_indx_cs <- min(which( temp_cumsum >0.95))
              max_indx_cs <- min(which( temp_cumsum >0.95))
              tcs[[l]]  <- order(temp, decreasing = TRUE)[1:max_indx_cs ]
              tpip[[l]] <- rep(1, lengths(susiF_obj$alpha)[[l]])-susiF_obj$alpha[[l]]
            }
            pip <- 1-  apply( do.call(rbind,tpip),2, prod)


            fitted_func <- list ()
            temp <- wd(rep(0, dim(Y_f)[2]))
            for ( l in 1:susiF_obj$L)
            {
              temp$D <-    (susiF_obj$alpha[[l]])%*%susiF_obj$fitted_wc[[l]][,-indx_lst[[length(indx_lst)]]]
              temp$C[length(temp$C)] <- (susiF_obj$alpha[[l]])%*%susiF_obj$fitted_wc[[l]][,indx_lst[[length(indx_lst)]]]
              fitted_func[[l]] <- wr(temp)
            }


            ind_fitted_func  <- matrix(0, ncol=dim(Y)[2], nrow=dim(Y)[1])
            for ( i in 1:dim(Y)[1])
            {
              ind_fitted_func[i,]  <- rep(0,dim(Y)[2])#fitted_baseline
              for ( l in 1:susiF_obj$L)
              {
                #add wavelet coefficient
                temp$D                         <-    (susiF_obj$alpha[[l]] *X[i,])%*%susiF_obj$fitted_wc[[l]][,-indx_lst[[length(indx_lst)]]]
                temp$C[length(temp$C)]         <-    (susiF_obj$alpha[[l]] *X[i,])%*%susiF_obj$fitted_wc[[l]][,indx_lst[[length(indx_lst)]]]
                #transform back
                ind_fitted_func[i,]  <-  ind_fitted_func[i,]+wr(temp)
              }
            }


            expect_equal(  update_cal_pip(susiF_obj)$pip                              ,pip)
            expect_equal(  update_cal_cs(susiF_obj)$cs                                ,tcs)
            expect_equal(  update_cal_indf(susiF_obj, Y=Y, X, indx_lst)$ind_fitted_func ,ind_fitted_func)
            expect_equal(  update_cal_fit_func(susiF_obj, indx_lst)$fitted_func       ,fitted_func)
            expect_equal(  out_prep(susiF_obj,Y=Y, X=X, indx_lst=indx_lst)$pip             ,pip)
            expect_equal(  out_prep(susiF_obj,Y=Y, X=X, indx_lst=indx_lst)$cs              ,tcs)
            expect_equal(  out_prep(susiF_obj,Y=Y, X=X, indx_lst=indx_lst)$ind_fitted_func ,ind_fitted_func)
            expect_equal(  out_prep(susiF_obj,Y=Y, X=X, indx_lst=indx_lst)$fitted_func     ,fitted_func)

          }
)



test_that("The precision of the fitted curves should be   ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )
            expect_equal(  sum( abs(unlist(update_cal_fit_func(susiF_obj, indx_lst)$fitted_func) -f1$sim_func)), 0, tol=0.004)


          }
)

outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
G_prior <- update_prior(G_prior,
                        tpi= outEM$tpi_k )

susiF_obj <-  update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )
susiF_obj <-  out_prep(susiF_obj,Y, X=X, indx_lst=indx_lst)

plot( unlist(susiF_obj$fitted_func), type="l", col="green")
lines(f1$sim_func, col="red")



test_that("SusiF performance should be",
          {
            set.seed(1)
            sim  <- simu_test_function(rsnr=2,is.plot = FALSE)
            Y <- sim$noisy.data
            X <- sim$G
            out <- susiF(Y,X,L=1, prior="mixture_normal_per_scale")
            expect_equal(  unlist( out$alpha) , c(1, rep(0,9)) , tol=1e-5)
            expect_equal(  sum( abs(unlist(out$fitted_func) -sim$f1)), 0, tol=0.2*length(sim$f1))

          }
)



test_that("The expected sum of square should be below and residual variance should be ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )
            expect_lte(  get_ER2( susiF_obj, Y_f,X ), 0.3)
            sigma2 <- estimate_residual_variance(susiF_obj,Y_f,X)
            susiF_obj <- update_residual_variance(susiF_obj, sigma2 = sigma2 )
            expect_equal(  susiF_obj$sigma2, sigma2)


          }
)


test_that("The KL of effect one ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst)
            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )
            KL_l <- cal_KL_l       (susiF_obj,l=1,Y=Y_f, X, D=W$D, C=W$C , indx_lst)
            susiF_obj <- update_KL ( susiF_obj,   Y=Y_f, X, D=W$D, C=W$C , indx_lst)
            expect_equal( susiF_obj$KL[1], KL_l)
            get_objective(susiF_obj,Y_f,X, D=W$D, C=W$C , indx_lst)

          }
)



test_that("SusiF performance should be",
          {
            set.seed(1)
            sim  <- simu_test_function(rsnr=2,pos2= 2 ,is.plot = FALSE)
            Y <- sim$noisy.data
            X <- sim$G
            out <- susiF(Y,X,L=2, prior="mixture_normal_per_scale")
            expect_equal(  Reduce("+", out$alpha) , c(1, 1,rep(0,8)) , tol=1e-5)
            expect_equal( max(  sum( abs(unlist(out$fitted_func[[1]]) -sim$f1)),
                                sum( abs(unlist(out$fitted_func[[1]]) -sim$f2))
            )
            , 0, tol=0.2*length(sim$f1))
            expect_equal( max(  sum( abs(unlist(out$fitted_func[[2]]) -sim$f1)),
                                sum( abs(unlist(out$fitted_func[[2]]) -sim$f2))
            )
            , 0, tol=0.2*length(sim$f1))

          }
)


