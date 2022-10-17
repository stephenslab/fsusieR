library(testthat)
library(ashr)
library(wavethresh)
library(mixsqp)
library(susiF.alpha)
control_mixsqp= list(verbose=FALSE)
set.seed(2)
f1 <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay =1.5)
lowc_wc=NULL
plot(f1$sim_func, type="l", ylab="y")
N=500
P=10
nullweight = 10/sqrt(N)
set.seed(23)
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1 <- 5
noisy.data  <- list()
greedy=TRUE
backfit=TRUE
rsnr=10
for ( i in 1:N)
{
  f1_obs <- f1$sim_func
  noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs +  rnorm(length(f1$sim_func), sd=  (1/  rsnr ) *sd(f1$sim_func))

}
noisy.data <- do.call(rbind, noisy.data)


Y <- noisy.data
X <- G
Y <- colScale(Y, scale=FALSE)
X <- colScale(X)
W <- DWT2(Y)
update_D <- W
Y_f <- cbind( W$D,W$C) #Using a column like phenotype
update_Y <-Y_f
v1 <- rep(1, dim(X)[2])
tt <-  cal_Bhat_Shat(Y_f,X,v1,lowc_wc = NULL)
indx_lst <- gen_wavelet_indx(9)
Bhat <- tt$Bhat
Shat <- tt$Shat
init_pi0_w=1
control_mixsqp =  list(verbose=FALSE)

### Test validity normal mixture per scale -----
G_prior <- init_prior(Y=Y_f,
                      X=X,
                      prior="mixture_normal_per_scale",
                      v1=v1,
                      indx_lst = indx_lst,
                      lowc_wc=NULL,
                      control_mixsqp = control_mixsqp,
                      nullweight     = nullweight )$G_prior

lBF <- log_BF (G_prior, tt$Bhat, tt$Shat , indx_lst,
               lowc_wc=NULL)

lBF
test_that("Max lBF should be in postion",
          {
            expect_equal(which.max(lBF),
                         pos1
            )
          }
)
susiF_obj <- init_susiF_obj(L_max =2,
                            G_prior,
                            Y,
                            X,
                            L_start   =2,
                            greedy = greedy,
                            backfit=backfit
                            )

susiF_obj$G_prior
test_that("susiF object pi are expected to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L_max =2,
                                        G_prior,
                                        Y,
                                        X,
                                        L_start   =2,
                                        greedy = greedy,
                                        backfit=backfit
            )

            expect_equal(get_pi(susiF_obj,1), get_pi_G_prior(G_prior)
            )
            expect_equal(get_pi(susiF_obj,2), get_pi_G_prior(G_prior)
            )
          }
)

susiF_obj <- init_susiF_obj(L_max =1,
                            G_prior,
                            Y,
                            X,
                            L_start   =4,
                            greedy = greedy,
                            backfit=backfit
)


test_that("susiF object pi are expected to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L_max =1,
                                        G_prior,
                                        Y,
                                        X,
                                        L_start   =1,
                                        greedy = greedy,
                                        backfit=backfit
            )


            expect_equal(get_pi(susiF_obj,1), get_pi_G_prior(G_prior)
            )
          }
)


test_that("correct expansion of susiF object",
{
  susiF.obj   <-  init_susiF_obj(L_max=10, G_prior, Y,X, L_start=3,
                                 greedy = greedy,
                                 backfit=backfit)
  expect_equal(susiF.obj$L_max, 10  )
  expect_equal(susiF.obj$L, 3  )
  susiF.obj <- expand_susiF_obj(susiF.obj,L_extra=7)
  expect_equal(susiF.obj$L_max, 10  )
  expect_equal(susiF.obj$L,length(susiF.obj$fitted_wc))
  expect_equal(susiF.obj$L,length(susiF.obj$alpha))
  expect_equal(susiF.obj$L,length(susiF.obj$fitted_wc2))
  expect_equal(susiF.obj$L,length(susiF.obj$G_prior))
  expect_equal(susiF.obj$L,length(susiF.obj$cs))
  expect_equal(susiF.obj$L,length(susiF.obj$est_pi))
  expect_equal(susiF.obj$L,length(susiF.obj$est_sd))
  expect_equal(susiF.obj$L,length(susiF.obj$lBF))
  expect_equal(susiF.obj$L,length(susiF.obj$cred_band))

}
)
test_that("susiF internal prior to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L_max=2, G_prior, Y,X, L_start=2,
                                        greedy = greedy,
                                        backfit=backfit)

            expect_equal(get_G_prior (susiF_obj ),  G_prior)

          }
)

test_that("Class of the prior is", {

  expect_equal(class(
    init_prior(Y=Y_f,
               X=X,
               prior="mixture_normal_per_scale",
               v1=v1,
               indx_lst = indx_lst,
               lowc_wc=NULL,
               control_mixsqp = control_mixsqp,
               nullweight     = nullweight )$G_prior

  )[1],
  "mixture_normal_per_scale"
  )
})
plot( Bhat,  post_mat_mean(G_prior,Bhat,Shat, indx_lst,lowc_wc) )
plot( Shat,  (post_mat_sd(G_prior,Bhat,Shat, indx_lst,lowc_wc) ))
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



L <- L_mixsq(G_prior, Bhat, Shat, indx_lst )


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

tpi <- m_step(L, zeta , indx_lst,
              init_pi0_w    = init_pi0_w,
              control_mixsqp = control_mixsqp,
              nullweight = nullweight)
class(tpi)

test_that("The output of the m_step function should of the class", {
  tpi <- m_step(L, zeta , indx_lst,
                init_pi0_w    = init_pi0_w,
                control_mixsqp = control_mixsqp,
                nullweight = nullweight)

  expect_equal( class(tpi)[1],"pi_mixture_normal_per_scale"
  )
})

test_that("The output of the m_step for the pi_0 should equal", {
  tpi <- m_step(L, zeta , indx_lst,
                init_pi0_w    = init_pi0_w,
                control_mixsqp = control_mixsqp,
                nullweight = nullweight)
  expect_equal( get_pi0(tpi = tpi), c(0,0.5, rep(1, 8)),
                tolerance = 0.02) #allow 1% error in the proportion estimation
})


G_update <- update_prior (G_prior, tpi)
test_that("Updated mixture proportion should be equal to provided input",
          {
            expect_equal(identical(get_pi_G_prior(G_update) ,tpi),
                         TRUE
            )
          }
)
outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst, control_mixsqp = control_mixsqp,
                lowc_wc=NULL,
                nullweight = nullweight)

test_that("The outputs of the EM_pi function should be  ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)

            expect_equal(class(outEM$tpi_k)[1] ,"pi_mixture_normal_per_scale")
            expect_equal(class(outEM$lBF)[1] ,"numeric")
            expect_equal(length(outEM$lBF)[1] ,dim(Bhat)[1])
            expect_equal( get_pi0(tpi = outEM$tpi_k), c(0,0.5, rep(1, 8)),
                          tolerance = 0.01)
          }
)

test_that("The mixture proportion of the updated susiF object should be equal to   ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)

            susiF_obj <- update_pi( susiF_obj, 1,  outEM$tpi_k)

            expect_equal( get_pi (susiF_obj , 1),outEM$tpi_k )
          }
)


test_that("The alpha value of  the update susiF object should be equal to   ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)

            new_alpha <- cal_zeta(outEM$lBF)
            susiF_obj <- update_alpha(susiF_obj, 1, new_alpha)
            expect_equal( get_alpha (susiF_obj , 1), new_alpha )
          }
)

test_that("The update susiF object should have its argument equal to    ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)

            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )

            expect_equal( susiF_obj$fitted_wc[[1]],post_mat_mean( G_prior , Bhat, Shat, indx_lst ,lowc_wc))
            expect_equal( susiF_obj$fitted_wc2[[1]],post_mat_sd  ( G_prior , Bhat, Shat , indx_lst,lowc_wc)^2)
            expect_equal( get_alpha (susiF_obj , 1), cal_zeta(outEM$lBF))
            expect_equal( get_G_prior(susiF_obj) ,G_prior)
            expect_equal(   susiF_obj$lBF[[1]]  , outEM$lBF )
          }
)

test_that("The partial residual should be    ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)

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



test_that("The precision of the fitted curves should be   ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)

            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )
            expect_equal(  sum( abs(unlist(update_cal_fit_func(susiF_obj, indx_lst)$fitted_funcfitted_func[[1]]) -f1$sim_func)), 0, tol=0.01)


          }
)

outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                init_pi0_w    = init_pi0_w,
                control_mixsqp = control_mixsqp,
                lowc_wc=NULL,
                nullweight = nullweight)

G_prior <- update_prior(G_prior,
                        tpi= outEM$tpi_k )

susiF_obj <-  update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst ,
                               lowc_wc=NULL)

susiF_obj <-  out_prep(susiF_obj,Y, X=X, indx_lst=indx_lst,filter.cs = FALSE)

plot( (unlist(susiF_obj$fitted_func[[1]])), type="l", col="green")
lines(f1$sim_func, col="red")



test_that("SusiF performance should be",
          {
            set.seed(1)
            sim  <- simu_test_function(rsnr=2,is.plot = FALSE)
            Y <- sim$noisy.data
            X <- sim$G
            out <- susiF(Y,X,L=1, prior="mixture_normal_per_scale", nullweight = .0)
            expect_equal(  unlist( out$alpha) , c(1, rep(0,9)) , tol=1e-5)
            expect_equal(  sum( abs(unlist(out$fitted_func) -sim$f1)), 0, tol=0.2*length(sim$f1))

          }
)



test_that("The expected sum of square should be below and residual variance should be ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)

            G_prior <- update_prior(G_prior,
                                    tpi= outEM$tpi_k )

            susiF_obj <- update_susiF_obj(susiF_obj, 1, outEM, Bhat, Shat, indx_lst )
            sigma2 <- estimate_residual_variance(susiF_obj,Y_f,X)
            susiF_obj <- update_residual_variance(susiF_obj, sigma2 = sigma2 )
            expect_equal(  susiF_obj$sigma2, sigma2)


          }
)


test_that("The KL of effect one ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)

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
            sim  <- simu_test_function(rsnr=1,pos2= 2 ,is.plot = FALSE)
            Y <- sim$noisy.data
            X <- sim$G
            out <- susiF(Y,X,L=2, prior="mixture_normal_per_scale" ,nullweight = 0, init_pi0_w = 1 )
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




test_that("Removing one wc coeef should lead to the followin results",
          {

            tt1 <- cal_Bhat_Shat (update_Y,X,v1 , lowc_wc=NULL  )
            tt2 <- cal_Bhat_Shat (update_Y,X,v1 , lowc_wc=1:10  )
            expect_equal(c(tt1$Bhat[,-c(1:10)] ),c(tt2$Bhat[,-c(1:10)]))
            expect_equal(c(tt1$Shat[,-c(1:10)] ),c(tt2$Shat[,-c(1:10)]))
            expect_equal(c(tt2$Bhat[, c(1:10)] ),rep(0, length(c(tt2$Bhat[, c(1:10)] ))))
            expect_equal(c(tt2$Shat[, c(1:10)] ),rep(1, length(c(tt2$Shat[, c(1:10)] ))))

          }
)


#out <- susiF(Y,X,L=2, prior="mixture_normal_per_scale", cal_obj = TRUE)
#out <- susiF(Y,X,L=2, prior="mixture_normal_per_scale", cal_obj = TRUE,quantile_trans = TRUE)

out <- susiF(Y,X,L=2, prior="mixture_normal_per_scale", cal_obj = FALSE,quantile_trans = TRUE, filter.cs = FALSE)
out$cs

