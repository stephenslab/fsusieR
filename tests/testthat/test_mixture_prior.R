library(testthat)
library(ashr)
library(wavethresh)
library(mixsqp)
library(fsusieR)
set.seed(2)
f1 <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay =1.5)
lowc_wc=NULL
plot(f1$sim_func, type="l", ylab="y")
N=500
P=10
set.seed(23)
nullweight=0.1
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1 <- 5
greedy=TRUE
backfit=TRUE
noisy.data  <- list()
control_mixsqp= list(verbose=FALSE)
rsnr=1
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
tt <-  cal_Bhat_Shat(Y_f,X,v1,resid_var=min(apply(Y_f,2,mad)),
                       lowc_wc=NULL)
indx_lst <- gen_wavelet_indx(9)
Bhat <- tt$Bhat
Shat <- tt$Shat
G_prior<- init_prior(Y=Y_f,
                     X=X,
                     prior="mixture_normal",
                     v1=v1,
                     indx_lst = indx_lst,
                     lowc_wc=NULL,
                     control_mixsqp = control_mixsqp,
                     nullweight     = nullweight,
                     max_SNP_EM=100 )$G_prior

lBF <- log_BF (G_prior, Bhat, Shat , indx_lst=indx_lst,lowc_wc=NULL)
lBF
tt <-  cal_Bhat_Shat(Y_f,X,v1,
                     lowc_wc=NULL)
Bhat <- tt$Bhat
Shat <- tt$Shat
init_pi0_w= 1

control_mixsqp =  list(verbose=FALSE)
test_that("Class of the prior is", {

  expect_equal(class(
    init_prior(Y=Y_f,
               X=X,
               prior="mixture_normal",
               v1=v1,
               indx_lst = indx_lst,
                lowc_wc=NULL,
               control_mixsqp = control_mixsqp,
               nullweight     = nullweight,
               max_SNP_EM=100 )
                    $G_prior),
              "mixture_normal"
               )
})

### Test validity normal mixture  -----
Bhat <- tt$Bhat
Shat <- tt$Shat
G_prior<- init_prior.default(Y=Y_f,
                X=X,
                prior="mixture_normal",
                v1=v1,
                indx_lst = indx_lst,
                lowc_wc=NULL,
                control_mixsqp = control_mixsqp,
                nullweight     = nullweight,
                max_SNP_EM=100 )$G_prior


plot( Bhat,post_mat_mean(G_prior,Bhat,Shat,
                         indx_lst=indx_lst,
                         lowc_wc=lowc_wc))



plot( Shat,  (post_mat_sd(G_prior,Bhat,Shat, indx_lst=indx_lst,
                          lowc_wc=lowc_wc) ))
get_pi_G_prior(G_prior)
get_sd_G_prior(G_prior)
L <- L_mixsq(G_prior, Bhat, Shat)

zeta <- cal_zeta(lBF)
tpi <- m_step(L, zeta , indx_lst,init_pi0_w, control_mixsqp ,
              nullweight = nullweight)
tpi
G_update <- update_prior (G_prior, tpi)
test_that("Updated mixture proportion should be equal to provided input",
          {
            expect_equal(identical(get_pi_G_prior(G_update) ,tpi),
                         TRUE
            )
          }
)

EM_pi(G_prior,Bhat,Shat, indx_lst=indx_lst,
      init_pi0_w    = init_pi0_w,
      control_mixsqp = control_mixsqp,
      lowc_wc=lowc_wc,
      nullweight = nullweight)



test_that("Max lBF should be in postion",
          {
            lBF <- log_BF (G_prior, Bhat, Shat ,
                           lowc_wc=NULL)
            expect_equal(which.max(lBF),
                         pos1
            )
          }
)


susiF_obj <- init_susiF_obj(L_max =1,
                            G_prior,
                            Y,
                            X,
                            L_start   =1,
                            greedy = greedy,
                            backfit=backfit
)

susiF_obj$csd_X
test_that("Susif object pi are expected to be equal to ",
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


test_that("Susif object pi are expected to be equal to ",
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


test_that("Susif internal prior to be equal to ",
          {
            susiF_obj <- init_susiF_obj(L_max =1,
                                        G_prior,
                                        Y,
                                        X,
                                        L_start   =1,
                                        greedy = greedy,
                                        backfit=backfit
            )

            expect_equal(get_G_prior (susiF_obj ),  G_prior)

          }
)

test_that("Class of the prior is", {

  expect_equal(class(
    init_prior(Y=Y_f,
               X=X,
               prior="mixture_normal",
               v1=v1,
               indx_lst = indx_lst,
               lowc_wc=NULL,
               control_mixsqp = control_mixsqp,
               nullweight     = nullweight,
               max_SNP_EM=100 )$G_prior
  ),
  "mixture_normal"
  )
})
plot( Bhat,  post_mat_mean(G_prior,Bhat,Shat, indx_lst=indx_lst,lowc_w= lowc_wc) )
plot( Shat,  (post_mat_sd(G_prior,Bhat,Shat, indx_lst=indx_lst,lowc_w= lowc_wc) ))
test_that("Class of the prior is", {

  expect_equal(class(get_pi_G_prior(G_prior))
               ,
               "pi_mixture_normal"
  )
})
test_that("Class of the standard deviations  is", {

  expect_equal(class(get_sd_G_prior(G_prior))
               ,
               "sd_mixture_normal"
  )
})
get_pi_G_prior(G_prior)
get_sd_G_prior(G_prior)


L <- L_mixsq(G_prior, Bhat, Shat, indx_lst)
test_that("The likelihood computed by L_mixsqp should be of class", {
  L <- L_mixsq(G_prior, Bhat, Shat, indx_lst)

  expect_equal(class(L), "lik_mixture_normal"
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
tpi

tpi <- m_step(L, zeta , indx_lst,
              init_pi0_w    = init_pi0_w,
              control_mixsqp = control_mixsqp,
              nullweight = 1)
tpi


test_that("The output of the m_step function should of the class", {
  tpi <- m_step(L, zeta , indx_lst,
                init_pi0_w    = init_pi0_w,
                control_mixsqp = control_mixsqp,
                nullweight = nullweight)

  expect_equal( class(tpi),"pi_mixture_normal"
  )
})

test_that("The estimated null proportion should greater or equal to", {
  tpi <- m_step(L, zeta , indx_lst,
                init_pi0_w    = 1,
                control_mixsqp = control_mixsqp,
                nullweight = nullweight)
  tol <- 0.01 #tolerance
  expect_gt( get_pi0(tpi = tpi), (1-1/2^9 ) -tol  )
})



G_update <- update_prior (G_prior, tpi)
test_that("Updated mixture proportion should be equal to provided input",
          {
            expect_equal(identical(get_pi_G_prior(G_update) ,tpi),
                         TRUE
            )
          }
)
EM_pi(G_prior,Bhat,Shat, indx_lst,
      init_pi0_w    = init_pi0_w,
      control_mixsqp = control_mixsqp,
      lowc_wc=NULL,
      nullweight = nullweight)


outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                init_pi0_w    = init_pi0_w,
                control_mixsqp = control_mixsqp,
                lowc_wc=NULL,
                nullweight = nullweight)
test_that("The outputs of the EM_pi function should be  ",
          {
            outEM <-  EM_pi(G_prior,Bhat,Shat, indx_lst,
                            init_pi0_w    = init_pi0_w,
                            control_mixsqp = control_mixsqp,
                            lowc_wc=NULL,
                            nullweight = nullweight)
            expect_equal(class(outEM$tpi_k) ,"pi_mixture_normal")
            expect_equal(class(outEM$lBF) ,"numeric")
            expect_equal(length(outEM$lBF) ,dim(Bhat)[1])
            tol <- 0.01 #tolerance
            expect_gt( get_pi0(tpi = tpi), (1-1/2^9 ) -tol  )
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
            susiF_obj <- update_alpha.susiF(susiF_obj, 1, new_alpha)
            expect_equal( get_alpha (susiF_obj , 1), new_alpha )
          }
)

test_that("The mixture proportion of the update susiF object shoudl be equal to   ",
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
            susiF_obj <- update_alpha.susiF(susiF_obj, 1, new_alpha)
            expect_equal( get_alpha (susiF_obj , 1), new_alpha )
          }
)


susiF_obj <- init_susiF_obj(L_max =1,
                            G_prior,
                            Y,
                            X,
                            L_start   =1,
                            greedy = greedy,
                            backfit=backfit
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

            expect_equal( susiF_obj$fitted_wc[[1]],post_mat_mean( G_prior , Bhat, Shat, lBF= outEM$lBF, indx_lst=indx_lst,lowc_w= lowc_wc))
            expect_equal( susiF_obj$fitted_wc2[[1]],post_mat_sd  ( G_prior , Bhat, Shat ,lBF= outEM$lBF, indx_lst=indx_lst,lowc_w= lowc_wc)^2)
            expect_equal( get_alpha (susiF_obj , 1), cal_zeta(outEM$lBF))
            expect_equal( get_G_prior(susiF_obj) ,G_prior)

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
              obj = susiF_obj,
              l         = 1,
              X         = X,
              D         = W$D,
              C         = W$C,
              indx_lst  = indx_lst
            )


            id_L <-1
            update_D  <-  W$D - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF_obj$alpha[[l]], rep.int(N,P))) %*% (susiF_obj$fitted_wc[[l]][,-dim(susiF_obj$fitted_wc[[l]])[2]])  ) )
            update_C  <-  W$C - Reduce("+", lapply  ( id_L, function(l) (X*rep(susiF_obj$alpha[[l]], rep.int(N,P))) %*% susiF_obj$fitted_wc[[l]][,dim(susiF_obj$fitted_wc[[l]])[2]] ) )
            manual_update <- cbind(  update_D, update_C)
            expect_equal(  update_T ,manual_update)

          }
)


susiF_obj <- init_susiF_obj(L_max =1,
                            G_prior,
                            Y,
                            X,
                            L_start   =1,
                            greedy = greedy,
                            backfit=backfit
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
            expect_equal(  sum( abs(unlist(update_cal_fit_func(susiF_obj, indx_lst=indx_lst,TI=FALSE)$fitted_func) -f1$sim_func)), 0, tol=0.03)

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
susiF_obj <-  out_prep(susiF_obj,Y=Y,
                       X=X,
                       indx_lst=indx_lst,
                       filter_cs = FALSE,
                       TI=TRUE,
                       outing_grid = 1:ncol(Y), lfsr_curve = 0.05 )

plot( unlist(susiF_obj$fitted_func) , type="l", col="green")
lines( susiF_obj$cred_band [[1]][1,])
lines( susiF_obj$cred_band [[1]][2,])

lines(f1$sim_func, col="red")




test_that("SusiF performance should be",
          {
            set.seed(1)
            sim  <- simu_test_function(rsnr=0.5,is.plot = FALSE)
            Y <- sim$noisy.data
            X <- sim$G
            out <- susiF(Y,X,L=1, prior="mixture_normal")
            expect_equal(  unlist( out$alpha) , c(1, rep(0,9)) , tol=1e-5)
            expect_equal(  sum( abs(unlist(out$fitted_func) -sim$f1)), 0, tol=0.2*length(sim$f1))

          }
)



test_that("SusiF performance should be",
          {
            set.seed(1)
            sim  <- simu_test_function(N=500,rsnr=0.5,pos2= 2 ,is.plot = FALSE)
            Y <- sim$noisy.data
            X <- sim$G
            out <- susiF(Y,X,L=2, prior="mixture_normal")
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


