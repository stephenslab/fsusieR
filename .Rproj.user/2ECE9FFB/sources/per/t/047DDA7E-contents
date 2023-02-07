library(susiF.alpha)

verbose=TRUE

gen_wavelet_indx(7)
#Example using curves simulated under the Mixture normal per scale prior
rsnr <- .5 #expected root signal noise ratio
N <- 1000    #Number of individuals
P <- 3     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 2   #Position of the causal covariate for effect 2
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
G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
beta2       <- 1
noisy.data  <- list()

for ( i in 1:N)
{
  f1_obs <- f1
  f2_obs <- f2
  noise <- rnorm(length(f1), sd=  (1/  rsnr ) * var(f1))
  noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs + beta2*G[i,pos2]*f2_obs + noise

}
noisy.data <- do.call(rbind, noisy.data)




plot( noisy.data[1,], type = "l", col=(G[1, pos1]*3+1),
      main="Observed curves \n colored by the causal effect", ylim= c(-40,40), xlab="")
for ( i in 2:N)
{
  lines( noisy.data[i,], type = "l", col=(G[i, pos1]*3+1))

}
legend(x=0.3,
       y=-10,
       lty = rep(1,3),
       legend= c("0", "1","2"),
       col=c("black","blue","yellow"))



Y <- noisy.data
X <- G
pt <- proc.time()
 prior ="mixture_normal_per_scale"

  nullweight <- 10#/(sqrt(nrow(X)))

  tol <-10^-3
  control_mixsqp =  list(verbose=FALSE,
                         eps = 1e-6,
                         numiter.em = 4
  )
  init_pi0_w=0.8
## Input error messages

  pos <- 1:dim(Y)[2]
  thresh_lowcount=0

#reshaping of the data




  outing_grid   <- pos
# centering and scaling covariate
X <- colScale(X)
# centering input
Y <- colScale(Y, scale=FALSE)
W <- DWT2(Y)
Y_f      <-  cbind( W$D,W$C)


#X <- matrix(X)
### Definition of some static parameters ---

indx_lst <-  gen_wavelet_indx(log2(length( outing_grid)))
#removing wc with variance 0 or below a certain level

lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)
if(verbose){
  print( paste("Discarding ", length(lowc_wc), "wavelet coefficients out of ", ncol(Y_f)))
}


v1       <-  rep(1, dim(X)[1])### used in fit_lm to add a column of 1 in the design matrix
# Wavelet transform of the inputs


update_D <- W
### Definition of some dynamic parameters ------

update_Y    <- cbind( W$D,W$C) #Using a column like phenotype, temporary matrix that will be regularly updated
temp        <- init_prior(Y              = update_Y,
                          X              = X,
                          prior          = prior ,
                          v1             = v1,
                          indx_lst       = indx_lst,
                          lowc_wc        = lowc_wc,
                          control_mixsqp = control_mixsqp,
                          nullweight     = nullweight,
                          gridmult       = sqrt(2))
G_prior     <- temp$G_prior
tt          <- temp$tt
init        <- TRUE

#Recycled for the first step of the while loop
EBmvFR.obj   <-  init_EBmvFR_obj(G_prior=G_prior,
                                 Y=Y,
                                 X=X)
lowc_wc=NULL
for(j in 1:3){

  EBmvFR.obj   <-  fit_effect.EBmvFR (EBmvFR.obj = EBmvFR.obj,
                                      j         = j,
                                      X         = X,
                                      D         = W$D,
                                      C         = W$C,
                                      indx_lst  = indx_lst,
                                      lowc_wc=NULL)
}


sigma2    <- estimate_residual_variance(EBmvFR.obj,Y=Y_f,X)

print(sigma2)
EBmvFR.obj <- update_residual_variance(EBmvFR.obj, sigma2 = sigma2 )

EBmvFR.obj <-  update_prior.EBmvFR ( EBmvFR.obj,
                                    max_step = 100,
                                    espsilon = 0.0001,
                                    init_pi0_w =1,
                                    control_mixsqp,
                                    lowc_wc,
                                    nullweight)


get_pi_G_prior(get_G_prior(EBmvFR.obj))



cal_obj=FALSE

maxit=10
check=10
tol =1e-5
iter =1
while( (check >tol & iter <maxit))
{
  for( j in 1:ncol(X))
  {

    if(verbose){
      print(paste("Fitting effect ", j,", iter" ,  iter ))
    }
    EBmvFR.obj   <-  fit_effect.EBmvFR (EBmvFR.obj = EBmvFR.obj,
                                        j         = j,
                                        X         = X,
                                        D         = W$D,
                                        C         = W$C,
                                        indx_lst  = indx_lst)

  }#end for l in 1:L  -----

  sigma2    <- estimate_residual_variance(EBmvFR.obj,Y=Y_f,X)
  print(sigma2)
  EBmvFR.obj <- update_residual_variance(EBmvFR.obj, sigma2 = sigma2 )

  EBmvFR.obj <- update_prior( EBmvFR.obj,
                              max_step       = 100,
                              espsilon       = 0.0001,
                              init_pi0_w     = init_pi0_w ,
                              control_mixsqp = control_mixsqp,
                              lowc_wc        = low_wc,
                              nullweight     = nullweight)

  EBmvFR.obj <- test_stop_cond(EBmvFR.obj = EBmvFR.obj,
                               check      = check,
                               cal_obj    = cal_obj,
                               Y          = Y_f,
                               X          = X,
                               D          = W$D,
                               C          = W$C,
                               indx_lst    = indx_lst)

  #print(EBmvFR.obj$alpha)
  #print(EBmvFR.obj$ELBO)
  check <- EBmvFR.obj$check



  iter <- iter +1


}#end while
EBmvFR.obj <-  update_cal_fit_func(EBmvFR.obj, indx_lst)
plot( EBmvFR.obj$fitted_func[1,],type="l", col="blue" )#fitted
lines( f1, col="green")#true
plot( EBmvFR.obj$fitted_func[2,],type="l", col="blue" )#fitted
lines( f2, col="green")#true
plot( EBmvFR.obj$fitted_func[3,],type="l", col="blue" )#fitted

res <- EBmvFR (Y=noisy.data,X=X)
plot( res$fitted_func[1,],type="l", col="blue" )#fitted
lines( f1, col="green")#true
plot( res$fitted_func[2,],type="l", col="blue" )#fitted
lines( f2, col="green")#true

sqrt(res$sigma2)
(1/  rsnr ) * var(f1)
