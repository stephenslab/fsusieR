library(profvis)
library(ashr)
library(susiF.alpha)
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
rsnr <- 0.2 #wished root signal noise ratio
N <- 100 #Number of individuals
P <- 10 # Number of covariates
pos1 <- 1#Position of the causal covariate
lev_res <- 7
temp_func <-  simu_IBSS_per_level(lev_res )
f1 <-  temp_func$sim_func
plot( f1, type ="l")
G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1

noisy.data  <- list()

for ( i in 1:N)
{
  f1_obs <- f1
  noise <- rnorm(length(f1), sd=  (1/  rsnr ) * var(f1))
  noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs +  noise

}
noisy.data <- do.call(rbind, noisy.data)

  L =2
pos = NULL
prior = "mixture_normal_per_scale"
verbose = TRUE
plot_out = TRUE
maxit = 100
tol = 1e-3
cov_lev = 0.95
min.purity=0.5
lfsr_curve = 0.05
filter.cs =TRUE

#testing if x is a wholenumber
#'
is.wholenumber <- function (x, tol = .Machine$double.eps^0.5)
  abs(x - round(x)) < tol

# Not in operator

'%!in%' <- function(x,y)!('%in%'(x,y))

#based on Rfast implementation
#'
fast_lm <- function(x,y)
{
  be <- solve(crossprod(x),crossprod(x,y))
  resid <-  y - x %*% be
  out <- list(be = be,
              residuals = resid)
  return(out)
}


#Circular permutation on vector
# Code adapted from https://mzuer.github.io
#'
shifter <- function(x, n = 1) {
  # if (n == 0) x else c(tail(x, -n), head(x, n))
  if (n == 0) x else c(tail(x, n), head(x, -n))
}

update_residual_variance  <- function(susiF.obj,sigma2, ...)
  UseMethod("update_residual_variance")

#' @rdname update_residual_variance
#'
#' @method update_residual_variance susiF
#'
#' @export update_residual_variance.susiF
#'
#' @export
#'

update_residual_variance.susiF <- function(susiF.obj,sigma2)
{
  susiF.obj$sigma2 <- sigma2
  return(susiF.obj)
}



Y <- noisy.data
X <- G

 profvis({
   if( prior %!in% c("normal", "mixture_normal", "mixture_normal_per_scale"))
   {
     stop("Error: provide valid prior input")
   }

   ## Input error messages

   if (is.null(pos))
   {
     pos <- 1:dim(Y)[2]
   }


   #reshaping of the data
   if ( !(length(pos)==dim(Y)[2])) #miss matching positions and number of observations
   {
     stop("Error: number of position provided different from the number of column of Y")
   }
   original_Y <-Y


   if(!is.wholenumber(log2(dim(Y)[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check wether dim(Y) not equal to 2^J or if the data are unevenly spaced
   {

     inter_pol.obj <- interpol_mat(Y, pos)
     Y             <- inter_pol.obj$Y
     bp            <- inter_pol.obj$bp
     outing_grid   <- inter_pol.obj$grid
     if(verbose)
     {
       message( "Response matrix dimensions not equal to nx 2^J \n or unevenly spaced data \n interpolation procedure used")
     }
   }  else{

     outing_grid <- 1:dim(Y)[2]
   }
   W <- DWT2(Y)

   ### Definition of some static parameters ---
   indx_lst <-  gen_wavelet_indx(log2(length( outing_grid)))
   v1       <-  rep(1, dim(X)[1])### used in fit_lm to add a column of 1 in the design matrix
   Y_f      <-  cbind( W$D,W$C)
   # Wavelet transform of the inputs


   update_D <- W
   ### Definition of some dynamic parameters ---

   update_Y    <-  cbind( W$D,W$C) #Using a column like phenotype, temporary matrix that will be regularly updated
   G_prior     <-  init_prior(update_Y,X,prior,v1 , indx_lst  )
   susiF.obj   <-  init_susiF_obj(L=L, G_prior, Y,X)

   # numerical value to check breaking condition of while
   check <- 1
   h     <- 0

   if(susiF.obj$L==1)
   {
     tt   <- cal_Bhat_Shat(update_Y,X,v1)
     Bhat <- tt$Bhat
     Shat <- tt$Shat #UPDATE. could be nicer
     tpi  <- get_pi(susiF.obj,1)
     G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

     EM_out  <- EM_pi(G_prior  = G_prior,
                      Bhat     = Bhat,
                      Shat     = Shat,
                      indx_lst = indx_lst
     )

     susiF.obj <-  update_susiF_obj(susiF.obj = susiF.obj ,
                                    l         = 1,
                                    EM_pi     = EM_out,
                                    Bhat      = Bhat,
                                    Shat      = Shat,
                                    indx_lst  = indx_lst
     )
     susiF.obj <- update_ELBO(susiF.obj,
                              get_objective( susiF.obj = susiF.obj,
                                             Y         = Y_f,
                                             X         = X,
                                             D         = W$D,
                                             C         = W$C,
                                             indx_lst  = indx_lst
                              )
     )

   }else{
     while(check >tol & (h/L) <maxit)
     {
       for( l in 1:susiF.obj$L)
       {

         h <- h+1
         tt <- cal_Bhat_Shat(update_Y,X,v1)
         Bhat <- tt$Bhat
         Shat <- tt$Shat #UPDATE. could be nicer
         tpi <-  get_pi(susiF.obj,l)
         G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

         EM_out  <- EM_pi(G_prior  = G_prior,
                          Bhat     =  Bhat,
                          Shat     =  Shat,
                          indx_lst =  indx_lst
         )

         susiF.obj <-  update_susiF_obj(susiF.obj = susiF.obj ,
                                        l         = l,
                                        EM_pi     = EM_out,
                                        Bhat      = Bhat,
                                        Shat      = Shat,
                                        indx_lst  = indx_lst
         )

         update_Y  <-  cal_partial_resid(
           susiF.obj = susiF.obj,
           l         = l,
           X         = X,
           D         = W$D,
           C         = W$C,
           indx_lst  = indx_lst
         )


       }#end for l in 1:L


       susiF.obj <- update_ELBO(susiF.obj,
                                get_objective( susiF.obj = susiF.obj,
                                               Y         = Y_f,
                                               X         = X,
                                               D         = W$D,
                                               C         = W$C,
                                               indx_lst  = indx_lst
                                )
       )

       sigma2    <- estimate_residual_variance(susiF.obj,Y=Y_f,X)
       susiF.obj <- update_residual_variance(susiF.obj, sigma2 = sigma2 )

       if(length(susiF.obj$ELBO)>1 )#update parameter convergence,
       {
         check <- diff(susiF.obj$ELBO)[(length( susiF.obj$ELBO )-1)]

       }
     }#end while
   }


   #preparing output
   susiF.obj <- out_prep(susiF.obj  = susiF.obj,
                         Y          = Y,
                         X          = X,
                         indx_lst   = indx_lst,
                         filter.cs  = filter.cs,
                         lfsr_curve = lfsr_curve
   )
 })
