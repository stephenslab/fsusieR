library(testthat)
library(ashr)
library(wavethresh)
library(mixsqp)
library(susiF.alpha)
set.seed(2)
f1 <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay =1.5)
lowc_wc=NULL
plot(f1$sim_func, type="l", ylab="y")
N=500
P=10
nullweight = 10/N
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

beta0       <- 0
beta1       <- 1



L=6

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
thresh_lowcount <-0.001
cal_obj=FALSE
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

which_lowcount <- function( Y_f, thresh_lowcount ){
  tt <- which( apply( abs(Y_f),2,median) <thresh_lowcount )

  if(length(tt)==0)
  {
    return(NULL)
  }else{
    return(tt)
  }
}



 profvis({
   if( prior %!in% c("normal", "mixture_normal", "mixture_normal_per_scale"))
   {
     stop("Error: provide valid prior input")
   }

     nullweight <- 10/nrow(X)


   ## Input error messages

     pos <- 1:dim(Y)[2]


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
   Y_f      <-  cbind( W$D,W$C)

   if(verbose){
     print("Starting initialization ")
   }

   ### Definition of some static parameters ---
   indx_lst <-  gen_wavelet_indx(log2(length( outing_grid)))

     lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)
     if(verbose){
       print( paste("Discarding ", length(lowc_wc), "wavelet coefficients out of ", ncol(Y_f)))
     }
     if(length(lowc_wc)> (ncol(Y_f )-3)){
       stop("almost all the wavelet coefficients are null, consider using univariate fine mapping")
     }

   if(quantile_trans)# important to do it after testing for lowcount
   {
     W$C <- Quantile_transform(W$C )

     W$D <- apply( W$D,2,  Quantile_transform )
     Y_f      <-  cbind( W$D,W$C)
   }

   v1       <-  rep(1, dim(X)[1])### used in fit_lm to add a column of 1 in the design matrix
   # Wavelet transform of the inputs


   update_D <- W
   ### Definition of some dynamic parameters ---

   update_Y    <-  cbind( W$D,W$C) #Using a column like phenotype, temporary matrix that will be regularly updated
   G_prior     <-  init_prior(update_Y,X,prior,v1, indx_lst, lowc_wc  )
   susiF.obj   <-  init_susiF_obj(L=L, G_prior, Y,X)

   # numerical value to check breaking condition of while
   check <- 1
   h     <- 0

   if(verbose){
     print("Initialization done")
   }

   if( L==1)
   {
     tt   <- cal_Bhat_Shat(update_Y,X,v1 , lowc_wc =lowc_wc )
     Bhat <- tt$Bhat
     Shat <- tt$Shat #UPDATE. could be nicer
     tpi  <- get_pi(susiF.obj,1)
     G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

     EM_out  <- EM_pi(G_prior        = G_prior,
                      Bhat           = Bhat,
                      Shat           = Shat,
                      indx_lst       = indx_lst,
                      init_pi0_w     = init_pi0_w,
                      control_mixsqp = control_mixsqp,
                      lowc_wc        = lowc_wc,
                      nullweight     = nullweight

     )

     susiF.obj <-  update_susiF_obj(susiF.obj = susiF.obj ,
                                    l         = 1,
                                    EM_pi     = EM_out,
                                    Bhat      = Bhat,
                                    Shat      = Shat,
                                    indx_lst  = indx_lst,
                                    lowc_wc   = lowc_wc
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
       for( l in 1:L)
       {
         if(verbose){
           print(paste("Fitting effect ", l,", iter" , (ceiling(h/L)+1)))
         }
         h <- h+1
         tt <- cal_Bhat_Shat(update_Y,X,v1, lowc_wc =lowc_wc )
         Bhat <- tt$Bhat
         Shat <- tt$Shat #UPDATE. could be nicer
         tpi <-  get_pi(susiF.obj,l)
         G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

         EM_out  <- EM_pi(G_prior        = G_prior,
                          Bhat           = Bhat,
                          Shat           = Shat,
                          indx_lst       = indx_lst,
                          init_pi0_w     = init_pi0_w,
                          control_mixsqp = control_mixsqp,
                          lowc_wc        = lowc_wc,
                          nullweight     = nullweight
         )
         #print(h)
         #print(EM_out$lBF)
         susiF.obj <-  update_susiF_obj(susiF.obj   = susiF.obj ,
                                        l           = l,
                                        EM_pi       = EM_out,
                                        Bhat        = Bhat,
                                        Shat        = Shat,
                                        indx_lst    = indx_lst,
                                        lowc_wc     = lowc_wc
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


       if( cal_obj){
         susiF.obj <- update_KL(susiF.obj,
                                X,
                                D= W$D,
                                C= W$C , indx_lst)

         if( h>1){
           susiF.obj <- update_ELBO(susiF.obj,
                                    get_objective( susiF.obj = susiF.obj,
                                                   Y         = Y_f,
                                                   X         = X,
                                                   D         = W$D,
                                                   C         = W$C,
                                                   indx_lst  = indx_lst
                                    )
           )

         }
         if(length(susiF.obj$ELBO)>1    )#update parameter convergence,
         {
           check <- abs(diff(susiF.obj$ELBO)[(length( susiF.obj$ELBO )-1)])

         }
       }
       else{
         if(ceiling(h/L)>1)#update parameter convergence, no ELBO for the moment
         {
           check <-0
           for( tt in 0:(susiF.obj$L-1))
           {
             check <-  check + var( susiF.obj$alpha_hist[[h-tt]] -susiF.obj$alpha_hist [[h-L-tt]])
           }
           check <- check/nrow(X)
           #print(check)
         }
       }

       sigma2    <- estimate_residual_variance(susiF.obj,Y=Y_f,X)
       susiF.obj <- update_residual_variance(susiF.obj, sigma2 = sigma2 )




     }#end while
   }


   #preparing output
   susiF.obj <- out_prep(susiF.obj  = susiF.obj,
                         Y          = Y,
                         X          = X,
                         indx_lst   = indx_lst,
                         filter.cs  = filter.cs
   )
 })
