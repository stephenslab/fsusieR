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
P=50
nullweight = 10/N
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

beta0       <- 0
beta1       <- 1



L=12

pos = NULL
prior =  "mixture_normal_per_scale"
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
quantile_trans=FALSE
L_start=10
gridmult =sqrt(5)
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
  resid <- 1# y - x %*% be
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

control_mixsqp=  list(verbose=FALSE,
                      eps = 1e-6,
                      numiter.em = 4
)
init_pi0_w=0.9
library(profvis)
 profvis({
   if( prior %!in% c("normal", "mixture_normal", "mixture_normal_per_scale"))
   {
     stop("Error: provide valid prior input")
   }
   if(missing(nullweight))
   {
     nullweight <- 10#/(sqrt(nrow(X)))
   }
   if(!cal_obj){
     tol <-10^-3
   }
   if(L_start >L)
   {
     L_start <- L
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


   map_data <- remap_data(Y=Y,
                          pos=pos,
                          verbose=verbose)

   outing_grid <- map_data$outing_grid
   Y           <- map_data$Y
   rm( map_data)
   # centering and scaling covariate
   X <- colScale(X)
   # centering input
   Y <- colScale(Y, scale=FALSE)
   W <- DWT2(Y)
   Y_f      <-  cbind( W$D,W$C)

   if(verbose){
     print("Starting initialization")
   }


   #X <- matrix(X)
   ### Definition of some static parameters ---

   indx_lst <-  gen_wavelet_indx(log2(length( outing_grid)))
   #removing wc with variance 0 or below a certain level

   lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)
   if(verbose){
     print( paste("Discarding ", length(lowc_wc), "wavelet coefficients out of ", ncol(Y_f)))
   }
   if(length(lowc_wc)> (ncol(Y_f )-3)){
     stop("almost all the wavelet coefficients are null/low variance, consider using univariate fine mapping")
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
   ### Definition of some dynamic parameters ------

   update_Y    <- cbind( W$D,W$C) #Using a column like phenotype, temporary matrix that will be regularly updated
   temp        <- init_prior(Y              = update_Y,
                             X              = X,
                             prior          = prior ,
                             v1             = v1,
                             indx_lst       = indx_lst,
                             lowc_wc        = lowc_wc,
                             control_mixsqp = control_mixsqp,
                             nullweight     = nullweight, max_SNP_EM=100,
                             gridmult       = gridmult )
   G_prior     <- temp$G_prior
   tt          <- temp$tt

   #Recycled for the first step of the while loop
   susiF.obj   <-  init_susiF_obj(L_max=L,
                                  G_prior=G_prior,
                                  Y=Y,
                                  X=X,
                                  L_start=L_start,
                                  greedy=greedy,
                                  backfit=backfit)
   if(verbose){
     print("Initialization done")
   }
   # numerical value to check breaking condition of while



   susiF.obj     <- susiF.workhorse(susiF.obj      = susiF.obj,
                                    W              = W,
                                    X              = X,
                                    tol            = tol,
                                    low_wc         = low_wc,
                                    init_pi0_w     = init_pi0_w ,
                                    control_mixsqp = control_mixsqp ,
                                    indx_lst       = indx_lst,
                                    lowc_wc        = lowc_wc,
                                    nullweight     = nullweight,
                                    cal_obj        = cal_obj,
                                    verbose        = verbose,
                                    cov_lev        = cov_lev,
                                    min.purity     = min.purity,
                                    maxit          = maxit,
                                    max_SNP_EM=100,
                                    tt             = tt)

 })
