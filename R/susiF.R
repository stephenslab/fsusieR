#' @title Sum of Single Function
#'
#' @description Implementation of the SuSiF method
#'
#' @details Implementation of the SuSiF method
#'
#'
#' @param Y functional phenotype, matrix of size N by size J. The
#'   underlying algorithm uses wavelet, which assumes that J is of the
#'   form J^2. If J is not a power of 2, susiF internally remaps the data
#'   into a grid of length 2^J
#'
#' @param X matrix of size n by p contains the covariates
#'
#' @param L upper bound on the number of effects to fit (if not specified, set to =2)
#'
#' @param pos vector of length J, corresponding to position/time pf
#' the observed column in Y, if missing, suppose that the observation
#' are evenly spaced
#'
#' @param prior specify the prior used in susiF. The two available choices are
#' available "mixture_normal_per_scale", "mixture_normal". Default "mixture_normal_per_scale",
#' if this susiF is too slow, consider using  "mixture_normal" (up to 40% faster), but this may result in
#' oversmoothing the estimated curves.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings are printed to the
#' console.
#'
#'
#' @param tol a small, non-negative number specifying the convergence
#' tolerance for the IBSS fitting procedure. The fitting procedure
#' will halt when the difference in the variational lower bound, or
#' \dQuote{ELBO} (the objective function to be maximized), is less
#' than \code{tol}.
#'
#' @param maxit Maximum number of IBSS iterations.
#'
#' @param cov_lev numeric between 0 and 1, corresponding to the
#' expected level of coverage of the cs if not specified set to 0.95
#'
#' @param min.purity minimum purity for estimated credible sets
#' @param filter.cs logical, if TRUE filter the credible set (removing low purity
#' cs and cs with estimated prior equal to 0). Set as TRUE by default.
#' @param init_pi0_w starting value of weight on null compoenent in mixsqp
#'  (between 0 and 1)
#' @param control_mixsqp list of parameter for mixsqp function see  mixsqp package
#' @param  cal_obj logical if set as TRUE compute ELBO for convergence monitoring
#' @param quantile_trans logical if set as TRUE perform normal quantile transform
#' on wavelet coefficients
#' @param L_start number of effect initialized at the start of the algorithm
#' @param nullweight numeric value for penalizing likelihood at point mass 0 (should be between 0 and 1)
#' (usefull in small sample size)
#' @param thresh_lowcount numeric, used to check the wavelet coefficients have
#'  problematic distribution (very low dispersion even after standardization).
#'  Basically check if the median of the absolute value of the distribution of
#'   a wavelet coefficient is below this threshold. If yes, the algorithm discard
#'   this wavelet coefficient (setting its estimate effect to 0 and estimate sd to 1).
#'   Set to 0 by default. It can be useful when analyzing sparse data from sequence
#'    based assay or small samples.
#' @param greedy logical, if TRUE allows greedy search for extra effect
#'  (up to L specified by the user). Set as TRUE by default
#' @param backfit logical, if TRUE allow discarding effect via backfitting.
#'  Set as true by default as TRUE. We advise keeping it as TRUE
#' @param gridmult numeric used to control the number of components used in the mixture prior (see ashr package
#'  for more details). From the ash function:  multiplier by which the default grid values for mixsd differ from one another.
#'   (Smaller values produce finer grids.). Increasing this value may reduce computational time
#'
#'
#' @examples
#'
#'library(ashr)
#'library(wavethresh)
#'set.seed(1)
#'#Example using curves simulated under the Mixture normal per scale prior
#'rsnr <- 0.2 #expected root signal noise ratio
#'N <- 100    #Number of individuals
#'P <- 10     #Number of covariates/SNP
#'pos1 <- 1   #Position of the causal covariate for effect 1
#'pos2 <- 5   #Position of the causal covariate for effect 2
#'lev_res <- 7#length of the molecular phenotype (2^lev_res)
#'f1 <-  simu_IBSS_per_level(lev_res )$sim_func#first effect
#'f2 <- simu_IBSS_per_level(lev_res )$sim_func #second effect
#'
#'plot( f1, type ="l", ylab="effect", col="blue")
#'abline(a=0,b=0)
#'lines(f2, type="l", col="green")
#'
#'legend(x=100,
#'       y=3,
#'       lty = rep(1,3),
#'       legend= c("effect 1", "effect 2" ),
#'       col=c("black","blue","yellow"))
#'G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
#'beta0       <- 0
#'beta1       <- 1
#'beta2       <- 1
#'noisy.data  <- list()
#'
#'for ( i in 1:N)
#'{
#'  f1_obs <- f1
#'  f2_obs <- f2
#'  noise <- rnorm(length(f1), sd=  (1/  rsnr ) * var(f1))
#'  noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs + beta2*G[i,pos2]*f2_obs + noise
#'
#'}
#'noisy.data <- do.call(rbind, noisy.data)
#'
#'
#'
#'
#'plot( noisy.data[1,], type = "l", col=(G[1, pos1]*3+1),
#'      main="Observed curves \n colored by the causal effect", ylim= c(-40,40), xlab="")
#'for ( i in 2:N)
#'{
#'  lines( noisy.data[i,], type = "l", col=(G[i, pos1]*3+1))
#'
#'}
#'legend(x=0.3,
#'       y=-10,
#'       lty = rep(1,3),
#'       legend= c("0", "1","2"),
#'       col=c("black","blue","yellow"))
#'
#'
#'
#'Y <- noisy.data
#'X <- G
#'#Running fSuSiE
#'
#'out <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale')
#'#the easiest way to visualize the result is to use the plot_susiF function
#'
#'plot_susiF(out)
#'
#'#You can also access the information directly in the output of susiF  as follow
#'par(mfrow=c(1,2))
#'
#'plot( f1, type="l", main="Estimated effect 1", xlab="")
#'lines(unlist(out$fitted_func[[1]]),col='blue' )
#'abline(a=0,b=0)
#'legend(x= 35,
#'       y=3,
#'       lty= rep(1,2),
#'       legend = c("effect 1"," fSuSiE est "),
#'       col=c("black","blue" )
#')
#'plot( f2, type="l", main="Estimated effect 2", xlab="")
#'lines(unlist(out$fitted_func[[2]]),col='green' )
#'abline(a=0,b=0)
#'legend(x= 20,
#'       y=-1.5,
#'       lty= rep(1,2),
#'       legend = c("effect 2"," fSuSiE est "),
#'       col=c("black","green" )
#')
#'
#'par(mfrow=c(1,1))
#'plot_susiF(out)
#'
#' @importFrom stats var
#'
#' @export
#'
susiF <- function(Y, X, L = 2,
                  pos = NULL,
                  prior = "mixture_normal_per_scale",
                  verbose = TRUE,
                  maxit = 100,
                  tol = 1e-3,
                  cov_lev = 0.95,
                  min.purity=0.5,
                  filter.cs =TRUE,
                  init_pi0_w= 1,
                  nullweight ,
                  control_mixsqp =  list(verbose=FALSE,
                                         eps = 1e-6,
                                         numiter.em = 4
                  ),
                  thresh_lowcount=0,
                  cal_obj=FALSE,
                  L_start = 3,
                  quantile_trans=FALSE,
                  greedy =TRUE,
                  backfit =TRUE,
                  gridmult= sqrt(2)
)
{


  ####Cleaning input -----
  pt <- proc.time()
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
  if(verbose){
    print("Data transform")
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
  if(verbose){
    print("Data transform done")
  }
  ### Definition of some dynamic parameters ------
  if(verbose){
    print("Initializing prior")
  }
  update_Y    <- cbind( W$D,W$C) #Using a column like phenotype, temporary matrix that will be regularly updated
  temp        <- init_prior(Y              = update_Y,
                            X              = X,
                            prior          = prior ,
                            v1             = v1,
                            indx_lst       = indx_lst,
                            lowc_wc        = lowc_wc,
                            control_mixsqp = control_mixsqp,
                            nullweight     = nullweight,
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
                                   tt             = tt)

  #preparing output
  susiF.obj <- out_prep(susiF.obj   = susiF.obj,
                        Y           = Y,
                        X           = X,
                        indx_lst    = indx_lst,
                        filter.cs   = filter.cs,
                        outing_grid = outing_grid
  )
  susiF.obj$runtime <- proc.time()-pt
  return(susiF.obj)
}
