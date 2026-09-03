#' @title Empirical Bayes multivariate functional regression
#'
#' @description Empirical Bayes multivariate functional regression
#'
#' @details Fits a factorized variational posterior by coordinate ascent in
#' wavelet space. Normal-mixture prior variances are on the coefficient scale
#' and are not multiplied by the residual variance. Mixture weights are
#' updated by averaging posterior component responsibilities, separately at
#' each wavelet scale when using \code{mixture_normal_per_scale}.
#' A shared residual variance is estimated from the expected residual sum of
#' squares over the retained wavelet columns.
#'
#'
#' @param Y functional phenotype, matrix of size N by size J. The
#'   underlying algorithm uses wavelets and expects a power-of-two grid.
#'   Otherwise the data are interpolated to a power-of-two grid. At least
#'   two individuals and four response positions are required.
#'
#' @param X finite numeric matrix of size n by p containing the covariates.
#'   Constant columns, including an explicit intercept, must be omitted;
#'   the response and predictors are centered internally.
#' @param adjust logical if set to TRUE (default FALSE), then the output contains the adjusted coeficients (usefull to correct for batch effect)
#' @param pos vector of length J, corresponding to position/time pf
#' the observed column in Y, if missing, suppose that the observation
#' are evenly spaced
#'
#' @param prior specify the prior used in susiF. The two available choices are
#' available "mixture_normal_per_scale", "mixture_normal". Default "mixture_normal_per_scale",
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings are printed to the
#' console.
#'
#'
#' @param tol Non-negative convergence tolerance. With \code{cal_obj=TRUE},
#'   stop when the absolute ELBO change is at most \code{tol}. Otherwise,
#'   check relative changes in coefficients and residual variance, and
#'   absolute changes in mixture probabilities. An explicitly supplied
#'   tolerance is always respected. When omitted, the default is 0.1 for
#'   ELBO monitoring and 0.001 for parameter monitoring.
#'
#' @param maxit Maximum number of complete coordinate-ascent sweeps.
#'
#' @param init_pi0_w,control_mixsqp,nullweight,max_step_EM,max_SNP_EM Legacy
#'   optimizer controls retained for call compatibility. These do not control
#'   the coordinate-ascent updates, which use one exact responsibility M-step
#'   per sweep and do not penalize the null mixture weight.
#' @param cal_obj If TRUE, record and monitor the ELBO, including the initial
#'   value. If FALSE, monitor parameter changes without storing ELBO history.
#' @param quantile_trans logical if set as TRUE perform normal quantile transform
#' on wavelet coefficients
#' @param thresh_lowcount numeric, used to check the wavelet coefficients have
#'  problematic distribution (very low dispersion even after standardization).
#'  Basically check if the median of the absolute value of the distribution of
#'   a wavelet coefficient is below this threshold. If yes, the algorithm discard
#'   this wavelet column from fitting, prior learning and noise estimation;
#'   its effect and posterior variance are set to zero.
#'   Set to 0 by default. It can be useful when analyzing sparse data from sequence
#'    based assay or small samples.
#' @param gridmult numeric used to control the number of components used in the mixture prior (see ashr package
#'  for more details). From the ash function:  multiplier by which the default grid values for mixsd differ from one another.
#'   (Smaller values produce finer grids.). Increasing this value may reduce computational time
#' @return An \code{EBmvFR} object. \code{fitted_func} contains effect curves
#'   in the original predictor units; \code{ind_fitted_func} contains the
#'   fitted centered response. \code{fitted_wc[[1]]} contains posterior means
#'   and \code{fitted_wc2[[1]]} posterior variances in the standardized design.
#'   \code{MLE_wc2[[1]]} retains standard errors, despite its historical name.
#'   \code{mixture_counts} stores summed responsibilities by predictor and
#'   prior group. \code{KL} contains augmented-posterior KL divergences to the
#'   current prior, \code{sigma2} the fitted noise variance, and
#'   \code{active_wc} the retained wavelet columns. \code{niter} counts
#'   completed sweeps and \code{converged} indicates whether the requested
#'   stopping criterion was reached. \code{ELBO} is populated only when
#'   \code{cal_obj=TRUE}.
#' @examples
#'
#'
#'library(ashr)
#'library(wavethresh)
#'set.seed(1)
#'#Example using curves simulated under the Mixture normal per scale prior
#'rsnr <- 1 #expected root signal noise ratio
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
#'#Running Empirical Bayes multivariate function regression
#'
#'out <- EBmvFR(Y,X )
#'#the easiest way to visualize the result is to use the plot_susiF function
#'
#'plot( f1, type="l", main="Estimated effect 1", xlab="")
#'lines(unlist(out$fitted_func[1,]),col='blue' )
#'abline(a=0,b=0)
#'legend(x= 60,
#'       y=3,
#'       lty= rep(1,2),
#'       legend = c("effect 1"," EBmvFR est "),
#'       col=c("black","blue" )
#')
#'plot( f2, type="l", main="Estimated effect 2", xlab="")
#'lines(unlist(out$fitted_func[5,]),col='green' )
#'abline(a=0,b=0)
#'legend(x= 85,
#'       y= 2,
#'       lty= rep(1,2),
#'       legend = c("effect 2"," EBmvFR est "),
#'       col=c("black","green" )
#')
#'
#'par(mfrow=c(1,1))
#' @importFrom stats var
#'
#' @export
#'
EBmvFR <- function(Y, X,
                   adjust=FALSE,
                   pos = NULL,
                   prior = "mixture_normal_per_scale",
                   verbose = TRUE,
                   maxit = 100,
                   tol = 0.1,
                   init_pi0_w= 1,
                   nullweight ,
                   control_mixsqp =  list(verbose=FALSE,
                                          eps = 1e-6,
                                          numiter.em = 4
                   ),
                   thresh_lowcount=0,
                   cal_obj=FALSE,
                   quantile_trans=FALSE,
                   gridmult= sqrt(2),
                   max_step_EM=1,

                   max_SNP_EM=100
)
{


  ####Cleaning input -----
  pt <- proc.time()
  if( prior %!in% c("mixture_normal", "mixture_normal_per_scale"))
  {
    stop("Error: provide valid prior input")
  }
  if(missing(nullweight))
  {
    nullweight <- 10#/(sqrt(nrow(X)))
  }
  if(missing(tol) && !cal_obj){
    tol <-10^-3
  }

  ## Input error messages
  if (!is.matrix(Y) || !is.numeric(Y) || !is.matrix(X) || !is.numeric(X) ||
      nrow(Y) != nrow(X) || nrow(Y) < 2L || ncol(X) < 1L || ncol(Y) < 4L ||
      any(!is.finite(Y)) || any(!is.finite(X))) {
    stop("Y and X must be finite numeric matrices with matching rows, at least two rows, four response columns and one predictor")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol < 0 ||
      length(maxit) != 1L || !is.finite(maxit) || maxit < 1 || maxit != floor(maxit)) {
    stop("tol must be non-negative and maxit must be a positive integer")
  }

  if (is.null(pos))
  {
    pos <- 1:dim(Y)[2]
  }

  #reshaping of the data
  if ( !(length(pos)==dim(Y)[2])) #miss matching positions and number of observations
  {
    stop("Error: number of position provided different from the number of column of Y")
  }



  if(!is.wholenumber(log2(dim(Y)[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check whether dim(Y) not equal to 2^J or if the data are unevenly spaced
  {

    inter_pol.obj <-interpol_mat(Y, pos)
    Y             <- inter_pol.obj$Y
    bp            <- inter_pol.obj$bp

    start_pos <- min( pos)
    end_pos <-max(pos)
    outing_grid   <- start_pos + (end_pos-start_pos)/(length(pos))*inter_pol.obj$grid
    if(verbose)
    {
      message( "Response matrix dimensions not equal to nx 2^J \n or unevenly spaced data \n interpolation procedure used")
    }
  }  else{

    outing_grid   <- pos
  }
  if(adjust){
    Y_org <- Y
    X_org <- X
  }
  # centering and scaling covariate
  X <- colScale(X)
  if (any(colSums(X^2) == 0)) stop("X must not contain constant columns")
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
  sigma2_init <- mean(update_Y[, setdiff(seq_len(ncol(update_Y)), lowc_wc),
                              drop = FALSE]^2)
  temp        <- init_prior(Y              = update_Y,
                            X              = X,
                            prior          = prior ,
                            v1             = v1,
                            indx_lst       = indx_lst,
                            lowc_wc        = lowc_wc,
                            control_mixsqp = control_mixsqp,
                            nullweight     = nullweight,
                            gridmult       = gridmult,
                            sigma2         = sigma2_init,
                            max_SNP_EM     = max_SNP_EM,
                            max_step_EM    = max_step_EM)
  G_prior     <- temp$G_prior


  #Recycled for the first step of the while loop
  obj   <-  init_EBmvFR_obj(G_prior=G_prior,
                                   Y=update_Y,
                                   X=X,
                                   lowc_wc=lowc_wc)

  obj   <- EBmvFR.workhorse(       obj            = obj,
                                   W              = W,
                                   X              = X,
                                   tol            = tol,
                                   lowc_wc        = lowc_wc,
                                   init_pi0_w     = init_pi0_w ,
                                   control_mixsqp = control_mixsqp ,
                                   indx_lst       = indx_lst,
                                   nullweight     = nullweight,
                                   cal_obj        = cal_obj,
                                   verbose        = verbose,
                                   maxit          = maxit,
                                   max_step_EM    = max_step_EM)
  #preparing output
  obj <- out_prep(       obj         = obj,
                         Y           = Y,
                         X           = X,
                         indx_lst    = indx_lst,
                         outing_grid = outing_grid
  )
  obj$runtime <- proc.time()-pt

  if(adjust){
    obj$Y_adjusted <-  Y_org-X_org%*%obj$fitted_func
  }
  return(obj)
}
