#' @title Sum of Single Functions
#'
#' @description Implementation of the SuSiF method
#'
#' @details Implementation of the SuSiF method
#'
#'
#' @param Y functional phenotype, matrix of size N by size J. The
#'   underlying algorithm uses wavelet which assume that J is of the
#'   form J^2. If J not a power of 2, susif internally remaps the data
#'   into grid of length 2^J
#'
#' @param X matrix of size n by p contains the covariates
#'
#' @param L upper bound on the number of effect to fit (if not specified set to =2)
#'
#' @param pos vector of length J, corresponding to position/time pf
#' the observed column in Y, if missing suppose that the observation
#' are evenly spaced
#'
#' @param prior specify the prior used in susif. Three choice are
#' available "mixture_normal_per_scale", "mixture_normal". Default "mixture_normal_per_scale",
#' if this susiF is too slow consider using  "mixture_normal" (up to 40% faster), but this may results in
#' oversmoothing the estimated curves
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings, are printed to the
#' console.
#'
#'
#' @param tol A small, non-negative number specifying the convergence
#' tolerance for the IBSS fitting procedure. The fitting procedure
#' will halt when the difference in the variational lower bound, or
#' \dQuote{ELBO} (the objective function to be maximized), is less
#' than \code{tol}.
#'
#' @param maxit Maximum number of IBSS iterations to perform.
#'
#' @param cov_lev numeric between 0 and 1, corresponding to the
#' expected level of coverage of the cs if not specified set to 0.95
#'
#' @param min.purity minimum purity for estimated credible sets
#' @param filter.cs logical, if TRUE filter the credible set (removing low purity cs and cs with estimated prior equal to 0)
#' @param init_pi0_w starting value of weight on null compoenent in mixsqp (between 0 and 1)
#' @param control_mixsqp list of parameter for mixsqp function see\link{\code{mixsqp}}
#' @param  cal_obj logical if set as true compute ELBO for convergence monitoring
#' @param quantile_trans logical if set as true perform normal quantile transform on wavelet coefficients
#' @param L_start number of effect initialized at the start of the algorithm
#' @param nullweight numeric value for penalizing likelihood at point mass 0 (should be between 0 and 1)
#' (usefull in small sample size)
#' @param greedy logical, if true allow greedy search for extra effect (up to L specify by the user). Set as TRUE by default
#' @param backfit logical, if true allow discarding effect via backfitting. Set as true by default as TRUE. We advise to keep it as TRUE
#'
#' @examples
#'
#'library(susiF.alpha)
#'library(ashr)
#'set.seed(1)
#'#Example using curves simulated under the Mixture normal per scale prior
#'rsnr <- 0.2 #expected root signal noise ratio
#'N <- 100    #Number of individuals
#'P <- 10     #Number of covariates
#'pos1 <- 1   #Position of the causal covariate for effect 1
#'pos2 <- 5   #Position of the causal covariate for effect 1
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
#'#the easiest way to vizualize the result is to use the plot_susiF function
#'
#'plot_susiF(out)
#'
#'#You can also acces the information directly in the out of susiF. as follow
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
                  tol = 10^-3,
                  cov_lev = 0.95,
                  min.purity=0.5,
                  filter.cs =TRUE,
                  init_pi0_w= 1,
                  nullweight ,
                  control_mixsqp =  list(verbose=FALSE,eps = 1e-6,numiter.em = 4),
                  thresh_lowcount,
                  cal_obj=FALSE,
                  L_start = 3,
                  quantile_trans=FALSE,
                  greedy =TRUE,
                  backfit =TRUE
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



  if(!is.wholenumber(log2(dim(Y)[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check whether dim(Y) not equal to 2^J or if the data are unevenly spaced
  {

    inter_pol.obj <- interpol_mat(Y, pos)
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
    start_pos <- min( pos)
    end_pos <-max(pos)
    outing_grid   <- start_pos + (end_pos-start_pos)/(length(pos))*inter_pol.obj$grid
  }
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

  if(!missing(thresh_lowcount)){
     lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)
     if(verbose){
       print( paste("Discarding ", length(lowc_wc), "wavelet coefficients out of ", ncol(Y_f)))
     }
     if(length(lowc_wc)> (ncol(Y_f )-3)){
       stop("almost all the wavelet coefficients are null, consider using univariate fine mapping")
     }
  }else{
    lowc_wc <- NULL
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
                            nullweight     = nullweight )
  G_prior     <- temp$G_prior
  tt          <- temp$tt
  init        <- TRUE
  #Recycled for the first step of the while loop
  susiF.obj   <-  init_susiF_obj(L_max=L,
                                 G_prior=G_prior,
                                 Y=Y,
                                 X=X,
                                 L_start=L_start,
                                 greedy=greedy,
                                 backfit=backfit)

  # numerical value to check breaking condition of while
  check <- 1
  # l=1
  #init        <- FALSE
  #susiF.obj$fitted_wc[[1]] <-tt$Bhat
  #update_D  <-  W$D -    (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% (susiF.obj$fitted_wc[[l]][,-indx_lst[[length(indx_lst)]]] )
  #update_C  <-  W$C -     (X*rep(susiF.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% susiF.obj$fitted_wc[[l]][,indx_lst[[length(indx_lst)]]]
  #update_Y  <- cbind(  update_D, update_C)
  ####Start while -----
  if(verbose){
    print("Initialization done")
  }

  if( susiF.obj$L_max==1)
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
    ##### Start While -----
    iter <- 1
    while( (check >tol & iter <maxit))
    {
      for( l in 1:susiF.obj$L)
      {

        #print(susiF.obj$alpha[[l]])
        update_Y  <-  cal_partial_resid(
          susiF.obj = susiF.obj,
          l         =  (l-1)  ,
          X         = X,
          D         = W$D,
          C         = W$C,
          indx_lst  = indx_lst
        )

        if(verbose){
          print(paste("Fitting effect ", l,", iter" ,  iter ))
        }

        if(init){#recycle operation used to fit the prior

            EM_out <- gen_EM_out (tpi_k= get_pi_G_prior(G_prior),
                                   lBF  = log_BF  (G_prior,
                                                   Bhat=tt$Bhat,
                                                   Shat=tt$Shat,
                                                   lowc_wc=lowc_wc,
                                                   indx_lst = indx_lst
                                                   )
                                  )
            init <- FALSE
        }else{
          tt <- cal_Bhat_Shat(update_Y,X,v1, lowc_wc =lowc_wc )

          tpi <-  get_pi(susiF.obj,l)
          G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

          EM_out  <- EM_pi(G_prior        = G_prior,
                           Bhat           = tt$Bhat,
                           Shat           = tt$Shat,
                           indx_lst       = indx_lst,
                           init_pi0_w     = init_pi0_w,
                           control_mixsqp = control_mixsqp,
                           lowc_wc        = lowc_wc,
                           nullweight     = nullweight
          )
        }


        #print(h)
        # print(EM_out$lBF[1:10])
        susiF.obj <-  update_susiF_obj(susiF.obj   = susiF.obj ,
                                       l           = l,
                                       EM_pi       = EM_out,
                                       Bhat        = tt$Bhat,
                                       Shat        = tt$Shat,
                                       indx_lst    = indx_lst,
                                       lowc_wc     = lowc_wc
        )

      }#end for l in 1:L  -----

      ####Check greedy/backfit and stopping condition -----
      susiF.obj <- greedy_backfit (susiF.obj,
                                  verbose    = verbose,
                                  cov_lev    = cov_lev,
                                  X          = X,
                                  min.purity = min.purity
                                  )

     susiF.obj <- test_stop_cond(susiF.obj = susiF.obj,
                                  check    = check,
                                 cal_obj   = cal_obj,
                                 Y         = Y_f,
                                 X         = X,
                                 D         = W$D,
                                 C         = W$C,
                                indx_lst  = indx_lst)
   #  print(susiF.obj$alpha)
     #print(susiF.obj$ELBO)
    check <- susiF.obj$check

    sigma2    <- estimate_residual_variance(susiF.obj,Y=Y_f,X)
    susiF.obj <- update_residual_variance(susiF.obj, sigma2 = sigma2 )

    iter <- iter +1


    }#end while
  }#end else in if(L==1)

  #preparing output
  susiF.obj <- out_prep(susiF.obj  = susiF.obj,
                        Y          = Y,
                        X          = X,
                        indx_lst   = indx_lst,
                        filter.cs  = filter.cs,
                        outing_grid=outing_grid
                        )
  susiF.obj$runtime <- proc.time()-pt
  return(susiF.obj)
}
