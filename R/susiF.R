#' @title Sum of Single Functions
#'
#' @description Implementation of the SuSiF method
#'
#' @details tbd
#'
#' @param Y functional phenotype, matrix of size N by size J. The
#'   underlying algorithm uses wavelet which assume that J is of the
#'   form J^2. If J not a power of 2, susif internally remaps the data
#'   into grid of length 2^J
#'
#' @param X matrix of size n by p contains the covariates
#'
#' @param L the number of effect to fit (if not specified set to =2)
#'
#' @param pos vector of length J, corresponding to position/time pf
#' the observed column in Y, if missing suppose that the observation
#' are evenly spaced
#'
#' @param prior specify the prior used in susif. Three choice are
#' available "normal", "mixture_normal", "mixture_normal_per_scale"
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings, are printed to the
#' console.
#'
#' @param plot_out If \code{plot_out = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings, are ploted.
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
                  plot_out = TRUE,
                  maxit = 100,
                  tol = 1e-3,
                  cov_lev = 0.95,
                  min.purity=0.5,
                  filter.cs =TRUE,
                  init_pi0_w = 1,
                  control_mixsqp =  list(
                                        eps = 1e-6,
                                        numiter.em = 40,
                                        verbose = FALSE
                                        )
)
{
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

  if( L==1)
  {
    tt   <- cal_Bhat_Shat(update_Y,X,v1)
    Bhat <- tt$Bhat
    Shat <- tt$Shat #UPDATE. could be nicer
    tpi  <- get_pi(susiF.obj,1)
    G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

    EM_out  <- EM_pi(G_prior       = G_prior,
                     Bhat          = Bhat,
                     Shat          = Shat,
                     indx_lst      = indx_lst,
                     init_pi0_w    = init_pi0_w,
                     control_mixsqp = control_mixsqp
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
      for( l in 1:L)
      {

        h <- h+1
        tt <- cal_Bhat_Shat(update_Y,X,v1)
        Bhat <- tt$Bhat
        Shat <- tt$Shat #UPDATE. could be nicer
        tpi <-  get_pi(susiF.obj,l)
        G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

        EM_out  <- EM_pi(G_prior       = G_prior,
                         Bhat          = Bhat,
                         Shat          = Shat,
                         indx_lst      = indx_lst,
                         init_pi0_w    = init_pi0_w,
                         control_mixsqp = control_mixsqp
        )

        susiF.obj <-  update_susiF_obj(susiF.obj   = susiF.obj ,
                                       l           = l,
                                       EM_pi       = EM_out,
                                       Bhat        = Bhat,
                                       Shat        = Shat,
                                       indx_lst    = indx_lst
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
      sigma2    <- estimate_residual_variance(susiF.obj,Y=Y_f,X)
      sigma2
      susiF.obj <- update_residual_variance(susiF.obj, sigma2 = sigma2 )
      susiF.obj <- update_KL(susiF.obj,
                             X,
                             D= W$D,
                             C= W$C , indx_lst)
       susiF.obj <- update_ELBO(susiF.obj,
                                 get_objective( susiF.obj = susiF.obj,
                                                Y         = Y_f,
                                                X         = X,
                                                D         = W$D,
                                                C         = W$C,
                                                indx_lst  = indx_lst
                                                )
                                          )



      if(length(susiF.obj$ELBO)>1 )#update parameter convergence,
      {
        check <- abs(diff(susiF.obj$ELBO)[(length( susiF.obj$ELBO )-1)])

      }
    }#end while
  }


  #preparing output
  susiF.obj <- out_prep(susiF.obj  = susiF.obj,
                        Y          = Y,
                        X          = X,
                        indx_lst   = indx_lst,
                        filter.cs  = filter.cs
                        )
  return(susiF.obj)
}
