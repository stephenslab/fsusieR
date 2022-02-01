#'@title Sum of Single Functions
#'
#'
#'@description Implementation of the SuSiF method
#'@details tbd
#'
#'@param Y functional phenotype, matrix of size N by size J. The underlying algorithm uses wavelet which assume that J is of the form J^2. If J not a power of 2, susif internally remaps the data into grid of length 2^J
#'@param X matrix of size n by p in
#'@param L the number of effect to fit (if not specified set to =2)
#'@param pos vector of length J, corresponding to position/time pf the observed column in Y, if missing suppose that the observation are evenly spaced
#'@param prior specify the prior used in susif. Three choice are available "normal", "mixture_normal", "mixture_normal_per_scale"
#'@param verbose If \code{verbose = TRUE}, the algorithm's progress, and a summary of the optimization settings, are printed to the console.
#'@param plot_out If \code{plot_out = TRUE}, the algorithm's progress, and a summary of the optimization settings, are ploted.
#'@param tol A small, non-negative number specifying the convergence tolerance for the IBSS fitting procedure. The fitting procedure will halt when the difference in the variational lower bound, or \dQuote{ELBO} (the objective function to be maximized), is less than \code{tol}. Currently checking the PIP
#'@param maxit Maximum number of IBSS iterations to perform.
#'@export
susif <- function( Y,X, L = 2,
                   pos = NULL,
                   prior = "mixture_normal_per_scale",
                   verbose = TRUE,
                   plot_out = TRUE,
                   maxit = 100,
                   tol = 1e-6

)
{

  if( prior %!in% c("normal", "mixture_normal", "mixture_normal_per_scale"))
  {
    stop("Error: provide valid prior input")
  }
  prior  = match.arg(prior)

  ## Input error messages

  if (is.null(pos))
  {
    pos <- 1:dim(Y)[2]
  }


  #reshaping of the data
  if ( !(length(pos)==dim(Y)[2])) #miss matching positions and number of observations
  {
    stop("Error: number of position provided different from number of column of Y")
  }
  orignal_Y <-Y

  i

  if(!is.wholenumber(log2(dim(Y)[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check wether dim(Y) not equal to 2^J or if the data are unevenly spaced
  {

    inter_pol.obj <- interpol_mat(Y, pos)
    Y             <- inter_pol.obj$Y
    bp            <-  inter_pol.obj$bp
    outing_grid   <- inter_pol.obj$grid
    if(verbose ==TRUE)
    {
      print( "Response matrix dimensions not equal to nx 2^J \n or unevenly spaced data \n interpolation procedure used")
    }
  }  else{

    outing_grid <- 1:dim(Y)[2]
  }


  ### Definition of some static parameters ---
  indx_lst <-  gen_wavelet_indx(log2(dim(Y_f)[2]))
  v1 <- rep(1, dim(X)[1])### used in fit_lm to add a column of 1 in the design matrix


  ### Definition of some dynamic parameters ---

  G_prior     <-  init_prior(Y_f,X,prior,v1 , indx_lst  )
  susiF.obj   <-  init_susiF_obj(L=L, G_prior, Y,X)


  #Wavelet transform of the inputs

  W <- DWT2(Y)
  update_D <- W
  Y_f <- cbind( W$D,W$C) #Using a column like phenotype
  update_Y <- Y_f # temporary matrix that will be regularly updated




  #IBSS   ---------
  while(check >tol & (h/L) <maxit)
  {
   for( l in 1:L)
   {
     tt <- cal_Bhat_Shat(update_Y,X,v1)
     Bhat <- tt$Bhat
     Shat <- tt$Shat #UPDATE. could be nicer


     tpi <-  get_pi(susiF.obj,l)
     G_prior <- update_prior(G_prior,
                             tpi= tpi )
     EM_out  <- EM_pi(G_prior= G_prior,
                      Bhat =  Bhat,
                      Shat = Shat,
                      indx_lst =  indx_lst
                      )
     susiF.obj <-  update_susiF_obj(susiF.obj = susiF.obj ,
                                    l = l,
                                    EM_pi = EM_out
                                    )



     G_prior <- update_prior(G_prior,
                             tpi= EM_out$tpi_k )

     fitted_wc_col[[l]]   <- post_mat_mean( G_prior , Bhat, Shat )
     fitted_wc_col2[[l]]  <- post_mat_sd  (G_prior , Bhat, Shat )
     lBF <- EM_out$lBF
     alpha_col[[l]] <-cal_zeta(lBF)
   }

  }
}
