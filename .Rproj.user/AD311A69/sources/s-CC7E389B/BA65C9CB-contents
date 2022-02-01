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
#'@param cov_lev numeric between 0 and 1, corresponding to the expected level of coverage of the cs if not specified set to 0.95
#'@export
susif <- function( Y,X, L = 2,
                   pos = NULL,
                   prior = "mixture_normal_per_scale",
                   verbose = TRUE,
                   plot_out = TRUE,
                   maxit = 100,
                   tol = 1e-6,
                   cov_lev = 0.95

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
  update_Y <- cbind( W$D,W$C) #Using a column like phenotype, temporary matrix that will be regularly updated




  #IBSS   ---------
  while(check >tol & (h/L) <maxit)
  {
   for( l in 1:L)
   {
     tt <- cal_Bhat_Shat(update_Y,X,v1)
     Bhat <- tt$Bhat
     Shat <- tt$Shat #UPDATE. could be nicer
     tpi <-  get_pi(susiF.obj,l)
     G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

     EM_out  <- EM_pi(G_prior  = G_prior,
                      Bhat     =  Bhat,
                      Shat     = Shat,
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
                                      susiF.obj = susiF_obj,
                                      l         = 1,
                                      X         = X,
                                      D         = W$D,
                                      C         = W$C,
                                      L         = 2,
                                      indx_lst  = indx_lst
                                    )


     if(L==1){
       break
     }
     if(ceiling(h/L)>2)#update parameter convergence, no ELBO for the moment
     {
       check <-0
       for( tt in 0:(L-1))
       {
         check <-  check + var( susiF.obj$alpha_hist[[h-tt]] -susiF.obj$alpha_histalpha_col_hist[[h-L-tt]])
       }
     }
   }#end for l in 1:L
  }#end while

  #preparing the output



  ##### Tidding output ----
  temp <- wd(rep(0, dim(Y_f)[2]))
  for ( l in 1:L)
  {


    temp$D        <-    (alpha_col[[l]])%*%fitted_wc_col[[l]][,-indx_lst[[length(indx_lst)]]]
    temp$C[length(temp$C)] <- (alpha_col[[l]])%*%fitted_wc_col[[l]][,indx_lst[[length(indx_lst)]]]
    fitted_func[[l]] <- wr(temp)


  }


  ########################
  #Individual prediction -----
  ########################


  #preparing output
  susiF.obj <- update_cal_pip     (susiF.obj)
  susiF.obj <- update_cal_cs      (susiF.obj, cov_lev)
  susiF.obj <- update_cal_indf    (susif.obj,Y,X)
  susiF.obj <- update_cal_fit_func(susif.obj)


  return(susiF.obj)
}
