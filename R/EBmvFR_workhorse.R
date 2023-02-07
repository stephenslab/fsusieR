EBmvFR.workhorse <- function(EBmvFR.obj,
                             W,
                             X,
                             tol,
                             low_wc,
                             init_pi0_w ,
                             control_mixsqp ,
                             indx_lst,
                             lowc_wc,
                             nullweight,
                             cal_obj ,
                             verbose,
                             maxit){

  Y_f      <-  cbind( W$D,W$C)
  # numerical value to check breaking condition of while
  check <- 3*tol

  if(verbose){
    print("Initialization done")
  }
  ##### Start While -----
  iter <- 1


  while( (check >tol & iter <maxit))
  {
    for( j in 1:ncol(X))
    {

      if(verbose){
        print(paste("Fitting effect ", j,", iter" ,  iter ))
      }
      EBmvFR.obj   <-  fit_effect.EBmvFR (EBmvFR.obj = EBmvFR.obj,
                                          j          = j,
                                          X          = X,
                                          D          = W$D,
                                          C          = W$C,
                                          indx_lst   = indx_lst,
                                          lowc_wc    = lowc_wc)

    }#end for l in 1:L  -----



    sigma2    <- estimate_residual_variance(EBmvFR.obj,Y=Y_f,X)
    print(sigma2)
    EBmvFR.obj <- update_residual_variance(EBmvFR.obj, sigma2 = sigma2 )

    EBmvFR.obj <- update_prior( EBmvFR.obj,
                                max_step       = 100,
                                espsilon       = 0.0001,
                                init_pi0_w     = init_pi0_w ,
                                control_mixsqp = control_mixsqp,
                                indx_lst       = indx_lst,
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


  return(EBmvFR.obj)
}
