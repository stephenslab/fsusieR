susiF.workhorse <- function(susiF.obj,
                            W,
                            X,
                            tol,
                            low_wc,
                            init_pi0_w ,
                            control_mixsqp ,
                            indx_lst,
                            lowc_wc,
                            nullweight,
                            cal_obj,
                            verbose,
                            cov_lev,
                            min.purity,
                            maxit,
                            tt){

  G_prior  <- get_G_prior(susiF.obj )
  Y_f      <-  cbind( W$D,W$C)
  update_Y <- Y_f
  # numerical value to check breaking condition of while
  # numerical value to check breaking condition of while
  check <- 3*tol
  init        <- TRUE
  v1       <-  rep(1, dim(X)[1])


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
          l         = (l-1)  ,
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
                                                Bhat     = tt$Bhat,
                                                Shat     = tt$Shat,
                                                lowc_wc  = lowc_wc,
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
      sigma2    <- estimate_residual_variance(susiF.obj,Y=Y_f,X)
      #print(sigma2)
      susiF.obj <- update_residual_variance(susiF.obj, sigma2 = sigma2 )
      susiF.obj <- test_stop_cond(susiF.obj = susiF.obj,
                                  check     = check,
                                  cal_obj   = cal_obj,
                                  Y         = Y_f,
                                  X         = X,
                                  D         = W$D,
                                  C         = W$C,
                                  indx_lst  = indx_lst)
      #  print(susiF.obj$alpha)
      #print(susiF.obj$ELBO)
      check <- susiF.obj$check



      iter <- iter +1


    }#end while
  }#end else in if(L==1)

 return(susiF.obj)
}
