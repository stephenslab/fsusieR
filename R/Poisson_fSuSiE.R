



#' @export


Pois_fSuSiE <- function(Y,
                        Z,
                        X,
                        L=3,
                        scaling= NULL,
                        L_start=3,
                        reflect =FALSE,
                        verbose=TRUE,
                        init_b_pm,
                        tol= 1e-3,
                        tol_vga_pois=1e-5,
                        control_mixsqp=  list(verbose=FALSE,
                                              eps = 1e-6,
                                              numiter.em = 4
                        ),
                        thresh_lowcount=1e-2,
                        prior_mv=  "mixture_normal_per_scale",
                        post_processing=c( "HMM","smash","TI","none"),
                        gridmult=sqrt(2),
                        nullweight.mrash=10,
                        init_pi0_w.mrash=10,
                        cov_lev=0.95,
                        min_purity     =0.5,
                        greedy=TRUE,
                        backfit=TRUE,
                        tol.mrash=1e-3,
                        verbose.mrash=TRUE,
                        maxit.mrash=10,
                        cal_obj.mrash=FALSE,
                        maxit.fsusie=50,
                        cal_obj.fsusie=FALSE,
                        max_SNP_EM     = 100,
                        max_step_EM    = 1,
                        cor_small=FALSE,
                        max.iter=3,
                        init_pi0_w=1,
                        nullweight_fsusie= .001,
                        print=FALSE
)
{
  ####Changer les calcul d'objective -----
  if(missing(X)&missing(Z)){
    stop("Please provide a Z or a X matrix")
  }
  
  fit_approach <- "both"
  if(missing(X)){
    print("No correlated covariate provided, the algorithm will perform penalized regression only")
    fit_approach <- "penalized"
  }
  if(missing(Z)){
    print("No Z matrix provided,  the algorithm will perform fine-mapping only")
    fit_approach <- "fine_mapping"
    
  }
  
  ##initiatilzation -----
  init=TRUE
  J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
  if(reflect){
    tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx #### indx of interest at the end
  }else{
    idx_out <- 1: ncol(Y)
  }
  #### to avoid 0 in Y_min to correct at the end
  Y <- Y
  
  
  if(fit_approach %in% c('both',"fine_mapping")){
    tidx <- which(apply(X,2,var)==0)
    if( length(tidx)>0){
      warning(paste("Some of the columns of X are constants, we removed" ,length(tidx), "columns"))
      X <- X[,-tidx]
    }
    X <- fsusieR:::colScale(X)
    names_colX <-  colnames(X)
  }
  
  if(fit_approach %in% c('both',"penalized")){
    tidx <- which(apply(Z,2,var)==0)
    if( length(tidx)>0){
      warning(paste("Some of the columns of Z are constants, we removed" ,length(tidx), "columns"))
      Z <- Z[,-tidx]
    }
  }
  
  if( is.null( scaling)){
    scaling = rep(1, nrow(Y))
  }else{
    if( ! (length(scaling)== nrow(Y))){
      warning(paste("scaling shoudl have a length equal to number of row  in Y"  ))
    }
  }
  
  indx_lst <-  fsusieR::gen_wavelet_indx(log2(ncol(Y)))
  
  
  
  
  init_val_pois<- c(log(Y+1))
  beta_pois <- 0* c(log(Y+1))
  sigma2_pois=1
  
  ##initiatilzation for count data -----
  Mu_pm<- Y
  iter=1
  beta_pois <- 0* c(log(Mu_pm +1))
  check <- 3*tol
  
  b_pm <- 0* Mu_pm
  fm_pm <- 0* Mu_pm
  
  
  while( check >tol & iter <=  max.iter ){
    
    
    if ( iter ==1 ){
      tt= ebpm_normal(c(Y),s= rep( scaling, ncol(Y)) )
      Mu_pm <- matrix( tt$posterior$mean_log,byrow = FALSE, ncol=ncol(Y))
      
    }else{
      
      
      tt <-    pois_mean_GP(x=c(Y),
                            prior_mean = c(Mu_pm_init),
                            s =  rep( scaling, ncol(Y)),
                            prior_var = sigma2_pois )
      Mu_pm <- matrix( tt$posterior$posteriorMean_latent,byrow = FALSE, ncol=ncol(Y))
      Mu_pv <- matrix( tt$posterior$posteriorVar_latent ,byrow = FALSE, ncol=ncol(Y))
    }
    
    
    
    if(verbose){
      print( paste('Posterior log intensity computed for iter ',iter))
    }
    
    
    
    
    
    
    if(init){
      
      tmp_Mu_pm <- fsusieR::colScale(Mu_pm, scale = FALSE)#potentially run smash on colmean
      lowc_wc <-  which_lowcount(tmp_Mu_pm,
                                 thresh_lowcount=thresh_lowcount)
      
      if( length(lowc_wc) > (ncol(tmp_Mu_pm)-10)){
        out <-NULL
        
        return(out)
      }
      
      
      W <- list( D = tmp_Mu_pm [, -ncol(tmp_Mu_pm )],
                 C = tmp_Mu_pm [,  ncol(tmp_Mu_pm )])
      if (fit_approach %in% c("both", "penalized")){
        temp <- fsusieR:: init_prior(Y              = tmp_Mu_pm,
                                     X              = Z ,
                                     prior          = prior_mv ,
                                     v1             = v1,
                                     indx_lst       = indx_lst,
                                     lowc_wc        = lowc_wc,
                                     control_mixsqp = control_mixsqp,
                                     nullweight     = nullweight.mrash,
                                     gridmult       = gridmult )
        G_prior     <- temp$G_prior
        
        
        #Recycled for the first step of the while loop
        EBmvFR.obj   <-  fsusieR::init_EBmvFR_obj(G_prior = G_prior,
                                                  Y       = Y,
                                                  X       = Z
        )
        print('Done initializing EBmvFR.obj')
      }
      if(fit_approach %in%c("both","fine_mapping")){
        
        temp <- fsusieR:: init_prior(    Y              = tmp_Mu_pm,
                                         X              = X ,
                                         prior          = prior_mv ,
                                         
                                         indx_lst       = indx_lst,
                                         lowc_wc        = lowc_wc,
                                         control_mixsqp = control_mixsqp,
                                         nullweight     = nullweight.mrash,
                                         gridmult       = gridmult )
        G_prior     <- temp$G_prior
        
        
        #Recycled for the first step of the while loop
        susiF.obj   <-  fsusieR::init_susiF_obj(L_max   = L,
                                                G_prior = G_prior,
                                                Y       = tmp_Mu_pm,
                                                X       = X,
                                                L_start = L_start,
                                                greedy  = greedy,
                                                backfit = backfit
        )
        print('Done initializing susiF.obj')
        
      }
      tmp_Mu_pm_pen <- 0*tmp_Mu_pm
      tmp_Mu_pm_fm  <- 0*tmp_Mu_pm
      init=FALSE
    }
    
    #### fit EBmvFR ----
    if(fit_approach%in% c("both", "penalized")){
      tmp_Mu_pm_pen <- Mu_pm  -  fm_pm#potentially run smash on colmean
      
      t_mean_EBmvFR <-  apply(tmp_Mu_pm_pen,2, mean )
      tmp_Mu_pm_pen <- fsusieR::colScale(tmp_Mu_pm_pen, scale=FALSE)
      W <- list( D = tmp_Mu_pm [, -ncol(tmp_Mu_pm_pen )],
                 C = tmp_Mu_pm [,  ncol(tmp_Mu_pm_pen )])
      
      
      ### TODO: Maybe use better restarting point for EBmvFR.obj
      EBmvFR.obj   <- fsusieR::EBmvFR.workhorse( obj     = EBmvFR.obj,
                                                 W              = W,
                                                 X              = Z,
                                                 tol            = tol.mrash,
                                                 lowc_wc        = lowc_wc  ,
                                                 init_pi0_w     = init_pi0_w.mrash ,
                                                 control_mixsqp = control_mixsqp ,
                                                 indx_lst       = indx_lst,
                                                 nullweight     = nullweight.mrash,
                                                 cal_obj        = cal_obj.mrash,
                                                 verbose        = FALSE,
                                                 maxit          = maxit.mrash,
                                                 max_step_EM =1
      )
      if(verbose){
        print( paste('Posterior of EB regression coefficient computed for iter ',iter))
      }
      b_pm <-   Z%*%  EBmvFR.obj$fitted_wc[[1]]
      
      if( fit_approach== "penalized")
        mat_mean <-   matrix( t_mean_EBmvFR , byrow = TRUE,
                              nrow=nrow(X), ncol=ncol(Y))
      
    }else{
      b_pm <- 0* tmp_Mu_pm_pen
      
    }
    
    if(fit_approach%in% c("both", "fine_mapping")){
      tmp_Mu_pm_fm <- Mu_pm -  b_pm#potentially run smash on colmean
      tmp_Mu_pm_fm <- fsusieR::colScale(tmp_Mu_pm_fm, scale=FALSE)
      susiF.obj     <- susiF (
        Y              =  tmp_Mu_pm_fm ,
        X               = X ,
        L               = L,
        tol             = tol,
        control_mixsqp  = control_mixsqp ,
        nullweight      = nullweight.mrash,
        cal_obj         = cal_obj.fsusie,
        verbose         = verbose,
        cov_lev         = cov_lev,
        min_purity      = min_purity,
        
        cor_small       = cor_small,
        maxit           = maxit.fsusie,
        post_processing = post_processing)
      
      
      
      
      
      
      fm_pm <- X%*%Reduce("+",lapply(1:length(susiF.obj$cs),
                                     function(l)
                                       t(susiF.obj$fitted_func[[l]]%*% t(susiF.obj$alpha[[l]]))
      )
      )
      mat_mean <-   matrix(Mu_pm -fm_pm , byrow = TRUE,
                           nrow=nrow(X), ncol=ncol(Y))
    }else{
      fm_pm <-0* tmp_Mu_pm_fm
      susiF.obj   <- NULL
    }
    
    
    resid <- Mu_pm -mat_mean -fm_pm-b_pm
    #not correct to work on later
    sigma2_pois <- var(c(resid ))
    #print(sigma2_pois)
    Mu_pm <- mat_mean +fm_pm+b_pm#update
    Mu_pm_init <-Mu_pm
    print(    susiF.obj$cs)
    iter=iter+1
    ##include mr.ash
    
    if (print){
      par (mfrow=c(1,2))
      
      plot ( Y[1,], col="blue")
      points ( exp(Mu_pm  [1,]))
      lines(exp(Mu_pm  [1,]), col="green")
      
      
      plot( Y[1,],exp(Mu_pm  [1,]))
      
      abline(a=0,b=1)
      par (mfrow=c(1,1))
    }
    
  }
  
  
  
  
  
  tt_all <-exp(Mu_pm   )
  
  
  
  if( fit_approach ==   "both" )
  {
    susiF.obj <- fsusieR::update_cal_pip(susiF.obj)
    out <- list( Mu_pm=Mu_pm,
                 susiF.obj=susiF.obj,
                 EBmvFR.obj=EBmvFR.obj,
                 fitted = tt_all[,idx_out] )
  }
  
  if( fit_approach ==   "fine_mapping" )
  {
    susiF.obj <- fsusieR::update_cal_pip(susiF.obj)
    out <- list( Mu_pm=Mu_pm,
                 susiF.obj=susiF.obj,
                 fitted = tt_all[,idx_out]  )
  }
  if( fit_approach ==   "penalized")
  {
    out <- list( Mu_pm=Mu_pm,
                 EBmvFR.obj=EBmvFR.obj,
                 fitted = tt_all[,idx_out]  )
  }
  return(out)
  
}
