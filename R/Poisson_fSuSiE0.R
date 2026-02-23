#' @export
Pois_fSuSiE0 <- function(Y,
                        X = NULL,
                        L = 3,
                        scaling = NULL,

                        reflect = FALSE,
                        verbose = TRUE,
                        tol = 1e-3,
                        maxit_outer = 2,
                        maxit_inner = 10 ,
                        control_mixsqp = list(verbose = FALSE,
                                              eps = 1e-6,
                                              numiter.em = 4),
                        nullweight = .10,
                        cov_lev = 0.95,
                        min_purity = 0.5,
                        cor_small = FALSE,
                        post_processing = "smash",
                        filter.number = 1,
                        family = 'DaubExPhase',
                        print=TRUE,
                        update_Mu_each_iter = TRUE,
                        True_intensity=NULL,
                        s2=1) {

  # Validate inputs
  if (is.null(X) ) {
    stop("Please provide X matrix")
  }
  if (is.null(Y) ) {
    stop("Please provide Y matrix")
  }
  obj_fn <- pois_mean_GP_opt_obj
  grad_fn <- pois_mean_GP_opt_obj_gradient
  ebnm_params = list(mode = 0)
  # Determine fit approach
  has_X <- !is.null(X)


  # Match arguments

  # Handle reflection for non-dyadic length
  J <- log2(ncol(Y))
  if ((J %% 1) != 0) reflect <- TRUE

  if (reflect) {
    tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i, ]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx
  } else {
    idx_out <- 1:ncol(Y)
  }

  # Remove constant columns

    tidx <- which(apply(X, 2, var) == 0)
    if (length(tidx) > 0) {
      warning(paste("Removed", length(tidx), "constant columns from X"))
      X <- X[, -tidx, drop = FALSE]
    }

    est_effect_fm <- matrix(0, nrow =ncol(X), ncol = ncol(Y))      # X effects




  # Setup scaling
  if (is.null(scaling)) {
    scaling <- rep(1, nrow(Y))
  } else {
    if (length(scaling) != nrow(Y)) {
      stop("scaling must have length equal to nrow(Y)")
    }
  }

  # ============================================================================
  # STEP 1 (INITIALIZATION): Compute initial Mu_pm via empirical Bayes
  # ============================================================================




    susiF.obj=susiF(log1p(Y),X=X)
    if (length(susiF.obj$cs) > 0) {
      est_effect_fm=Reduce("+", lapply(1:length(susiF.obj$cs), function(l) {
        t(susiF.obj$fitted_func[[l]] %*% t(susiF.obj$alpha[[l]]))
      }))


      Eb_pm <- X %*% est_effect_fm
      B_pm <- X %*% est_effect_fm
    } else {
      tt <- ebpm_normal(c(Y), s = rep(scaling, ncol(Y)))
      Eb_pm <- matrix(tt$posterior$mean_log, byrow = FALSE, ncol = ncol(Y))
    }






  # Initialize components
  alpha_0 <- rowMeans(Eb_pm)  # Intercept per sample
  Theta_pm <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))  # Z effects

  sigma2 <-  .1  # Initial variance

  # ============================================================================
  # COORDINATE ASCENT ALGORITHM (Algorithm 5 from supplement)
  # ============================================================================
  converged <- FALSE
  iter <- 1
  elbo_hist <- numeric(maxit_outer)
  W = (t(wavethresh::GenW(ncol(Y),filter.number,family)))[-1,]
  # Store old values for convergence
  alpha_0_old <- alpha_0
  Theta_pm_old <- Theta_pm
  B_pm_old <- est_effect_fm
  Mu_pm <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  Mu_pv <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))
  Eb_pm <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  Eb_pv <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))
  KL    <- vector("list",  nrow(Y) )
  while (!converged && iter <= maxit_outer) {

    if (verbose) cat("\n=== Outer iteration", iter, "===\n")


    for (i in 1:nrow(Y)) {

      n = ncol(Y)
      mu_pm = rep(0,n)
      mu_pv = rep(1/n,n)
      # get m, s^2
      opt = optim(c(mu_pm, log(mu_pv)),
                  fn = obj_fn,
                  gr = grad_fn,
                  lower = c(rep(-30, n), rep(-20, n)),
                  upper = c(rep(30, n), rep(10, n)),
                  x =Y[i, ],
                  s = scaling[i],
                  beta =B_pm[i, ],
                  sigma2 = sigma2, n = n,
                  method ='L-BFGS-B')

      Mu_pm[i,] = opt$par[1:n]
      Mu_pv[i,] = exp(opt$par[(n+1):(2*n)])
      qb = smash_dwt(x=Mu_pm[i,],
                     sigma=sqrt(sigma2),
                     filter.number=filter.number,
                     family=family,
                     ebnm_params=ebnm_params,W=W)


      Eb_pm[i, ] = qb$mu.est
      Eb_pv[i, ] = qb$mu.est.var + Eb_pm[i, ]^2
      KL[i]=qb$dKL
    }




    susiF.obj <- susiF(
      Y = Mu_pm,
      X = X,
      L = L,
      tol = tol,
      control_mixsqp = control_mixsqp,
      nullweight = nullweight,
      cal_obj = FALSE,
      verbose = FALSE,
      cov_lev = cov_lev,
      min_purity = min_purity,
      maxit = maxit_inner,
      cor_small = cor_small,
      post_processing = post_processing
    )
    if (length(susiF.obj$cs) > 0) {
      est_effect_fm=Reduce("+", lapply(1:length(susiF.obj$cs), function(l) {
        t(susiF.obj$fitted_func[[l]] %*% t(susiF.obj$alpha[[l]]))
      }))


      B_pm <- X %*% est_effect_fm
    }

    sigma2 = mean(Mu_pm^2+Mu_pv+Eb_pv-2*Mu_pm*Eb_pm)
    sigma2 = max(sigma2, 1e-6)

    Eb_pm=B_pm

    obj_temp=0
    for ( i in 1:nrow(Y)){
      obj_temp=obj_temp+pois_smooth_split_obj(x=Y[i,],s=scaling[i],
                                              m =Mu_pm[i,],
                                              s2=Mu_pv[i,],
                                              Eb=Eb_pm[i,],
                                              Eb2=Eb_pv[i,],
                                              sigma2=sigma2,
                                              KLb=  KL[i][[1]])
    }
    elbo_hist [iter+1] = obj_temp



    #

    if (print){

      if (is.null(True_intensity)){
        par (mfrow=c(3,1))
        plot( log1p(Y ),  (Mu_pm  ))
        abline(a=0,b=1)
        plot(susiF.obj$fitted_func[[1]])
        if ( length(susiF.obj$fitted_func )>1){

          plot(susiF.obj$fitted_func[[2]])
        }


        par (mfrow=c(1,1))
      }else{
        par (mfrow=c(2,2))
        plot( log1p(Y ),  (Mu_pm  ))

        abline(a=0,b=1)
        trmse= mean((c(Eb_pm) - c(True_intensity))^2)


        plot(True_intensity,  (Eb_pm  ),
             main=  paste( " MSE" , trmse ))

        abline(a=0,b=1)
        plot(susiF.obj$fitted_func[[1]])
        if ( length(susiF.obj$fitted_func )>1){

          plot(susiF.obj$fitted_func[[2]])
        }



        par (mfrow=c(1,1))
      }

    }

  #  converged <- max_diff < tol
    iter <- iter + 1
  }

  if (!converged && verbose) {
    warning("Algorithm did not converge in ", maxit_outer, " iterations")
  }

  # ============================================================================
  # Finalize and return results
  # ============================================================================

    susiF.obj <- fsusieR::update_cal_pip(susiF.obj)


  # Compute fitted intensities
  fitted_log_intensity <- Eb_pm
  fitted_intensity <- exp(fitted_log_intensity)

  # Prepare output
  out <- list(
    Mu_pm = Mu_pm,
    Eb_pm = Eb_pm,
    sigma2 = sigma2,
    fitted = fitted_intensity[, idx_out],
    converged = converged,
    n_iter = iter - 1,
    elbo_hist=  elbo_hist
  )

  out$susiF.obj <- susiF.obj

  return(out)
}
