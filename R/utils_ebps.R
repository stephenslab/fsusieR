#'@title empirical Bayes Poisson smoothing
#'@param x,s data vector and scaling factor. s can be a vector of the same length as x, or a scalar.
#'@param g_init a list of initial value of sigma2, and g_smooth. g_smooth is the initial prior g of the smoothing method. Can be NULL.
#'@param q_init a list of initial values of m, smooth. m is the posterior mean of mu, smooth the posterior mean of b. See the details below.
#'@param init_control See function ebps_init_control_default
#'@param general_control See function ebps_general_control_default
#'@param smooth_control See function ebps_smooth_control_default
#'@examples
#' set.seed(12345)
#' n=2^9
#' sigma=0.5
#' mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
#' x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
#' fit = ebps(x)
#' plot(x,col='grey80')
#' lines(fit$posterior$mean_smooth)
#' fit$sigma2
#' plot(fit$elbo_trace)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\lambda_i,}
#'\deqn{\lambda_i = \exp(\mu_i)),}
#'\deqn{\mu_i\sim N(b_i,\sigma^2),}
#'\deqn{\b_i\sim g(.).}
#'
#'The \code{init_control} argument is a list in which any of the following
#'named components will override the default algorithm settings (as
#'defined by \code{ebps_init_control_default}):
#'
#'\describe{
#'\item{\code{m_init_method}}{'vga' or 'smash_poi'}
#'}
#'
#'The \code{general_control} argument is a list in which any of the following
#'named components will override the default algorithm settings (as
#'defined by \code{ebps_general_control_default}):
#'
#'\describe{
#'\item{\code{est_sigma2}}{whether estiamte sigma2 or fix it}
#'\item{\code{maxiter}}{max iteration of the main algorithm, default is 100}
#'\item{\code{maxiter_vga}}{max iteration of the vga step}
#'\item{\code{vga_tol}}{tolerance for vga step stopping}
#'\item{\code{verbose}}{print progress?}
#'\item{\code{tol}}{tolerance for stopping the main algorithm}
#'\item{\code{convergence_criteria}}{'objabs' or 'nugabs'}
#'\item{\code{make_power_of_2}}{'reflect' or 'extend'}
#'\item{\code{plot_updates}}{internal use only}
#'}
#'
#'The \code{smooth_control} argument is a list in which any of the following
#'named components will override the default algorithm settings (as
#'defined by \code{ebps_smooth_control_default}):
#'
#'\describe{
#'\item{\code{wave_trans}}{'dwt' or 'ndwt'}
#'\item{\code{ndwt_method}}{'smash' or 'ti.thresh'}
#'\item{\code{ebnm_params}}{parameters for ebnm used in wavelet smoothing}
#'\item{\code{warmstart}}{init posterior using last iteration's results}
#'\item{\code{W}}{DWT matrix for non-haar wavelet basis}
#'}
#'@import wavethresh
#'@import smashr
#'@export




ebps = function(x,
                s = NULL,
                g_init = NULL,
                q_init = NULL,
                init_control = list(),
                general_control = list(),
                smooth_control = list()
){
  t_start = Sys.time()
  n_orig = length(x)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n_orig)
  }

  init_controls = modifyList(ebps_init_control_default(),init_control,keep.null = TRUE)
  general_controls = modifyList(ebps_general_control_default(),general_control,keep.null = TRUE)
  smooth_controls = modifyList(ebps_smooth_control_default(),smooth_control,keep.null = TRUE)

  init_val = ebps_init(x,s,
                       general_controls$make_power_of_2,
                       g_init,
                       q_init,
                       init_controls$m_init_method)

  sigma2 = init_val$g_init$sigma2
  m = init_val$q_init$m
  x = init_val$x
  s = init_val$s
  idx = init_val$idx
  Eb = init_val$q_init$smooth
  if(is.null(Eb)){
    Eb = rep(mean(m),length(x))
  }
  n = length(x)


  const = sum(lfactorial(x))
  v = rep(sigma2/2,length(x))

  if(smooth_controls$wave_trans=='ndwt'&smooth_controls$ndwt_method=='ti.thresh' | smooth_controls$robust){
    general_controls$convergence_criteria = 'nugabs'
  }
  if(smooth_controls$wave_trans=='dwt'&is.null(smooth_controls$W)&(smooth_controls$filter.number != 1 | smooth_controls$family != 'DaubExPhase')){
    smooth_controls$W = (t(GenW(n,filter.number,family)))[-1,]
  }


  obj = -Inf
  s_update = list(Eb=Eb,
                  qb = list(fitted_g = init_val$g_init$g_smooth))

  Eb_old = Inf
  sigma2_trace = c(sigma2)

  for(iter in 1:general_controls$maxiter){

    s_update = ebps_smooth_update(m             = m,
                                  sigma2        = sigma2,
                                  filter.number = smooth_controls$filter.number,
                                  family        = smooth_controls$family,
                                  wave_trans    = smooth_controls$wave_trans,
                                  ndwt_method   = smooth_controls$ndwt_method,
                                  warmstart     = smooth_controls$warmstart,
                                  ebnm_params   = smooth_controls$ebnm_params,
                                  qb            = s_update$qb,
                                  W             = smooth_controls$W,
                                  robust        = smooth_controls$robust,
                                  Eb            = s_update$Eb)









    if(general_controls$plot_updates){
      plot(m,col='grey80')
      lines(s_update$Eb)
    }

    # opt = vga_pois_solver(m,x,s,Eb,sigma2,tol=vga_tol)
    if(general_controls$maxiter_vga==1){
      opt = vga_pois_solver_newton_1iter(m,v,x,s,s_update$Eb+s_update$E_out,sigma2)
    }else{
      opt = vga_pois_solver(m,x,s,s_update$Eb+s_update$E_out,sigma2,tol=general_controls$vga_tol,maxiter = general_controls$maxiter_vga)
    }

    m = opt$m
    v = opt$v

    # get sigma2
    if(general_controls$est_sigma2){
      sigma2_new = mean(m^2+v+s_update$Eb2+s_update$E_out2+2*s_update$Eb*s_update$E_out-2*m*(s_update$Eb+s_update$E_out))
      sigma2_trace = c(sigma2_trace,sigma2_new)
      if(general_controls$convergence_criteria=='nugabs'){
        if(general_controls$verbose){
          if(iter%%general_controls$printevery==0){
            print(paste("Done iter",iter,"sigma2 =",sigma2_new))
          }
        }
        if(abs(sigma2_new-sigma2)<general_controls$tol){
          sigma2 = sigma2_new
          break
        }
      }
      #print(sigma2_new)
      sigma2 = sigma2_new
    }else{
      if(general_controls$convergence_criteria=='nugabs'){
        if(sqrt(mean((s_update$Eb-s_update$Eb_old)^2))<general_controls$tol){
          break
        }
        Eb_old = s_update$Eb
      }
    }


    # calc obj
    if(general_controls$convergence_criteria=='objabs'){
      # check this  ----
      obj[iter+1] = pois_smooth_split_obj(x=x, # here is the potential pb
                                          s=s,
                                          m=m,
                                          s2=v,
                                          Eb=s_update$Eb,
                                          Eb2=s_update$Eb2,
                                          sigma2 = sigma2,
                                          KLb =s_update$qb$dKL)#,const)
      if(general_controls$verbose){
        if(iter%%general_controls$printevery==0){
          print(paste("Done iter",iter,"obj =",obj[iter+1]))
        }
      }

      if((obj[iter+1]-obj[iter])/n <general_controls$tol){
        break
      }
    }

  }
  t_end = Sys.time()
  if(smooth_controls$wave_trans=='dwt'){
    return(list(posterior=list(mean=exp(m+v/2)[idx],
                               mean_log = m[idx],
                               mean_smooth = exp(s_update$Eb)[idx],
                               mean_log_smooth=s_update$Eb[idx],
                               var_log = v[idx],
                               var_log_smooth = (s_update$Eb2-s_update$Eb^2)[idx]),
                fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace,g = s_update$qb$fitted_g),
                elbo=obj[length(obj)]/n*n_orig,
                elbo_trace = obj/n*n_orig,
                H = (s_update$qb$dKL + sum(log(2*pi*v)/2-log(2*pi*sigma2)/2-(m^2+v-2*m*s_update$Eb+s_update$Eb2)/2/sigma2))/n*n_orig,
                log_likelihood = obj[length(obj)]/n*n_orig,
                run_time = difftime(t_end,t_start,units='secs')))
  }else{
    if(smooth_controls$ndwt_method=='smash'){
      return(list(posterior=list(mean = exp(m+v/2)[idx],
                                 mean_log = m[idx],
                                 mean_smooth = exp(s_update$Eb)[idx],
                                 mean_log_smooth=s_update$Eb[idx]),
                  log_likelihood = obj[length(obj)]/n*n_orig,
                  elbo_trace = obj/n*n_orig,
                  fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace),
                  run_time = difftime(t_end,t_start,units='secs')))
    }else{
      return(list(posterior=list(mean = exp(m+v/2)[idx],
                                 mean_log = m[idx],
                                 mean_smooth = exp(s_update$Eb)[idx],
                                 mean_log_smooth=s_update$Eb[idx]),
                  log_likelihood = NULL,
                  obj_trace = obj/n*n_orig,
                  fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace),
                  run_time = difftime(t_end,t_start,units='secs')))
    }
  }
}

ebps_smooth_update = function(m,sigma2,
                              filter.number,
                              family,
                              wave_trans,
                              ndwt_method,
                              warmstart,
                              ebnm_params,
                              qb,
                              W,
                              robust,
                              Eb){
  #browser()
  if(robust){
    #print(range(m-Eb))
    res = ebnm(m-Eb,sqrt(sigma2),
               prior_family = 'point_laplace',
               mode = 0,
               scale =  "estimate",
               g_init = NULL,
               fix_g = FALSE)
    E_out = res$posterior$mean
    E_out2 = res$posterior$sd^2 + E_out^2
  }else{
    E_out = 0
    E_out2 = 0
  }
  m = m - E_out

  if(wave_trans=='dwt'){
    if(warmstart){
      qb = suppressWarnings( smash_dwt(m,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=list(g_init=qb$fitted_g),W=W))
    }else{
      qb = smash_dwt(m,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=ebnm_params,W=W)
    }
    Eb = qb$mu.est
    Eb2 = Eb^2+qb$mu.est.var
    qb$dKL = qb$loglik - Eloglik(m,rep(sqrt(sigma2),length(m)),Eb,Eb2)
  }
  if(wave_trans=='ndwt'){
    if(ndwt_method=='smash'){
      qb = smashr::smash.gaus(m,sqrt(sigma2),filter.number=filter.number,family=family,
                              ebnm_param=ebnm_params,
                              post.var = TRUE,
                              return.loglr = T)
      Eb = qb$mu.est
      Eb2 = Eb^2+qb$mu.est.var
      qb$dKL = qb$loglik - Eloglik(m,rep(sqrt(sigma2),length(m)),Eb,Eb2)
    }
    if(ndwt_method=='ti.thresh'){
      Eb = smashr::ti.thresh(m,sqrt(sigma2),filter.number=filter.number,family=family)

      Eb2 = Eb^2
    }
  }
  return(list(Eb=Eb,
              Eb2=Eb2,
              E_out=E_out,
              E_out2=E_out2,
              qb=qb))
}

#'@title Default parameters for ebps initialization
#'@param m_init_method 'vga' or 'smash_poi'
#'@export
ebps_init_control_default = function(){
  return(list(m_init_method = 'vga'))
}

#'@title Default parameters for ebps
#'@export
ebps_general_control_default = function(){
  return(list(est_sigma2 = TRUE,
              maxiter = 100,
              maxiter_vga = 100,
              vga_tol = 1e-5,
              verbose=FALSE,
              tol=1e-5,
              printevery = 10,
              convergence_criteria = 'objabs',
              make_power_of_2 = 'reflect',
              plot_updates = FALSE))
}

#'@title Default parameters for ebps smoothing function
#'@param ndwt_method 'smash' or 'ti.thresh'
#'@export
ebps_smooth_control_default = function(){
  return(list(filter.number = 1,
              family = 'DaubExPhase',
              wave_trans='dwt',
              ndwt_method='ti.thresh',
              ebnm_params=list(),
              warmstart=TRUE,
              W=NULL,
              robust = FALSE
  ))
}


# pois_smooth_split_obj = function(x,s,m,s2,Eb,Eb2,sigma2,KLb,const){
#   return(sum(x*m-s*exp(m+s2/2)+log(s2)/2-log(sigma2)/2-(m^2+s2-2*m*Eb+Eb2)/2/sigma2)+KLb-const)
# }
#
# extend = function(x){
#   n = length(x)
#   J = log2(n)
#   if ((J%%1) == 0) {
#     return(list(x = x, idx = 1:n))
#   }else {
#     n.ext = 2^ceiling(J)
#     lnum = round((n.ext - n)/2)
#     rnum = n.ext - n - lnum
#     if (lnum == 0) {
#       x.lmir = NULL
#     }else {
#       x.lmir = x[lnum:1]
#     }
#     if (rnum == 0) {
#       x.rmir = NULL
#     }else {
#       x.rmir = x[n:(n - rnum + 1)]
#     }
#     x.ini = c(x.lmir, x, x.rmir)
#     return(list(x = x.ini, idx = (lnum + 1):(lnum + n)))
#   }
#
# }

ebnm_params_default <- function(){
  return(list(prior_family='point_laplace',
              mode='estimate',
              scale = "estimate",
              g_init = NULL,
              fix_g = FALSE,
              output = output_default(),
              optmethod = NULL))
}
output_default=  ebnm:::ebnm_output_default

