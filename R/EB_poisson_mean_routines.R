#copied from https://github.com/DongyueXie/vebpm/

#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian prior, Gaussian posterior in Poisson mean problem.
#'@param x data vector
#'@param s scaling vector
#'@param g_init a list of mean, and var. Can be NULL for both parameters.
#'@param fix_g Whether fix g at g_init. If only fix either mean, or var, fix_g can be a length 2 boolean vector.
#'@param q_init a list of init value of m_init(posterior mean) and v_init(posterior var).
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@param conv_type convergence criteria, default to be elbo
#'@param return_sigma2_trace whether return the trace of sigma2 estiamtes
#'@return a list of
#'  \item{posterior:}{mean_log/var_log is the posterior mean/var of mu, mean is the posterior of exp(mu)}
#'  \item{fitted_g:}{estimated prior}
#'  \item{obj_value:}{objective function values}
#'@examples
#'\dontrun{
#'n = 1000
#'mu = rnorm(n)
#'x = rpois(n,exp(mu))
#'ebpm_normal(x)
#'}
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(\beta,\sigma^2).}
#'@export


ebpm_normal = function(x,
                       s = NULL,
                       g_init = NULL,
                       fix_g = FALSE,
                       q_init = NULL,
                       maxiter = 20,
                       tol = 1e-5,
                       vga_tol=1e-5,
                       conv_type='sigma2abs',
                       return_sigma2_trace=FALSE){
  
  # init the posterior mean and variance?
  n = length(x)
  
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  
  if(is.null(q_init)){
    m = log(x/s+1)
    v = rep(1/n,n)
  }else{
    if(is.null(q_init$m_init)){
      m = log(1+x/s)
    }else{
      m = q_init$m_init
    }
    if(is.null(q_init$v_init)){
      v = rep(1/n,n)
    }else{
      v = q_init$v_init
    }
  }
  if(length(v)==1){
    v = rep(v,n)
  }
  
  const = sum((x-1)*log(s)) - sum(lfactorial(x))
  #
  t_start = Sys.time()
  
  if(is.null(g_init)){
    prior_mean = NULL
    prior_var = NULL
  }else{
    prior_mean = g_init$mean
    prior_var = g_init$var
  }
  
  if(length(fix_g)==1){
    est_prior_mean = !fix_g
    est_prior_var = !fix_g
  }else if(length(fix_g)==2){
    est_prior_mean = !fix_g[1]
    est_prior_var = !fix_g[2]
  }else{
    stop('fix_g can be either length 1 or 2')
  }
  
  sigma2_trace = prior_var
  if(est_prior_mean | est_prior_var){
    
    if(is.null(prior_mean)){
      est_prior_mean = TRUE
      beta = mean(m)
    }else{
      beta = prior_mean
    }
    if(is.null(prior_var)){
      est_prior_var=TRUE
      sigma2 = mean(m^2+v-2*m*beta+beta^2)
    }else{
      sigma2=prior_var
    }
    
    obj = rep(0,maxiter+1)
    obj[1] = -Inf
    sigma2_trace = sigma2
    for(iter in 1:maxiter){
      sigma2_old = sigma2
      m = vga_pois_solver(m,x,s,beta,sigma2,tol=vga_tol)
      v =  m$v
      m = m$m
      
      if(est_prior_mean){
        beta = mean(m)
      }
      if(est_prior_var){
        sigma2 = mean(m^2+v-2*m*beta+beta^2)
        if(return_sigma2_trace){
          sigma2_trace[iter+1] = sigma2
        }
      }
      
      if(conv_type=='elbo'){
        obj[iter+1] = ebpm_normal_obj(x,s,beta,sigma2,m,v,const)
        if((obj[iter+1] - obj[iter])/n <tol){
          obj = obj[1:(iter+1)]
          if((obj[iter+1]-obj[iter])<0){
            warning('An iteration decreases ELBO. This is likely due to numerical issues.')
          }
          break
        }
      }
      if(conv_type=='sigma2abs'){
        obj[iter+1] = abs(sigma2-sigma2_old)
        if(obj[iter+1] <tol){
          obj = obj[1:(iter+1)]
          break
        }
      }
      
    }
    
  }else{
    beta = prior_mean
    sigma2 = prior_var
    m = vga_pois_solver(m,x,s,beta,sigma2,tol=vga_tol)
    v = m$v
    m = m$m
    obj = ebpm_normal_obj(x,s,prior_mean,prior_var,m,v,const)
    
  }
  t_end = Sys.time()
  
  return(list(posterior = list(mean_log = m,
                               var_log = v,
                               mean = exp(m + v/2)),
              fitted_g = list(mean = beta, var=sigma2),
              elbo=ebpm_normal_obj(x,s,beta,sigma2,m,v,const),
              obj_trace = obj,
              sigma2_trace=sigma2_trace,
              run_time = difftime(t_end,t_start,units='secs')))
  
}


ebpm_normal_obj = function(x,s,beta,sigma2,m,v,const){
  return(sum(x*m-s*exp(m+v/2)-log(sigma2)/2-(m^2+v-2*m*beta+beta^2)/2/sigma2+log(v)/2)+const)
}
  






#'@title Optimize vga poisson problem
#'@description This function tries Newton's method. If not working, then use bisection.
#'@param init_val initial value for posterior mean
#'@param x,s data and scale factor
#'@param beta,sigma2 prior mean and variance. Their length should be equal to n=length(x)
#'@export
vga_pois_solver = function(init_val,x,s,beta,sigma2,maxiter=1000,tol=1e-5,method = 'newton'){
  
  n = length(x)
  if(length(sigma2)==1){
    sigma2 = rep(sigma2,n)
  }
  if(length(beta)==1){
    beta = rep(beta,n)
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  if(method=='newton'){
    # use Newton's method fist
    res = try(vga_pois_solver_Newton(init_val,x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent = TRUE)
    if(class(res)=='try-error'){
      # If Newton failed, use bisection
      res = try(vga_pois_solver_bisection(x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent=TRUE)
      if(class(res)=='try-error'){
        # If bisection also failed, return initial  values with a warning.
        warnings('Both Newton and Bisection methods failed. Returning initial values.')
        return(list(m = init_val,v = sigma2/(sigma2*x+beta+1-init_val)))
      }else{
        return(res)
      }
    }else{
      return(res)
    }
  }else if(method=='bisection'){
    res = try(vga_pois_solver_bisection(x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent=TRUE)
    if(class(res)=='try-error'){
      return(list(m = init_val,v = sigma2/(sigma2*x+beta+1-init_val)))
    }else{
      return(res)
    }
  }else{
    stop('Only Newton and bisection are supported.')
  }
  
  
}

#'@title Optimize vga poisson problem 1 iteration.
#'@description This function only performs vga for 1 iteration
#'@export
vga_pois_solver_newton_1iter = function(m,v,y,s,beta,sigma2){
  if(is.null(v)){
    v = sigma2/2
  }
  #for(i in 1:maxiter){
  # newton for M
  temp = 1/(s*exp(m+v/2))
  m = m - (y*temp-1-(m-beta)/sigma2*temp)/(-1-1/sigma2*temp)
  # newton for V
  temp = 1/(s*exp(m+v/2)*v/2)
  v = v/exp((-1-v/2/sigma2*temp + 0.5*temp)/(-(1 + v/2) - v/2/sigma2*temp))
  #}
  return(list(m=m,v=v))
}

h_v = function(v,x,s,beta,sigma2){
  val = (-s*exp(sigma2*x+beta+1-sigma2/v + v/2) - 1/sigma2 + 1/v)
  return(-val)
}

# h_m = function(v,x,s,beta,sigma2){
#   m = sigma2*x + beta + 1 - sigma2/v
#   val = x - s*exp(m+v/2) - (m-beta)/sigma2
#   return(-val)
# }


# vga_pois_solver_m = function(init_val,x,s,beta,sigma2,maxiter=1000,tol=1e-8){
#
#   n = length(x)
#   if(length(sigma2)==1){
#     upper = rep(sigma2,n)
#   }else if(length(sigma2)==n){
#     upper = sigma2
#   }else{
#     stop('check length of sigma2')
#   }
#   return(vga_optimize_m(init_val,x,s,beta,sigma2))
#
# }

#'@export
vga_pois_solver_bisection = function(x,s,beta,sigma2,maxiter=1000,tol=1e-5){
  n = length(x)
  if(length(sigma2)==1){
    upper = rep(sigma2,n)
  }else if(length(sigma2)==n){
    upper = sigma2
  }else{
    stop('check length of sigma2')
  }
  v = bisection(h_v,
                lower = rep(0,n),upper = upper,
                x=x,s=s,beta=beta,sigma2=sigma2,
                auto_adjust_interval = FALSE,
                maxiter=maxiter,tol=tol)
  m = sigma2*x + beta + 1 - sigma2/v
  return(list(m=m,v=v))
}

#'@export
vga_pois_solver_Newton = function(m,x,s,beta,sigma2,maxiter=1000,tol=1e-5){
  
  const0 = sigma2*x+beta + 1
  const1 = 1/sigma2
  const2 = sigma2/2
  const3 = beta/sigma2
  
  # make sure m < sigma2*x+beta
  m = pmin(m,const0-1)
  # idx = (m>(const0-1))
  # if(sum(idx)>0){
  #   m[idx] =suppressWarnings(vga_pois_solver_bisection(x[idx],s[idx],beta[idx],sigma2[idx],maxiter = 10)$m)
  # }
  
  
  for(i in 1:maxiter){
    
    temp = (const0-m)
    sexp = s*exp(m+const2/temp)
    # f = x - sexp - (m-beta)/sigma2
    f = x - sexp - m*const1 + const3
    if(max(abs(f))<tol){
      break
    }
    # f_grad = -sexp*(1+const2/temp^2)-const1
    m = m - f/(-sexp*(1+const2/temp^2)-const1)
  }
  if(i>=maxiter){
    warnings('Newton method not converged yet.')
  }
  return(list(m=m,v=sigma2/temp))
  
}
 


#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian prior, Gaussian posterior in Poisson mean problem.
#'@param x data vector
#'@param s scaling vector
#'@param w prior weights
#'@param prior_mean prior mean
#'@param prior_var prior variance
#'@param optim_method optimization method in `optim` function
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{posteriorMean:}{posterior mean}
#'  \item{posteriorVar:}{posterior variance}
#'  \item{obj_value:}{objective function values}
#'  \item{prior_mean:}{prior mean}
#'  \item{prior_var:}{prior variance}
#'  @example
#'  n = 10000
#'  mu = rnorm(n)
#'  x = rpois(n,exp(mu))
#'  pois_mean_GG(x)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(\beta,\sigma^2).}
#'@export
pois_mean_GP = function(x,
                        s = NULL,
                        prior_mean = NULL,
                        prior_var=NULL,
                        optim_method = 'L-BFGS-B',
                        maxiter = 1000,
                        tol = 1e-5){
  
  # init the posterior mean and variance?
  n = length(x)
  m = log(x+0.1)
  v = rep(1/sqrt(n),n)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  #
  if(is.null(prior_mean) | is.null(prior_var)){
    
    if(is.null(prior_mean)){
      est_beta = TRUE
    }else{
      est_beta = FALSE
      beta = prior_mean
    }
    if(is.null(prior_var)){
      est_sigma2=TRUE
    }else{
      est_sigma2 = FALSE
      sigma2=prior_var
    }
    
    obj = rep(0,maxiter+1)
    obj[1] = -Inf
    for(iter in 1:maxiter){
      if(est_beta){
        beta = mean(m)
      }
      if(est_sigma2){
        sigma2 = mean(m^2+v-2*m*beta+beta^2)
      }
      # for(i in 1:n){
      #   temp = pois_mean_GG1(x[i],s[i],beta,sigma2,optim_method,m[i],v[i])
      #   m[i] = temp$m
      #   v[i] = temp$v
      # }
      opt = optim(c(m,log(v)),
                  fn = pois_mean_GP_opt_obj,
                  gr = pois_mean_GP_opt_obj_gradient,
                  x=x,
                  s=s,
                  beta=beta,
                  sigma2=sigma2,
                  n=n,
                  method = optim_method)
      m = opt$par[1:n]
      v = exp(opt$par[(n+1):(2*n)])
      obj[iter+1] = pois_mean_GG_obj(x,s,beta,sigma2,m,v)
      if((obj[iter+1] - obj[iter])<tol){
        obj = obj[1:(iter+1)]
        break
      }
    }
    
  }else{
    beta = prior_mean
    sigma2 = prior_var
    # for(i in 1:n){
    #   temp = pois_mean_GG1(x[i],s[i],prior_mean,prior_var,optim_method,m[i],v[i])
    #   m[i] = temp$m
    #   v[i] = temp$v
    # }
    opt = optim(c(m,log(v)),
                fn = pois_mean_GP_opt_obj,
                gr = pois_mean_GP_opt_obj_gradient,
                x=x,
                s=s,
                beta=beta,
                sigma2=sigma2,
                n=n,
                method = optim_method)
    m = opt$par[1:n]
    v = exp(opt$par[(n+1):(2*n)])
    obj = pois_mean_GG_obj(x,s,prior_mean,prior_var,m,v)
    
  }
  
  return(list(posterior = list(posteriorMean_latent = m,
                               posteriorVar_latent = v,
                               posteriorMean_mean = exp(m + v/2)),
              fitted_g = list(mean = beta, var=sigma2),
              obj_value=obj))
  
  #return(list(posteriorMean=m,priorMean=beta,priorVar=sigma2,posteriorVar=v,obj_value=obj))
  
}
#'calculate objective function
pois_mean_GP_opt_obj = function(theta,x,s,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  return(-sum(x*m-s*exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+v/2))
}
#'calculate gradient vector
pois_mean_GP_opt_obj_gradient = function(theta,x,s,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  g1 = -(x-s*exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  g2 = -(-exp(v)/2*s*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  return(c(g1,g2))
}


pois_mean_GG_obj = function(x,s,beta,sigma2,m,v){
  return(sum(x*m-s*exp(m+v/2)-log(sigma2)/2-(m^2+v-2*m*beta+beta^2)/2/sigma2+log(v)/2))
}
