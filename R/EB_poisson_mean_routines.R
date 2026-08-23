
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

vga_pois_solution_is_valid = function(res, x, s, beta, sigma2,
                                      tol = 1e-4) {
  n = length(x)
  recycle_parameter = function(value) {
    if(length(value) == 1) rep(value, n) else value
  }
  s = recycle_parameter(s)
  beta = recycle_parameter(beta)
  sigma2 = recycle_parameter(sigma2)

  if(is.null(res) || is.null(res$m) || is.null(res$v) ||
     length(res$m) != n || length(res$v) != n ||
     length(s) != n || length(beta) != n || length(sigma2) != n ||
     any(!is.finite(res$m)) || any(!is.finite(res$v)) ||
     any(!is.finite(s)) || any(!is.finite(beta)) ||
     any(!is.finite(sigma2)) || any(s <= 0) || any(sigma2 <= 0) ||
     any(res$v <= 0) || any(res$v > sigma2 * (1 + tol))) {
    return(FALSE)
  }

  log_rate = log(s) + res$m + res$v/2
  if(any(!is.finite(log_rate)) ||
     any(log_rate > log(.Machine$double.xmax))) {
    return(FALSE)
  }
  rate = exp(log_rate)
  mean_score = x - rate - (res$m - beta)/sigma2
  variance_score = 1/res$v - rate - 1/sigma2
  mean_scale = 1 + x + rate + abs((res$m - beta)/sigma2)
  variance_scale = 1 + 1/res$v + rate + 1/sigma2

  all(is.finite(mean_score)) && all(is.finite(variance_score)) &&
    max(abs(mean_score)/mean_scale) <= tol &&
    max(abs(variance_score)/variance_scale) <= tol
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
  if(length(init_val) != n || length(sigma2) != n || length(beta) != n ||
     length(s) != n || any(!is.finite(init_val)) || any(!is.finite(x)) ||
     any(!is.finite(s)) || any(!is.finite(beta)) ||
     any(!is.finite(sigma2)) || any(x < 0) || any(s <= 0) ||
     any(sigma2 <= 0)) {
    stop('Invalid Poisson VGA inputs.', call. = FALSE)
  }
  if(method=='newton'){
    # Use safeguarded Newton first. A finite return value is not sufficient:
    # both variational score equations must also be satisfied.
    res = try(vga_pois_solver_Newton(init_val,x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent = TRUE)
    if(inherits(res,'try-error') ||
       !vga_pois_solution_is_valid(res,x,s,beta,sigma2,
                                   tol=max(10*tol,1e-8))){
      # If Newton failed, use bisection
      res = try(vga_pois_solver_bisection(x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent=TRUE)
      if(inherits(res,'try-error') ||
         !vga_pois_solution_is_valid(res,x,s,beta,sigma2,
                                     tol=max(10*tol,1e-8))){
        stop('Both Newton and Bisection methods failed to solve the Poisson VGA equations.',
             call. = FALSE)
      }else{
        return(res)
      }
    }else{
      return(res)
    }
  }else if(method=='bisection'){
    res = try(vga_pois_solver_bisection(x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent=TRUE)
    if(inherits(res,'try-error') ||
       !vga_pois_solution_is_valid(res,x,s,beta,sigma2,
                                   tol=max(10*tol,1e-8))){
      stop('Bisection failed to solve the Poisson VGA equations.',
           call. = FALSE)
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

h_m = function(m,x,s,beta,sigma2){
  const0 = sigma2*x+beta+1
  temp = const0-m
  log_rate = log(s)+m+sigma2/(2*temp)
  rate = rep(Inf,length(x))
  safe = is.finite(log_rate) & log_rate <= log(.Machine$double.xmax)
  rate[safe] = exp(log_rate[safe])
  # Negative mean score; this is monotone increasing in m.
  rate+(m-beta)/sigma2-x
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
  recycle_parameter = function(value) {
    if(length(value) == 1) rep(value,n) else value
  }
  s = recycle_parameter(s)
  beta = recycle_parameter(beta)
  sigma2 = recycle_parameter(sigma2)
  if(length(s) != n || length(beta) != n || length(sigma2) != n){
    stop('check Poisson VGA parameter lengths')
  }

  upper = beta+sigma2*x
  lower = pmin(beta,log1p(x)-log(s))-1
  lower_step = rep(1,n)
  lower_value = h_m(lower,x,s,beta,sigma2)
  for(k in 1:100){
    bad_lower = is.na(lower_value) | lower_value >= 0
    if(!any(bad_lower)) break
    lower[bad_lower] = lower[bad_lower]-lower_step[bad_lower]
    lower_step[bad_lower] = 2*lower_step[bad_lower]
    lower_value = h_m(lower,x,s,beta,sigma2)
  }
  if(any(is.na(lower_value)) || any(lower_value >= 0)){
    stop('vga_pois_solver_bisection: failed to bracket the root',
         call. = FALSE)
  }

  m = bisection(h_m,
                lower = lower,upper = upper,
                x=x,s=s,beta=beta,sigma2=sigma2,
                auto_adjust_interval = FALSE,
                maxiter=maxiter,tol=tol)
  v = sigma2/(sigma2*x+beta+1-m)
  return(list(m=m,v=v))
}

#'@export
vga_pois_solver_Newton = function(m,x,s,beta,sigma2,maxiter=1000,tol=1e-5){

  n = length(x)
  recycle_parameter = function(value) {
    if(length(value) == 1) rep(value, n) else value
  }
  s = recycle_parameter(s)
  beta = recycle_parameter(beta)
  sigma2 = recycle_parameter(sigma2)
  if(length(m) != n || length(s) != n || length(beta) != n ||
     length(sigma2) != n || any(!is.finite(m)) || any(!is.finite(x)) ||
     any(!is.finite(s)) || any(!is.finite(beta)) ||
     any(!is.finite(sigma2)) || any(x < 0) || any(s <= 0) ||
     any(sigma2 <= 0)) {
    stop('vga_pois_solver_Newton: invalid inputs', call. = FALSE)
  }

  const0 = sigma2*x+beta + 1
  upper = beta + sigma2*x
  if(any(!is.finite(const0)) || any(!is.finite(upper))) {
    stop('vga_pois_solver_Newton: non-finite search interval', call. = FALSE)
  }

  evaluate_score = function(value) {
    temp = const0-value
    log_rate = log(s) + value + sigma2/(2*temp)
    rate = rep(Inf, n)
    safe = is.finite(log_rate) & log_rate <= log(.Machine$double.xmax)
    rate[safe] = exp(log_rate[safe])
    score = x - rate - (value-beta)/sigma2
    gradient = -rate*(1+sigma2/(2*temp^2))-1/sigma2
    list(score=score,gradient=gradient,rate=rate)
  }

  # At upper = beta + sigma2*x the score is strictly negative. Move the lower
  # bound down until its score is positive, giving each coordinate a guaranteed
  # bracket around the unique root.
  lower = pmin(m,beta,log1p(x)-log(s))-1
  lower_step = rep(1,n)
  lower_score = evaluate_score(lower)$score
  for(k in 1:100){
    bad_lower = is.na(lower_score) | lower_score <= 0
    if(!any(bad_lower)) break
    lower[bad_lower] = lower[bad_lower]-lower_step[bad_lower]
    lower_step[bad_lower] = 2*lower_step[bad_lower]
    if(any(!is.finite(lower[bad_lower]))) {
      stop('vga_pois_solver_Newton: failed to bracket the root', call. = FALSE)
    }
    lower_score = evaluate_score(lower)$score
  }
  if(any(is.na(lower_score)) || any(lower_score <= 0)) {
    stop('vga_pois_solver_Newton: failed to bracket the root', call. = FALSE)
  }

  m = pmin(pmax(m,lower),upper)
  converged = rep(FALSE,n)

  for(i in 1:maxiter){
    evaluated = evaluate_score(m)
    f = evaluated$score
    f_grad = evaluated$gradient
    score_scale = 1+x+evaluated$rate+abs((m-beta)/sigma2)
    newly_converged = is.finite(f) & is.finite(score_scale) &
      abs(f) <= tol*score_scale
    converged = converged | newly_converged
    if(all(converged)) break

    active = !converged
    positive = active & is.finite(f) & f > 0
    lower[positive] = m[positive]
    upper[active & !positive] = m[active & !positive]

    candidate = m-f/f_grad
    midpoint = lower+(upper-lower)/2
    use_newton = active & is.finite(candidate) &
      candidate > lower & candidate < upper
    m[active] = midpoint[active]
    m[use_newton] = candidate[use_newton]
  }
  if(!all(converged)){
    stop('vga_pois_solver_Newton: method did not converge', call. = FALSE)
  }
  temp = const0-m
  v = sigma2/temp
  res = list(m=m,v=v)
  if(!vga_pois_solution_is_valid(res,x,s,beta,sigma2,
                                 tol=max(10*tol,1e-8))){
    stop("vga_pois_solver_Newton: invalid stationary point", call. = FALSE)
  }
  return(res)

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
#'  pois_mean_GP(x)
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
      #   temp = pois_mean_GP1(x[i],s[i],beta,sigma2,optim_method,m[i],v[i])
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
      obj[iter+1] = pois_mean_GP_obj(x,s,beta,sigma2,m,v)
      if((obj[iter+1] - obj[iter])<tol){
        obj = obj[1:(iter+1)]
        break
      }
    }

  }else{
    beta = prior_mean
    sigma2 = prior_var
    # for(i in 1:n){
    #   temp = pois_mean_GP1(x[i],s[i],prior_mean,prior_var,optim_method,m[i],v[i])
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
    obj = pois_mean_GP_obj(x,s,prior_mean,prior_var,m,v)

  }

  return(list(posterior = list(posteriorMean_latent = m,
                               posteriorVar_latent = v,
                               posteriorMean_mean = exp(m + v/2)),
              fitted_g = list(mean = beta, var=sigma2),
              obj_value=obj))

  #return(list(posteriorMean=m,priorMean=beta,priorVar=sigma2,posteriorVar=v,obj_value=obj))

}
#'calculate objective function


# Modify pois_mean_GP_opt_obj
pois_mean_GP_opt_obj = function(theta,x,s,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]

  # Cap the values to avoid Inf
  m = pmin(pmax(m, -20), 20)
  log_v = pmin(pmax(v, -20), 10)

  # Objective calculation with capped eta
  eta = m + exp(log_v)/2
  eta = pmin(eta, 25) # Exp(25) is large but finite

  return(-sum(x*m - s*exp(eta) - (m^2 + exp(log_v) - 2*m*beta)/2/sigma2 + log_v/2))
}
#'calculate gradient vector
pois_mean_GP_opt_obj_gradient = function(theta,x,s,beta,sigma2,n){
  #m = theta[1:n]
  #v = theta[(n+1):(2*n)]
  #g1 = -(x-s*exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  #g2 = -(-exp(v)/2*s*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  v = pmin(pmax(v, log(1e-6)), log(1e2))   # constrain log(v)
  sigma2 = max(sigma2, 1e-6)
  eta = m + exp(v)/2
  eta[eta > 20] <- 20
  g1 = -(x - s*exp(eta) - m/sigma2 + beta/sigma2)
  g2 = -(-exp(v)/2*s*exp(eta) - exp(v)/(2*sigma2) + 1/2)



  return(c(g1,g2))
}


pois_mean_GP_obj = function(x,s,beta,sigma2,m,v){
  return(sum(x*m-s*exp(m+v/2)-log(sigma2)/2-(m^2+v-2*m*beta+beta^2)/2/sigma2+log(v)/2))
}







#'@title Smooth over-dispersed Poisson sequence
#'@param x data vector
#'@param maxiter,tol max iteration and tolerance for stopping it.
#'@param Eb_init,sigma2_init initial values of smooth mean and nugget effect.
#'@examples
#' set.seed(12345)
#' n=2^9
#' sigma=0.5
#' mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
#' x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
#' fit = pois_smooth_split(x,maxiter=30)
#' plot(x,col='grey80')
#' lines(exp(fit$Eb))
#' fit$sigma2
#' plot(fit$obj)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(b_i,\sigma^2),}
#'\deqn{\b_i\sim g(.).}
#'@export

pois_smooth_split <- function(x,
                                                s = NULL,
                                                Eb_init = NULL,
                                                sigma2_init = NULL,
                                                est_sigma2 = TRUE,
                                                maxiter = 100,
                                                tol = 1e-5,
                                                filter.number = 1,
                                                family = 'DaubExPhase',
                                                verbose = FALSE,
                                                printevery = 10,
                                                ebnm_params = list(mode = 0),
                                                optim_method = 'L-BFGS-B',
                                                link = c("log", "log1p"))

  {
  link <- match.arg(link)
  n = length(x)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  if(is.null(Eb_init)){
    Eb = log(runmed(x/s,1 + 2 * min((n-1)%/% 2, ceiling(0.1*n)))+0.01)
  }else{
    Eb = Eb_init
  }
  if(is.null(sigma2_init)){
    sigma2 = var(log(x/s+0.01)-Eb)
  }else{
    sigma2 = sigma2_init
  }

  W = (t(wavethresh::GenW(n,filter.number,family)))[-1,]

  mu_pm = rep(0,n)
  mu_pv = rep(1/n,n)
  obj = -Inf
  if (link == "log") {
    obj_fn <- pois_mean_GP_opt_obj
    grad_fn <- pois_mean_GP_opt_obj_gradient
  } else if (link == "log1p") {
    obj_fn <- pois_mean_GP_log1p_opt_obj
    grad_fn <- pois_mean_GP_log1p_opt_obj_gradient
  }

  for(iter in 1:maxiter){
    # get m, s^2
    opt = optim(c(mu_pm, log(mu_pv)),
                fn = obj_fn,
                gr = grad_fn,
                lower = c(rep(-30, n), rep(-20, n)),
                upper = c(rep(30, n), rep(10, n)),
                x = x, s = s, beta = Eb,
                sigma2 = sigma2, n = n,
                method = optim_method)

    mu_pm = opt$par[1:n]
    mu_pv = exp(opt$par[(n+1):(2*n)])
    qb = smash_dwt(x=mu_pm,
                   sigma=sqrt(sigma2),
                   filter.number=filter.number,
                   family=family,
                   ebnm_params=ebnm_params,W=W)
    Eb = qb$mu.est
    Eb2 = qb$mu.est.var + Eb^2
    # get sigma2
    if(est_sigma2){
      sigma2 = mean(mu_pm^2+mu_pv+Eb2-2*mu_pm*Eb)
      sigma2 = max(sigma2, 1e-6)
    }


    # calc obj
    obj[iter+1] = pois_smooth_split_obj(x,s,mu_pm,mu_pv,Eb,Eb2,sigma2,qb$dKL)

    if(verbose){
      if(iter%%printevery==0){
        print(paste("Done iter",iter,"obj =",obj[iter+1]))
      }

    }

    if((obj[iter+1]-obj[iter])<tol){
      break
    }

  }
  return(list(Emean = exp(mu_pm+mu_pv/2),
              Vmean = exp(mu_pv-1)*exp(2*mu_pm+mu_pv),
              Emu=mu_pm,
              Vmu=mu_pv,
              Eb=Eb,
              Eb2=Eb2,
              sigma2=sigma2,
              obj=obj,
              H = qb$dKL + sum(log(2*pi*mu_pv)/2-log(2*pi*sigma2)/2-(mu_pm^2+mu_pv-2*mu_pm*Eb+Eb2)/2/sigma2)))
}

pois_smooth_split_obj = function(x,s,m,s2,Eb,Eb2,sigma2,KLb){
  return(sum(x*m-s*exp(m+s2/2)+log(s2)/2-log(sigma2)/2-(m^2+s2-2*m*Eb+Eb2)/2/sigma2)+KLb)
}



#'@title Empirical Bayes wavelet smoothing via DWT
#'@description Smooth homogeneous Gaussian data.
#'@param x data
#'@param sigma known standard error
#'@param filter.number,family wavelet family and filter number as in wavethresh package
#'@param ebnm_params a list of `ebnm` parameters
#'@param W the dwt matrix for calc posterior variance. Remove the first row which is all 1/sqrt(n)
#'@return a list of
#'  \item{mu.est:}{posterior mean}
#'  \item{mu.est.var:}{posterior variance}
#'  \item{loglik:}{log likelihood}
#'  \item{dKL:}{KL divergence between g(the prior) and q(the posterior)}
#'@import wavethresh
#'@import ebnm
#'@export
smash_dwt = function(x,sigma,filter.number=1,
                     family="DaubExPhase",
                     ebnm_params=list(mode=0),W=NULL){

  n = length(x)
  J = log(n,2)
  if(ceiling(J)!=floor(J)){
    stop('Length of x must be power of 2')
  }
  if(is.null(ebnm_params)){
    ebnm_params = ebnm_params_default()
  }else{
    temp = ebnm_params_default()
    for(i in 1:length(ebnm_params)){
      temp[[names(ebnm_params)[i]]] = ebnm_params[[i]]
    }
    ebnm_params = temp
  }
  tsum = sum(x)/sqrt(n)
  x.w = wavethresh::wd(x, filter.number = filter.number,
           family = family, type = "wavelet")

  data.var = sigma^2
  if(length(data.var==1)){
    data.var = rep(data.var,n)
  }

  if(is.null(W)){
    W = (t(GenW(n,filter.number,family)))[-1,]
  }

  if(length(sigma)==1){
    x.w.v = rep(sigma^2,n-1)
    tsum.var = sigma^2
  }else{
    x.w.v =  data.var
    tsum.var = x.w.v[1]
    x.w.v = x.w.v[-1]
  }

  dKL = 0
  loglik.scale = c()
  x.w.v.s = rep(0, 2^J-1)
  for (j in 0:(J - 1)) {
    x.pm = rep(0, 2^j)
    #index = (((J - 1) - j) * n + 1):((J - j) * n)
    index = (n-2^(j+1)+1):(n-2^j)
    x.w.j = accessD(x.w, j)
    x.w.v.j = x.w.v[index]
    ind.nnull = (x.w.v.j != 0)

    a = ebnm(x.w.j[ind.nnull],sqrt(x.w.v.j[ind.nnull]),
             mode=ebnm_params$mode,
             prior_family=ebnm_params$prior_family,
             scale = ebnm_params$scale,
             g_init = ebnm_params$g_init,
             fix_g = ebnm_params$fix_g,
             output = ebnm_params$output,
             optmethod = ebnm_params$optmethod)

    dKL = dKL + a$log_likelihood - Eloglik_GP(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),a$posterior$mean, a$posterior$mean^2+a$posterior$sd^2)
    x.pm[ind.nnull] = a$posterior$mean
    x.pm[!ind.nnull] = 0
    x.w = putD(x.w, j, x.pm)
    loglik.scale[j + 1] = a$log_likelihood
    x.w.v.s[index[ind.nnull]] = a$posterior$sd^2
    x.w.v.s[index[!ind.nnull]] = 0
  }
  mu.est = wr(x.w)
  loglik = sum(loglik.scale)
  #x.w.v.s = c(tsum.var,x.w.v.s)
  mu.est.var = colSums(W^2*x.w.v.s)
  return(list(mu.est=mu.est,mu.est.var=mu.est.var,loglik = loglik,dKL = dKL))
}



Eloglik_GP = function(x, s, Et, Et2) {
  # Deal with infinite SEs:
  idx = is.finite(s)
  x = x[idx]
  s = s[idx]
  Et = Et[idx]
  Et2 = Et2[idx]
  return(-0.5 * sum(log(2*pi*s^2) + (1/s^2) * (Et2 - 2*x*Et + x^2)))
}
output_default <- ebnm:::ebnm_output_default

ebnm_params_default = function(){
  return(list(prior_family='point_laplace',
              mode='estimate',
              scale = "estimate",
              g_init = NULL,
              fix_g = FALSE,
              output = output_default(),
              optmethod = NULL))
}
