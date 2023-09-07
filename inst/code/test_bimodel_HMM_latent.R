fit_hmm <- function (x,sd, halfK=20,mult=2, smooth=FALSE,thresh=0.001){

  K = halfK+1
  sd=  sd
  X <-x

  pos <- seq(  0, 1 ,   length.out=K)

  #define the mean states
  mu <- (pos^(1/mult))*1.5*max(abs(X)) # put 0 state at



  P <- diag(0.99,K) #this ensure that the HMM can only "transit via null state"
  P[1,-1] <- 0.1
  P[-1,1] <- 0.1


  pi = rep( 1/length(mu), length(mu)) #same initial guess

  emit = function(k,x,t){
    0.5*dnorm(x,mean=mu[k],sd=sd[t]  )+0.5*dnorm(x,mean=-mu[k],sd=sd[t]  )
  }


  alpha = matrix(nrow = length(X),ncol=K)
  for(k in 1:K){
    alpha[1,k] = pi[k] * emit(k=k,x=X[1],t=1)
  }



  # Initialize alpha[1,]
  for(k in 1:K){
    alpha[1,k] = pi[k] * emit(k=k,x=X[1],t=1)
  }

  # Forward algorithm
  for(t in 1:(length(X)-1)){
    m = alpha[t,] %*% P
    for(k in 1:K){
      alpha[t+1,k] = m[k]*emit(k=k,x=X[t+1], t= t+1 )
    }
  }

  beta = matrix(nrow =  length(X),ncol=K)

  # Initialize beta
  for(k in 1:K){
    beta[ length(X),k] = 1
  }

  # Backwards algorithm
  for(t in ( length(X)-1):1){
    for(k in 1:K){
      beta[t,k] = sum(beta[t+1,]*P[k,]*emit(1:K,X[t+1],t=t+1))
    }
  }
  ab = alpha*beta
  prob = ab/rowSums(ab)
  #plot(apply(prob[,-1],1, sum), type='l')


  idx_comp <- which( apply(prob, 2, mean) >thresh )
  ash_obj_pos <- list()
  ash_obj_neg <- list()
  x_post <- 0*x
  for ( i in 1:length(idx_comp)){
    mu_ash_pos <-  mu[idx_comp[i]]

    mu_ash_neg <-  -mu[idx_comp[i]]
    weight <- prob[,idx_comp[i]]

    ash_obj_pos[[i]]  <- ash(x,sd,weight=weight,mode=mu_ash_pos)
    ash_obj_neg[[i]]  <- ash(x,sd,weight=weight,mode=mu_ash_neg)
    assig_latent_pos <-   dnorm( x, mean=mu_ash_pos,sd= sd)/(
      dnorm( x, mean=mu_ash_pos,sd= sd)+
        dnorm( x, mean=mu_ash_neg,sd= sd))
    x_post <-  x_post +weight*(assig_latent_pos*ash_obj_neg[[i]]$result$PosteriorMean+
                                 (1-assig_latent_pos)*ash_obj_pos[[i]]$result$PosteriorMean)

  }





  out <- list( prob =prob,
               fitted= prob%*%mu,
               x_post = x_post,
               mu= mu)




  return( out)
}
