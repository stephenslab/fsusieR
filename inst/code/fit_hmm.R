fit_hmm <- function (x,sd, halfK=20,mult=2){

  K = 2*halfK-1
  sd=  sd
  X <-x

  pos <- seq(  0, 1 ,   length.out=halfK)

  #define the mean states
  mu <- (pos^(1/mult))*max(abs(X)) # put 0 state at
  #the firstplace
  mu <- c(mu, -mu[-1] )



  P <- diag(0.99,K) #this ensure that the HMM can only "transit via null state"
  P[1,-1] <- 0.01
  P[-1,1] <- 0.01

  pi = rep( 1/length(mu), length(mu)) #same initial guess

  emit = function(k,x,t){
    dnorm(x,mean=mu[k],sd=sd[t])
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




  out <- list( prob =prob,
               fitted= prob%*%mu,
               mu= mu)




  return( out)
}
