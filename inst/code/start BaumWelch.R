


### basic forward back + P update.
x  = res$Bhat[idx[2  ],]
sd = res$Shat[idx[2],]
halfK=50
mult=2
K = 2*halfK-1
sd=  sd
X <-x
thresh=1e-5
pos <- seq(  0, 1 ,   length.out=halfK)

#define the mean states
mu <- (pos^(1/mult))*1.5*max(abs(X)) # put 0 state at
#the firstplace
mu <- c(mu, -mu[-1] )



P <- diag(0.5 ,K) #this ensure that the HMM can only "transit via null state"
P[1,-1] <- 0.5
P[-1,1] <- 0.5




pi = rep( 1/length(mu), length(mu)) #same initial guess

emit = function(k,x,t){
  dnorm(x,mean=mu[k],sd=sd[t]  )
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
if ( !(1%in% idx_comp) ){ #ensure 0 is in the model

  idx_comp<- c(1, idx_comp)
}
ash_obj <- list()
x_post <- 0*x
for ( i in 1:length(idx_comp)){
  mu_ash <-  mu[idx_comp[i]]

  weight <- prob[,idx_comp[i]]

  ash_obj[[i]]  <- ash(x,sd,weight=weight,mode=mu_ash)
  x_post <-  x_post +weight*ash_obj[[i]]$result$PosteriorMean

}


list_z_nz <- list() # transition from 0 to non zero state
list_nz_z <- list() # transition from  non zero state  to 0
list_self <- list()# stay in the same state
list_tt   <- list()


for ( t in 1:(length(X)-1)){
  tt_z_nz <-c()
  tt_nz_z <-c()
  tt_self <-c()
  for(  j in idx_comp[-1]){
    tt_z_nz <-c(tt_z_nz,
                (alpha[t, 1]*P[1,j]*beta[t+1,j ]*emit(k=j,x=X[t+1], t= t+1 )  )  )
    # transition from 0 to non zero state
  }
  for(  j in idx_comp[-1]){
    tt_nz_z  <-c(tt_nz_z , (alpha[t, j]*P[1,j]*beta[  t+1,1 ]*emit(k=1,x=X[t+1], t= t+1 )  )  )
    # transition from  non zero state  to 0
  }
  for(  j in idx_comp){
    tt_self <-c(tt_self, (alpha[t, j]*P[j,j]*beta[  t+1,j ]*emit(k=j,x=X[t+1], t= t+1 )  )  )
    # stay in the same state
  }


  n_c <-  (sum(tt_z_nz) +sum(tt_nz_z )+ sum(tt_self)  )
  list_z_nz[[t]] <- tt_z_nz/ n_c  # transition from 0 to non zero state
  list_nz_z[[t]] <- tt_nz_z/  n_c   # transition from  non zero state  to 0
  list_self[[t]] <- tt_self/  n_c# stay in the same state

  list_tt[[t]] <-  c(tt_z_nz/ n_c,
                     tt_nz_z/  n_c ,
                     tt_self/  n_c )




}
expect_number_obs_state <- apply(prob[-nrow(  prob),  idx_comp ],2,sum)

image (t(do.call(cbind,list_tt)))
plot(apply(prob[,-1],1, sum), type='l')

transition <- apply( do.call(cbind,list_tt),1,function(j)  mean(j, na.rm=TRUE))# add log sum exp tric


#Formula from Baum Welch update Wikipedi pagfe

diag_P <- apply(do.call( cbind,list_self),1 ,sum) /expect_number_obs_state
z_nz  <- apply(do.call( cbind,list_z_nz),1 ,sum) /expect_number_obs_state[1]
nz_z  <- apply(do.call( cbind,list_nz_z ),1 ,sum) /rev( expect_number_obs_state[-1])
nz_z  <- apply(do.call( cbind,list_nz_z ),1 ,sum) /  expect_number_obs_state[-1]

P <- matrix(0, ncol= length(idx_comp),nrow=length(idx_comp))
P <- P + diag(c( diag_P ) )
P[1,-1] <- z_nz
P[-1,1]<- nz_z
apply(P ,1,sum)
#normalization necessary due to removing some dist
col_s <- 1/ apply(P,1,sum)
P <- P*col_s
apply(P ,1,sum)

### try to write it it function form


iter_forward_backward <- function(X,sd,mu ,P,pi,thresh =1e-4 ){

  K <- ncol(P)
  if (missing(pi)){
    pi = rep( 1/length(mu), length(mu)) #same initial guess

  }
  emit = function(k,x,t){
    dnorm(x,mean=mu[k],sd=sd[t]  )
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

  print( idx_comp)
  ash_obj <- list()
  x_post <- 0*x
  for ( i in 1:length(idx_comp)){
    mu_ash <-  mu[idx_comp[i]]

    weight <- prob[,idx_comp[i]]

    ash_obj[[i]]  <- ash(x,sd,weight=weight,mode=mu_ash)
    x_post <-  x_post +weight*ash_obj[[i]]$result$PosteriorMean

  }
   prob  <- prob[ ,idx_comp]
   mu    <- mu  [  idx_comp]
   P     <- P   [idx_comp ,idx_comp]

   alpha <- alpha[ ,idx_comp]
   beta  <- beta [ ,idx_comp]

   out <- list(prob   = prob,
               fitted = prob%*%mu,
               x_post = x_post,
               mu     = mu,
               P      = P,
               alpha  = alpha,
               beta   = beta,
               pi     = prob[1, ]
               )



}





Baum_Welch <-  function(X,sd,mu,P, prob, alpha, beta ){

  list_z_nz <- list() # transition from 0 to non zero state
  list_nz_z <- list() # transition from  non zero state  to 0
  list_self <- list()# stay in the same state
  list_tt   <- list()

  emit = function(k,x,t){
    dnorm(x,mean=mu[k],sd=sd[t]  )
  }

  for ( t in 1:(length(X)-1)){
    tt_z_nz <-c()
    tt_nz_z <-c()
    tt_self <-c()
    for(  j in  2:ncol(P)){
      tt_z_nz <-c(tt_z_nz,
                  (alpha[t, 1]*P[1,j]*beta[t+1,j ]*emit(k=j,x=X[t+1], t= t+1 )  )  )
      # transition from 0 to non zero state
    }
    for(  j in  2:ncol(P)){
      tt_nz_z  <-c(tt_nz_z , (alpha[t, j]*P[1,j]*beta[  t+1,1 ]*emit(k=1,x=X[t+1], t= t+1 )  )  )
      # transition from  non zero state  to 0
    }
    for(  j in  1:ncol(P)){
      tt_self <-c(tt_self, (alpha[t, j]*P[j,j]*beta[  t+1,j ]*emit(k=j,x=X[t+1], t= t+1 )  )  )
      # stay in the same state
    }


    n_c <-  (sum(tt_z_nz) +sum(tt_nz_z )+ sum(tt_self)  )
    list_z_nz[[t]] <- tt_z_nz/ n_c  # transition from 0 to non zero state
    list_nz_z[[t]] <- tt_nz_z/  n_c   # transition from  non zero state  to 0
    list_self[[t]] <- tt_self/  n_c# stay in the same state

    list_tt[[t]] <-  c(tt_z_nz/ n_c,
                       tt_nz_z/  n_c ,
                       tt_self/  n_c )




  }
  expect_number_obs_state <- apply(prob[-nrow(  prob),    ],2,sum)

  #image (t(do.call(cbind,list_tt)))
  #plot(apply(prob[,-1],1, sum), type='l')

  transition <- apply( do.call(cbind,list_tt),1,function(j)  mean(j, na.rm=TRUE))# add log sum exp tric


  #Formula from Baum Welch update Wikipedi pagfe

  diag_P <- apply(do.call( cbind,list_self),1 ,sum) /expect_number_obs_state
  z_nz  <- apply(do.call( cbind,list_z_nz),1 ,sum) /expect_number_obs_state[1]
  nz_z  <- apply(do.call( cbind,list_nz_z ),1 ,sum) /rev( expect_number_obs_state[-1])
  nz_z  <- apply(do.call( cbind,list_nz_z ),1 ,sum) /  expect_number_obs_state[-1]

  P <- matrix(0, ncol= length(diag_P),nrow=length(diag_P))
  P <- P + diag(c( diag_P ) )
  P[1,-1] <- z_nz
  P[-1,1]<- nz_z
 # apply(P ,1,sum)
  #normalization necessary due to removing some dist
  col_s <- 1/ apply(P,1,sum)
  P <- P*col_s


  return( P)

}






for (o in 1:4){

  tt <-  iter_forward_backward (X  = res$Bhat[idx[2],],
                                sd = res$Shat[idx[2],],
                                mu = mu  ,
                                P=P,
                                pi = pi,
                                thresh =1e-4 )


  dim(tt$P)
  length(tt$mu)
  dim(tt$prob)

  pi <- tt$pi
  P <- Baum_Welch  (X    = res$Bhat[idx[2],],
                    sd   = res$Shat[idx[2],],
                    mu   = tt$mu,
                    P    =tt$P,
                    prob = tt$prob,
                    alpha= tt$alpha,
                    beta = tt$beta )
}
