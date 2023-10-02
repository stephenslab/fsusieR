#x=res$Bhat[idx[2],]
#sd=res$Shat[idx[2],]
#halfK=20
#mult=2
#thresh=0.00001


fit_hmm <- function (x,sd,
                     halfK=20,
                     mult=2,
                     smooth=FALSE,
                     thresh=0.00001,
                     prefilter=TRUE,
                     thresh_prefilter=1e-3){


K = 2*halfK-1
sd=  sd
X <-x

pos <- seq(  0, 1 ,   length.out=halfK)

#define the mean states
mu <- (pos^(1/mult))*1.5*max(abs(X)) # put 0 state at
#the firstplace
mu <- c(mu, -mu[-1] )


if( prefilter){
  tt <- apply(
    do.call( rbind, lapply(1:length(x), function( i){
      tt <-dnorm(x[i],mean = mu, sd=sd[i])
      return( tt/ sum(tt))
    } )),

    2,
    mean)
  temp_idx <- which(tt > thresh_prefilter)
  if( 1 %!in% temp_idx){
    temp_idx <- c(1,temp_idx)
  }
  mu <- mu[temp_idx]
  K <- length(mu)
}

P <- diag(0.999,K) #this ensure that the HMM can only "transit via null state"
P[1,-1] <- 0.001
P[-1,1] <- 0.001



pi = rep( 1/length(mu), length(mu)) #same initial guess

emit = function(k,x,t){
  dnorm(x,mean=mu[k],sd=sd[t]  )
}


alpha_hat = matrix(nrow = length(X),ncol=K)
alpha_tilde = matrix(nrow = length(X),ncol=K)
G_t <- rep(NA, length(X))
for(k in 1:K){
  alpha_hat[1,k] = pi[k] * emit(k=k,x=X[1],t=1)
  alpha_tilde[1,k] = pi[k] * emit(k=k,x=X[1],t=1)
}





# Forward algorithm
for(t in 1:(length(X)-1)){
  m = alpha_hat[t,] %*% P
  for(k in 1:K){
    alpha_tilde[t+1,k] = m[k]*emit(k=k,x=X[t+1], t= t+1 )
  }
  G_t[t+1] <- sum( alpha_tilde[t+1,])
  alpha_hat[t+1,] <-  alpha_tilde[t+1,]/ ( G_t[t+1])
}

beta_hat = matrix(nrow =  length(X),ncol=K)

beta_tilde = matrix(nrow =  length(X),ncol=K)
C_t <- rep(NA, length(X))
# Initialize beta
for(k in 1:K){
  beta_hat[ length(X),k] = 1
  beta_tilde [ length(X),k] = 1
}

# Backwards algorithm
for(t in ( length(X)-1):1){
  for(k in 1:K){
    beta_tilde [t,k] = sum(beta_hat[t+1,]*P[k,]*emit(1:K,X[t+1],t=t+1))
  }
  C_t[t] <- max(beta_tilde[t,])
  beta_hat[t,] <-  beta_tilde [t, ] /C_t[t]
}




ab = alpha_hat*beta_hat
prob = ab/rowSums(ab)


#image(prob)#plot(apply(prob[,-1],1, sum), type='l')
#plot(x)
#lines(1-prob[,1])








#Baum_Welch <-  function(X,sd,mu,P, prob, alpha, beta ){

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
                  (alpha_hat[t, 1]*P[1,j]*beta_hat[t+1,j ]*emit(k=j,x=X[t+1], t= t+1 )  )  )
      # transition from 0 to non zero state
    }
    for(  j in  2:ncol(P)){
      tt_nz_z  <-c(tt_nz_z , (alpha_hat[t, j]*P[1,j]*beta_hat[  t+1,1 ]*emit(k=1,x=X[t+1], t= t+1 )  )  )
      # transition from  non zero state  to 0
    }
    for(  j in  1:ncol(P)){
      tt_self <-c(tt_self, (alpha_hat[t, j]*P[j,j]*beta_hat[  t+1,j ]*emit(k=j,x=X[t+1], t= t+1 )  )  )
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


  #image(P)#return( P)

#}





  emit = function(k,x,t){
    dnorm(x,mean=mu[k],sd=sd[t]  )
  }


  alpha_hat = matrix(nrow = length(X),ncol=K)
  alpha_tilde = matrix(nrow = length(X),ncol=K)
  G_t <- rep(NA, length(X))
  for(k in 1:K){
    alpha_hat[1,k] = pi[k] * emit(k=k,x=X[1],t=1)
    alpha_tilde[1,k] = pi[k] * emit(k=k,x=X[1],t=1)
  }





  # Forward algorithm
  for(t in 1:(length(X)-1)){
    m = alpha_hat[t,] %*% P
    for(k in 1:K){
      alpha_tilde[t+1,k] = m[k]*emit(k=k,x=X[t+1], t= t+1 )
    }
    G_t[t+1] <- sum( alpha_tilde[t+1,])
    alpha_hat[t+1,] <-  alpha_tilde[t+1,]/ ( G_t[t+1])
  }

  beta_hat = matrix(nrow =  length(X),ncol=K)

  beta_tilde = matrix(nrow =  length(X),ncol=K)
  C_t <- rep(NA, length(X))
  # Initialize beta
  for(k in 1:K){
    beta_hat[ length(X),k] = 1
    beta_tilde [ length(X),k] = 1
  }

  # Backwards algorithm
  for(t in ( length(X)-1):1){
    for(k in 1:K){
      beta_tilde [t,k] = sum(beta_hat[t+1,]*P[k,]*emit(1:K,X[t+1],t=t+1))
    }
    C_t[t] <- max(beta_tilde[t,])
    beta_hat[t,] <-  beta_tilde [t, ] /C_t[t]
  }




  ab = alpha_hat*beta_hat
  prob = ab/rowSums(ab)


  #image(prob)#plot(apply(prob[,-1],1, sum), type='l')
  #plot(x)
  #lines(1-prob[,1])




  idx_comp <- which( apply(prob, 2, mean) >thresh )
  if ( !(1%in% idx_comp) ){ #ensure 0 is in the model

    idx_comp<- c(1, idx_comp)
  }
  ash_obj <- list()
  x_post <- 0*x



  for ( i in 2:length(idx_comp)){
    mu_ash <- sum(prob[,idx_comp[i]]*X)/(sum(prob[,idx_comp[i]])) #M step for the mean

    weight <- prob[,idx_comp[i]]

    ash_obj[[i]]  <- ash(x,sd,
                         weight=weight,
                         mode=mu_ash,
                         mixcompdist = "normal"
                         )
    x_post <-  x_post +weight*ash_obj[[i]]$result$PosteriorMean

  }




  out <- list( prob =prob,
               fitted= prob%*%mu,
               x_post = x_post,
               mu= mu)



}

