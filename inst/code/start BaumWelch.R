halfK=50
K = 2*halfK-1
sd=  sd
X <-x

pos <- seq(  0, 1 ,   length.out=halfK)

#define the mean states
mu <- (pos^(1/mult))*1.5*max(abs(X)) # put 0 state at
#the firstplace
mu <- c(mu, -mu[-1] )



P <- diag(0.999,K) #this ensure that the HMM can only "transit via null state"
P[1,-1] <- 0.001
P[-1,1] <- 0.001

if( smooth){
  P <- diag(0.95,K) #this ensure that the HMM can only "transit via null state"
  P[1,-1] <- 0.05
  P[-1,1] <- 0.05
  for ( j in ceiling(halfK/2):(halfK+1)){

    P[1,j] <-  0.025
    P[j,1]  <- 0.025
    P[j,(j+1)] <- 0.025
    P[(j+1),j] <-0.025

  }

  for ( j in (halfK+1 +ceiling(halfK/2)):( K -1)){

    P[1,j] <-  0.025
    P[j,1]  <- 0.025
    P[j,(j+1)] <- 0.025
    P[(j+1),j] <- 0.025

  }
}



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


list_tt <- list()
for ( t in 1:(length(X)-1)){
  tt <-c()

  for(  j in idx_comp){
    tt <-c(tt, (alpha[t, 1]*P[1,j]*beta[t+1,j ]*emit(k=j,x=X[t+1], t= t+1 )  )  )
    # transition from 0 to non zero state
  }
  norm_const <-  Reduce("+", lapply( 2:K, function( j)
    (alpha[t, j]*P[1,j]*beta[t+1,1 ]*emit(k=1,x=X[t+1], t= t+1 )  )
  )
  )
  print( c(norm_const,sum(tt)+norm_const ) )
  list_tt[[t]] <- tt/ (sum(tt)+norm_const  )




}
image (do.call(cbind,list_tt))
transition <- apply( do.call(cbind,list_tt),1,function(j)  mean(j, na.rm=TRUE))# add log sum exp tric





P <- matrix(0, ncol= length(idx_comp),nrow=length(idx_comp))
P <- P + diag(c(   1- transition ) )
P[1,-1] <- transition[-1]
P[-1,1]<- transition[-1]



alpha = matrix(nrow = length(X),ncol=length(idx_comp))
alpha[1, ] =  prob[1,idx_comp]

