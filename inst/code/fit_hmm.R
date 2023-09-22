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
                     thresh_prefilter=1e-3,
                     maxiter=10){


  K = 2*halfK-1
  sd=  sd
  X <-x

  pos <- seq(  0, 1 ,   length.out=halfK)

  #define the mean states
  mu <- (pos^(1/mult))*1.5*max(abs(X)) # put 0 state at
  #the firstplace
  mu <- c(mu, -mu[-1] )


  min_delta <- abs(mu[2]-mu[1])
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
    alpha_hat[1, ] = pi* emit(1:K,x=X[1],t=1)
    alpha_tilde[1, ] = pi* emit(1:K,x=X[1],t=1)
  }




  # Forward algorithm -----
  for(t in 1:(length(X)-1)){
    m = alpha_hat[t,] %*% P

    alpha_tilde[t+1, ] = m *emit(1:K,x=X[t+1], t= t+1 )

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

  # Backwards algorithm-----
  for(t in ( length(X)-1):1){



    emissio_p <- emit(1:K,X[t+1],t=t+1)





    beta_tilde [t, ] = apply( sweep( P,2, beta_hat[t+1,]*emissio_p ,"*" ),1,sum)


    C_t[t] <- max(beta_tilde[t,])
    beta_hat[t,] <-  beta_tilde [t, ] /C_t[t]
  }




  ab = alpha_hat*beta_hat
  prob = ab/rowSums(ab)


  #image(prob)#plot(apply(prob[,-1],1, sum), type='l')
  #plot(x)
  #lines(1-prob[,1])







  #Baum_Welch-----
  #Baum_Welch <-  function(X,sd,mu,P, prob, alpha, beta ){

  list_z_nz <- list() # transition from 0 to non zero state
  list_nz_z <- list() # transition from  non zero state  to 0
  list_self <- list()# stay in the same state



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
    if( n_c==0){
      list_z_nz[[t]] <- tt_z_nz*0  # transition from 0 to non zero state
      list_nz_z[[t]] <- tt_nz_z*0    # transition from  non zero state  to 0
      list_self[[t]] <- tt_self*0 # stay in the same state
    }else{
      list_z_nz[[t]] <- tt_z_nz/ n_c  # transition from 0 to non zero state
      list_nz_z[[t]] <- tt_nz_z/  n_c   # transition from  non zero state  to 0
      list_self[[t]] <- tt_self/  n_c# stay in the same state

    }




  }
  expect_number_obs_state <- apply(prob[-nrow(  prob),    ],2,sum)

  #image (t(do.call(cbind,list_tt)))
  #plot(apply(prob[,-1],1, sum), type='l')



  #Formula from Baum Welch update Wikipedi pagfe

  diag_P <- apply(do.call( cbind,list_self),1 ,sum) /expect_number_obs_state
  z_nz  <- apply(do.call( cbind,list_z_nz),1 ,sum) /expect_number_obs_state[1]
  nz_z  <- apply(do.call( cbind,list_nz_z ),1 ,sum) /  expect_number_obs_state[-1]

  P <- matrix(0, ncol= length(diag_P),nrow=length(diag_P))
  P <- P + diag(c( diag_P ) )
  P[1,-1] <- z_nz
  P[-1,1]<- nz_z
  # apply(P ,1,sum)
  #normalization necessary due to removing some dist
  col_s <- 1/ apply(P,1,sum)
  P <- P*col_s



  idx_comp <- which( apply(prob, 2, mean) >thresh )
  if ( !(1%in% idx_comp) ){ #ensure 0 is in the model

    idx_comp<- c(1, idx_comp)
  }

  ash_obj <- list()
  x_post <- 0*x


  for ( i in 2:length(idx_comp)){
    mu_ash <-mu[idx_comp[i] ]
    weight <- prob[,idx_comp[i]]

    ash_obj[[i]]  <- ash(x,sd,
                         weight=weight,
                         mode=mu_ash,
                         mixcompdist = "normal"
    )
    x_post <-  x_post +weight*ash_obj[[i]]$result$PosteriorMean

  }
  P <-P[idx_comp, idx_comp]
  K <- length( idx_comp)
  mu <- mu[idx_comp]





  iter =1
  plot( X)

  while( iter <maxiter){


    alpha_hat = matrix(nrow = length(X),ncol=K)
    alpha_tilde = matrix(nrow = length(X),ncol=K)
    G_t <- rep(NA, length(X))

    data0 <-  set_data(X[1],sd[1])

    pi <- (rep(1/K,K))

    alpha_hat[1, ] = pi  *c(dnorm(X[1], mean=0, sd=sd[1]),
                            sapply( 2:K, function( k) exp(calc_loglik(ash_obj[[k]],
                                                                      data0)
                            )
                            )
    )
    alpha_tilde[1, ] = pi  *c(dnorm(X[1], mean=0, sd=sd[1]),
                              sapply( 2:K, function( k) exp(calc_loglik(ash_obj[[k]],
                                                                        data0)
                              )
                              )
    )





    # Forward algorithm
    for(t in 1:(length(X)-1)){
      m = alpha_hat[t,] %*% P
      data0 <-  set_data(X[t],sd[t])



      alpha_tilde[t+1, ] = m  *c(dnorm(X[t], mean=0, sd=sd[t]),
                                 sapply( 2:K, function( k) exp(calc_loglik(ash_obj[[k]],
                                                                           data0)
                                 )
                                 )
      )
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

      data0 <-  set_data(X[t+1],sd[t+1])
      emissio_p <- c(dnorm(X[t+1], mean=0, sd=sd[t+1]),
                     sapply( 2:K, function( k) exp(calc_loglik(ash_obj[[k]],
                                                               data0)
                     )
                     )
      )



      beta_tilde [t, ] = apply( sweep( P,2, beta_hat[t+1,]*emissio_p ,"*" ),1,sum)

      C_t[t] <- max(beta_tilde[t,])
      beta_hat[t,] <-  beta_tilde [t, ] /C_t[t]
    }




    ab = alpha_hat*beta_hat
    prob = ab/rowSums(ab)


    # image(prob)#plot(apply(prob[,-1],1, sum), type='l')
    #plot(x)
    #lines(1-prob[,1])

    ash_obj <- list()
    x_post <- 0*x

    for ( k in 2:K){
      # mu_ash <- sum(prob[,k]*X)/(sum(prob[,k])) #M step for the mean
      mu_ash <-mu[k ]



      weight <- prob[,k]

      ash_obj[[k]]  <- ash(x,sd,
                           weight=weight,
                           mode=mu_ash,
                           mixcompdist = "normal"
      )
      x_post <-  x_post +weight*ash_obj[[k]]$result$PosteriorMean

    }

    #Baum_Welch-----
    #Baum_Welch <-  function(X,sd,mu,P, prob, alpha, beta ){

    list_z_nz <- list() # transition from 0 to non zero state
    list_nz_z <- list() # transition from  non zero state  to 0
    list_self <- list()# stay in the same state



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
      if( n_c==0){
        list_z_nz[[t]] <- tt_z_nz*0  # transition from 0 to non zero state
        list_nz_z[[t]] <- tt_nz_z*0    # transition from  non zero state  to 0
        list_self[[t]] <- tt_self*0 # stay in the same state
      }else{
        list_z_nz[[t]] <- tt_z_nz/ n_c  # transition from 0 to non zero state
        list_nz_z[[t]] <- tt_nz_z/  n_c   # transition from  non zero state  to 0
        list_self[[t]] <- tt_self/  n_c# stay in the same state

      }





    }
    expect_number_obs_state <- apply(prob[-nrow(  prob),    ],2,sum)

    #image (t(do.call(cbind,list_tt)))
    #plot(apply(prob[,-1],1, sum), type='l')



    #Formula from Baum Welch update Wikipedi pagfe

    diag_P <- apply(do.call( cbind,list_self),1 ,sum) /expect_number_obs_state
    z_nz  <- apply(do.call( cbind,list_z_nz),1 ,sum) /expect_number_obs_state[1]
    nz_z  <- apply(do.call( cbind,list_nz_z ),1 ,sum) /  expect_number_obs_state[-1]

    P <- matrix(0, ncol= length(diag_P),nrow=length(diag_P))
    P <- P + diag(c( diag_P ) )
    P[1,-1] <- z_nz
    P[-1,1]<- nz_z
    # apply(P ,1,sum)
    #normalization necessary due to removing some dist
    col_s <- 1/ apply(P,1,sum)
    P <- P*col_s
    P[is.na(P)] <- 0
    iter =iter +1
    lines( x_post, col=iter)


    print( sum(log(G_t[-1])))
  }
  out <- list( prob =prob,
               fitted= prob%*%mu,
               x_post = x_post,
               mu= mu)



}

