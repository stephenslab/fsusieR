#' @title Simulate data under the mixture normal prior
#'
#' @description Add description here.
#'
#' @param lev_res numerical corresponds to the resolution of the simulated function (idealy between 3 and 10)
#'
#' @param length_grid vector numerical corresponds to the length of the grid of sigma for mixture component(cf ash)
#'
#' @param pi0 vector numerical , contain a digit  between 0 and 1, which corresponds to the null proportion ( non assocatied wavelet coefficients)
#'
#' @importFrom stats rchisq
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom wavethresh wd
#' @importFrom ashr normalmix
#'
#' @export
#'
#' @examples
#'
#'out <- simu_IBSS_ash_vanilla(lev_res=8, length_grid= 10, pi0= 0.85)
#'plot(out$sim_func, type="l", ylab="y")
#'out$emp_pi0
simu_IBSS_ash_vanilla <- function( lev_res=7, length_grid= 10, pi0= 0.85)
{
  #define mixutre components
  grid <- c(0,cumsum(rchisq(length_grid -1,df=1) )) # grid of sd for the ash mixture
  tt <- runif(length_grid-1)
  tt <- (1-pi0)*tt/sum(tt)
  pi_sim <- c(pi0,tt )# proportion of each of the mixture component, first element is prop of the null

  true_g = ashr::normalmix(pi_sim,rep(0, length_grid),grid) # define normalmix object used for ash

  #generating a set of wavelet coefficients under this model
  tem_func <- rep(0, 2^lev_res)
  twav <- wavethresh::wd(tem_func)
  while ( sum(twav$D ==0 ) == length(twav$D))#to ensure that there is at least one coefficient different from 0
  {
    clust  <-  c()
    for ( i in 1:length(twav$D))
    {
      clust <- c( clust, sample (1:length_grid, prob=pi_sim, size=1))
    }
    for( i in 1:length(twav$D))
    {
      twav$D[i] <-  ifelse(clust[i]==1,0,rnorm(1, mean=0, sd= grid[clust]))
    }
  }
  sim_func <- wavethresh::wr(twav)
  out <- list( sim_func  = sim_func,
               true_coef =twav$D,
               true_g=true_g,
               emp_pi0 =length(which(twav$D==0))/(2^lev_res))
  return(out)
}


#' @title Simulate data under the mixture normal prior
#'
#' @description Add description here.
#'
#' @param lev_res numerical corresponds to the resolution of the simulated function (idealy between 3 and 10)
#'
#' @param length_grid vector numerical corresponds to the length of the grid of sigma for mixture component(cf ash)
#'
#' @param pi0 vector numerical , contain a digit  between 0 and 1, which corresponds to the null proportion ( non assocatied wavelet coefficients)
#'
#' @param alpha numeric >0, control smoothness of the curves, should be positive and up 4 in particular d_sl ~  pi_{0,sl}  delta_0 + sum_k  pi_k N(0, 2^{- alpha * s}   sigma_k^2)
#'
#' @param prop_decay numeric >0, control the proportion of non zero wavelet coefficient per scale, pi_{0,sl} = 1- exp(-prop_decay*s)
#'
#' @importFrom stats rchisq
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom wavethresh wd
#' @importFrom wavethresh accessD
#' @importFrom ashr normalmix
#'
#'
#'
#' @export
#'
#' @examples
#'out <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay = 0.5)
#'plot(out$sim_func, type="l", ylab="y")
#'out$emp_pi0
#'temp_func <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay = 0)
#'print( temp_func$emp_pi0)
#'
simu_IBSS_per_level  <-function( lev_res=7,
                                 length_grid= 10,
                                 pi0,
                                 alpha=0.8,
                                 prop_decay=0.1)
{
  if(missing(pi0))
  {
    pi0 <- 1-exp(- (  prop_decay*(1:lev_res)))
  }
  #usefull to parse the indexes of wavelt coefficients

  indx_lst <- list()
  indx_lst[[1]] <- 1 #coefficient
  for ( s in 1:(lev_res-1))
  {

    indx  <- (2^s):(2^((s+1))-1)

    indx_lst[[s+1]] <- indx
  }
  #for the sake of ease same grid for each scale
  grid <- c(0,cumsum(rchisq(length_grid -1,df=3) ))# grid of sd for the ash mixture
  G_level <-list() #list of the mixture per scale
  for ( i in 1:lev_res)
  {

    tt <- runif(length_grid-1)
    tt <- (1-pi0[i])*tt/sum(tt)
    pi_sim <- c(pi0[i],tt )# proportion of each of the mixture component, first element is prop of the null


    G_level[[i]] <- ashr::normalmix(pi_sim,rep(0, length_grid),grid)# define normalmix object used for ash at lev res i
  }
  tem_func <- rep(0, 2^lev_res)
  twav <- wavethresh::wd(tem_func)

  while ( sum(twav$D ==0 ) == length(twav$D))#to ensure that ther is at least one coefficient different from 0
  {
    for ( i in 1:lev_res)
    {
      for ( k in unlist(indx_lst[i]))
      {

        clust <- sample (1:length_grid, prob=G_level[[i]]$pi, size=1)
        twav$D[k] <-   rnorm(1, mean=0, sd= grid[clust]/    ( 2^( alpha*i)))
        #print(twav$D[k])
      }

    }
    tt <- twav
    twav$D <- rev(twav$D)
  }

  #plot(accessD(twav,level=6), rev( tt$D[unlist(indx_lst[7])]) )
  sim_func <- wavethresh::wr(twav)
  emp_pi0 <- rep( 0, lev_res)
  for ( i in 0:(lev_res-1))
  {
    if(length(which(wavethresh::accessD(twav,level=i)==0)) ==0)
    {
      emp_pi0[i+1] <- 0
    }else{
      emp_pi0[i+1] <-  length(which(wavethresh::accessD(twav,level=i)==0))/(length(wavethresh::accessD(twav,level=i)))

    }
    # print( emp_pi0[i+1])

  }


  out <- list( sim_func  = sim_func,
               true_coef =twav$D,
               mix_per_scale=G_level,
               emp_pi0=emp_pi0)
  return(out)
}

# @title Simulation using Donoho and Johnstone test functions
#
# @param N integer number of sample to simulate
#
# @param P number of covariate
#
# @param lev_res control length of the generated function (length function = 2^lev_res)
#
# @param rsnr root signal noise ratio for the noise
#
# @param is.plot logical if set to TRUE plot underying function
#
# @param pos1 position of the first active covariate
#
# @param pos2 position of the first active covariate (optional)
#
# @importFrom wavethresh DJ.EX
# @importFrom graphics lines
# @importFrom graphics legend
# @importFrom graphics plot
#
# @export
#
simu_test_function <- function(N=50, P=10,lev_res=7, rsnr=2,is.plot=TRUE, pos1 =1, pos2)
{

  if(!missing(pos2)){
    if( pos1==pos2)
    {
      stop("Error pos2 and pos1 should be different byu default pos=1")
    }
  }

  G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
  beta0       <- 0
  beta1       <- 1
  beta2       <- ifelse(missing(pos2), 0,1)
  noisy.data  <- list()
  idx <- sample( size =3, 1:4)#sample at random the different function for basaline/effect
  if(!missing(pos2)){
    for ( i in 1:N)
    {
      test_func <- wavethresh::DJ.EX(n = 128, rsnr = rsnr, noisy = TRUE )
      f0        <- beta0*test_func[[idx[1]]] #Baseline
      f1        <- test_func[[idx[2]]]
      f2        <- test_func[[idx[3]]]
      noisy.data [[i]] <-  beta0*f0 +  beta1*G[i,pos1]*f1 + beta2*G[i,pos2]*f2

    }
  }else{
    for ( i in 1:N)
    {
      test_func <- wavethresh::DJ.EX(n = 2^lev_res, rsnr = rsnr, noisy = TRUE )
      f0        <- beta0*test_func[[idx[1]]] #Baseline
      f1        <- test_func[[idx[2]]]

      noisy.data [[i]] <-  beta0*f0 +  beta1*G[i,pos1]*f1

    }
  }

  noisy.data <- do.call(rbind, noisy.data)
  test_func <- wavethresh::DJ.EX(n = 2^lev_res,   noisy = FALSE )
  f0        <- beta0*test_func[[idx[1]]] #Baseline
  f1        <- test_func[[idx[2]]]
  f2        <- test_func[[idx[3]]]

  if( is.plot)
  {
    plot(f0-0.1,
         type="l",
         main="Underlying function depending on the SNP",
         ylim=c(3*min(f0,f1,f2),3*max(f0,f1,f2)),
         ylab="y",
         xlab="time"
    )
    lines(f1+0.1+f0, col="red")
    if( !missing(pos2)){
      lines(f2-0.2+f0, col="green")
      lines(f1+f2+f0+0.3, col="blue")
      legend(x = c(0),
             y= 100,
             c("0,0", "1,0", "0,1","1,1"),
             col= c("black", "red","green", "blue"),
             lty = rep(1,4)
      )
    }else{

      legend(x = c(0),
             y= 100,
             c("0 ", "1 " ),
             col= c("black", "red" ),
             lty = rep(1,2)
      )
    }





    if (!missing( pos2)){
      plot( noisy.data[1,],
            col= "black",
            type ="l",
            ylim=c(-150,150),
            ylab="y",
            xlab="time",
            main="Observed noisy curves"
      )
      for ( i  in 2: N)
      {
        if( G[i, pos1]==0  & G[i,pos2]==0)
        {
          my_col <- "black"
        }

        if( G[i, pos1]>0  & G[i,pos2]==0)
        {
          my_col <- "red"
        }
        if( G[i, pos1]==0  & G[i,pos2]>0)
        {
          my_col <- "green"
        }
        if( G[i, pos1]>0  & G[i,pos2]>0)
        {
          my_col <- "blue"
        }

        lines( noisy.data[i,], col=   my_col)
        legend(x = c(0),
               y= -5,
               c("0,0", "1,0", "0,1","1,1"),
               col= c("black", "red","green", "blue"),
               lty = rep(1,4)
        )
      }

    }else{
      plot( noisy.data[1,],
            col= "black",
            type ="l",
            ylim=c(-150,150),
            ylab="y",
            xlab="time",
            main="Observed noisy curves"
      )
      for ( i  in 2: N)
      {
        if( G[i, pos1]==0  )
        {
          my_col <- "black"
        }

        if( G[i, pos1]>0  )
        {
          my_col <- "red"
        }


        lines( noisy.data[i,], col=   my_col)


      }
      legend(x = c(0),
             y= -5,
             c("0", "1" ),
             col= c("black", "red"),
             lty = rep(1,2)
      )
    }
  }
  if(!missing(pos2))
  {
    out <- list( G=G,
                 noisy.data = noisy.data,
                 pos1       = pos1,
                 pos2       = pos2,
                 f0         = f0,
                 f1         = f1,
                 f2         = f2)
  }else{
    out <- list( G=G,
                 noisy.data = noisy.data,
                 pos1       = pos1,
                 f0         = f0,
                 f1         = f1 )
  }
  return(out)
}




