
#lev_res numerical corresponds to the resolution of the simulated function (idealy between 3 and 10)
#length_grid vector numerical corresponds to the length of the grid of sigma for mixture component(cf ash)
#piO vector numerical , contain a digit  between 0 and 1, which corresponds to the null proportion ( non assocatied wavelet coefficients)
simu_IBSS_ash_vanilla <- function( lev_res=7, length_grid= 10, pi0= 0.85)
{
  #define mixutre components
  grid <- c(0,cumsum(rchisq(length_grid -1,df=1) )) # grid of sd for the ash mixture
  tt <- runif(length_grid-1)
  tt <- (1-pi0)*tt/sum(tt)
  pi_sim <- c(pi0,tt )# proportion of each of the mixture component, first element is prop of the null

  true_g = normalmix(pi_sim,rep(0, length_grid),grid) # define normalmix object used for ash

  #generating a set of wavelet coefficients under this model
  tem_func <- rep(0, 2^lev_res)
  twav <- wd(tem_func)
  while ( sum(twav$D ==0 ) == length(twav$D))#to ensure that ther is at least one coefficient different from 0
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
  sim_func <- wr(twav)
  out <- list( sim_func  = sim_func,
               true_coef =twav$D,
               true_g=true_g,
               emp_pi0 =length(which(twav$D==0))/(2^lev_res))
  return(out)
}

not_run <- FALSE
if(not_run)
{

  out <- simu_IBSS_ash_vanilla(lev_res=8, length_grid= 10, pi0= 0.85)
  plot(out$sim_func, type="l", ylab="y")
  out$emp_pi0
}
#lev_res numerical corresponds to the resolution of the simulated function (idealy between 3 and 10)
#piO vector of length lev_res, containing digits between 0 and 1, which correspond to the mixture component for each level of resolution to be zero.
# by default set as 1-1/2^(0.1*(1:lev_res))
#length_grid vector numerical corresponds the length of the grid of sigma for mixture component(cf ash)
#alpha numeric >0, control smoothness of the curves, should be positive and up 4 in particular d_sl ~ \pi_{0,sl} \delta_0 + \sum_k \pi_k N(0, 2^{-\alpha * s}  \sigma_k^2)
#prop_decay numeric >0, control the proportion of non zero wavelet coefficient per scale, \pi_{0,sl} = 1- exp(-prop_decay*s)

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


    G_level[[i]] <- normalmix(pi_sim,rep(0, length_grid),grid)# define normalmix object used for ash at lev res i
  }
  tem_func <- rep(0, 2^lev_res)
  twav <- wd(tem_func)

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
  sim_func <- wr(twav)
  emp_pi0 <- rep( 0, lev_res)
  for ( i in 0:(lev_res-1))
  {
    if(length(which(accessD(twav,level=i)==0)) ==0)
    {
      emp_pi0[i+1] <- 0
    }else{
      emp_pi0[i+1] <-  length(which(accessD(twav,level=i)==0))/(length(accessD(twav,level=i)))

    }
    # print( emp_pi0[i+1])

  }


  out <- list( sim_func  = sim_func,
               true_coef =twav$D,
               mix_per_scale=G_level,
               emp_pi0=emp_pi0)
  return(out)
}

not_run <- FALSE
if(not_run)
{
  out <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay = 0.5)
  plot(out$sim_func, type="l", ylab="y")
  out$emp_pi0

  temp_func <- simu_IBSS_per_level(lev_res=9, alpha=1, prop_decay = 0.)
  print( temp_func$emp_pi0)



}
