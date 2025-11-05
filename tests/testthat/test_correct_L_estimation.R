#test if the null weight is not high enough then fails
library(susieR)
library(fsusieR)
set.seed(1)
N <- 100 # Sample size
P <- 100 # Number of SNPs
L <- 2 # Number of effects 
  
# Genotype matrix
G <- N3finemapping$X[sample(1:nrow(N3finemapping$X),size=N),]
if ( length(which(apply(G,2,var)<0.0001))>0){
  G=G[,-which(apply(G,2,var)<0.0001)]
}
list_lev_res <- list(5, 6)
# Two functional phenotypes, one of length 2^5 and one of length 2^6
n_univ <- 3 # 3 univariate phenotypes
eff <- list()
for (l in 1:L) { # Simulate the multi-trait effect
  eff[[l]] <- simu_effect_multfsusie(list_lev_res=list_lev_res,
                                     n_univ=n_univ, output_level=2)
}

Y_f1 <- matrix(rnorm((length(c(eff[[l]]$func_effect[[1]]$sim_func,eff[[l]]$func_effect[[2]]$sim_func[10:30]))) * N, sd=1), nrow=N)
Y_f2 <- matrix(rnorm((length(c(eff[[l]]$func_effect[[2]]$sim_func,eff[[l]]$func_effect[[1]]$sim_func[10:30]))) * N, sd=1), nrow=N)
Y_u <- matrix(rnorm(n_univ * N, sd=1), nrow=N)

true_pos <- sample(1:ncol(G), L) # Actually causal columns/SNPs

for (i in 1:N) {
  for (l in 1:L) {
    Y_f1[i,] <- Y_f1[i,] + c(eff[[l]]$func_effect[[1]]$sim_func,eff[[l]]$func_effect[[2]]$sim_func[10:30]) * G[i, true_pos[[l]]]
    Y_f2[i,] <- Y_f2[i,] + c(eff[[l]]$func_effect[[2]]$sim_func,eff[[l]]$func_effect[[1]]$sim_func[10:30])  * G[i, true_pos[[l]]]
    Y_u[i,]  <- Y_u[i,]  + eff[[l]]$univ_effect * G[i, true_pos[[l]]]
  }
}

Y_f <- list(Y_f1, Y_f2)
Y_f_1 <- list(Y_f1 )
Y_f_2 <- list(Y_f2 )
Y <- list(Y_f=Y_f, Y_u=Y_u) # Preparing data
Y1 <- list(Y_f=Y_f_1, Y_u=NULL) # Preparing data
Y2 <- list(Y_f=Y_f_2, Y_u=NULL) # Preparing data

pos <- list(pos1=1:ncol(Y$Y_f[[1]]), pos2=1:ncol(Y$Y_f[[2]])

  m1 <-   susiF(Y=Y$Y_f[[1]], X=G,L=3, verbose=FALSE,
                        control_mixsqp =  list(verbose=FALSE,
                                               eps = 1e-6,
                                               numiter.em = 40
                        ),
                        nullweight= 10
  )
  
  
  m2 <-  susiF(Y=Y$Y_f[[2]], X=G,L=3, verbose=FALSE,
                        control_mixsqp =  list(verbose=FALSE,
                                               eps = 1e-6,
                                               numiter.em = 40
                        ),
                        nullweight= 10 )
  
  
 
 
  n_false_fsusie1 = Reduce("+",sapply(1:length(m1$cs), function(k)
    ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
  ))
  n_false_fsusie2 = Reduce("+",sapply(1:length(m2$cs), function(k)
    ifelse( length(which(true_pos%in%m2$cs[[k]] ))==0, 1,0)
  ))
  
  test_that("susiF correctly estimate number of effects", {
    expect_equal(n_false_fsusie1,0)
    expect_lte(n_false_fsusie2, 0)
  })

  
  
  
  N <- 100 # Sample size
  P <- 100 # Number of SNPs
  L <- 2 # Number of effects 
  
  # Genotype matrix
  G <- N3finemapping$X[sample(1:nrow(N3finemapping$X),size=N),]
  if ( length(which(apply(G,2,var)<0.0001))>0){
    G=G[,-which(apply(G,2,var)<0.0001)]
  }
  list_lev_res <- list(5, 6)
  # Two functional phenotypes, one of length 2^5 and one of length 2^6
  n_univ <- 3 # 3 univariate phenotypes
  eff <- list()
  for (l in 1:L) { # Simulate the multi-trait effect
    eff[[l]] <- simu_effect_multfsusie(list_lev_res=list_lev_res,
                                       n_univ=n_univ, output_level=2)
  }
  
  Y_f1 <- matrix(rnorm((length(c(eff[[l]]$func_effect[[1]]$sim_func,eff[[l]]$func_effect[[2]]$sim_func[10:30]))) * N, sd=1), nrow=N)
  Y_f2 <- matrix(rnorm((length(c(eff[[l]]$func_effect[[2]]$sim_func,eff[[l]]$func_effect[[1]]$sim_func[10:30]))) * N, sd=1), nrow=N)
  Y_u <- matrix(rnorm(n_univ * N, sd=1), nrow=N)
  
  true_pos <- sample(1:ncol(G), L) # Actually causal columns/SNPs
  
  for (i in 1:N) {
    for (l in 1:L) {
      Y_f1[i,] <- Y_f1[i,] + c(eff[[l]]$func_effect[[1]]$sim_func,eff[[l]]$func_effect[[2]]$sim_func[10:30]) * G[i, true_pos[[l]]]
      Y_f2[i,] <- Y_f2[i,] + c(eff[[l]]$func_effect[[2]]$sim_func,eff[[l]]$func_effect[[1]]$sim_func[10:30])  * G[i, true_pos[[l]]]
      Y_u[i,]  <- Y_u[i,]  + eff[[l]]$univ_effect * G[i, true_pos[[l]]]
    }
  }
  
  Y_f <- list(Y_f1, Y_f2)
  Y_f_1 <- list(Y_f1 )
  Y_f_2 <- list(Y_f2 )
  Y <- list(Y_f=Y_f, Y_u=Y_u) # Preparing data
  Y1 <- list(Y_f=Y_f_1, Y_u=NULL) # Preparing data
  Y2 <- list(Y_f=Y_f_2, Y_u=NULL) # Preparing data
  
  pos <- list(pos1=1:ncol(Y$Y_f[[1]]), pos2=1:ncol(Y$Y_f[[2]])
              
              m1 <-   susiF(Y=Y$Y_f[[1]], X=G,L=3, verbose=FALSE,
                            control_mixsqp =  list(verbose=FALSE,
                                                   eps = 1e-6,
                                                   numiter.em = 40
                            ),
                            nullweight= 10
              )
              
              
              m2 <-  susiF(Y=Y$Y_f[[2]], X=G,L=3, verbose=FALSE,
                           control_mixsqp =  list(verbose=FALSE,
                                                  eps = 1e-6,
                                                  numiter.em = 40
                           ),
                           nullweight= 10 )
              
              
              
              
              n_false_fsusie1 = Reduce("+",sapply(1:length(m1$cs), function(k)
                ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
              ))
              n_false_fsusie2 = Reduce("+",sapply(1:length(m2$cs), function(k)
                ifelse( length(which(true_pos%in%m2$cs[[k]] ))==0, 1,0)
              ))
              
              test_that("susiF correctly estimate number of effects", {
                expect_equal(n_false_fsusie1,0)
                expect_lte(n_false_fsusie2, 0)
              })
              