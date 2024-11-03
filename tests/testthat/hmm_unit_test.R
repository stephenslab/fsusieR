##unit test HMM



# sharp transition and low variance
library(testthat)

devtools::load_all(".")

library(fsusieR)
library(ashr)
library(wavethresh)
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
sd_noise <- 0.1 #expected root signal noise ratio
N <- 100    #Number of individuals
P <- 10     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 5   #Position of the causal covariate for effect 2
lev_res <-7#length of the molecular phenotype (2^lev_res)
f1 <-  rep(0, 2^lev_res)
f1[ 20:25] <-2
f1[ 50:55] <-1
f1[ 20:25] <-2

f1 <-  0*simu_IBSS_per_level(lev_res )$sim_func
f1[60:length(f1)] <-0
f1[ 70] <- -1
#first effect)
f2 <-  0.1*DJ.EX(128)$blocks

beta0       <- 0
beta1       <- 1
beta2       <- 1
noisy.data  <- list()
#f2 <-f1
G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
X <-G
for ( i in 1:N)
{
  noise <- rnorm(length(f1), sd=  sd_noise)
  noisy.data [[i]] <-   beta1*G[i,pos1]*f1 + beta1*G[i,pos2]*f2  + noise

}
noisy.data <- do.call(rbind, noisy.data)
Y <- noisy.data

out <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale', filter.number =8  )


X <- colScale(X)
# centering input
Y <- colScale(Y, scale=FALSE)
susiF.obj <- out


idx <- do.call( c, lapply( 1:length(susiF.obj$cs),
                           function(l){
                             tp_id <-  which.max( susiF.obj$pip[susiF.obj$cs[[l]]])
                             susiF.obj$cs[[l]][tp_id]
                           }
)
)

temp_Y <- Y
res <- cal_Bhat_Shat(temp_Y,X )






fitted_trend <- list()
est_prob <- list()
for ( k in 1:2){
  for (j in 1:length(idx)){
    res <- cal_Bhat_Shat(temp_Y,X )

    s =fit_hmm(x=res$Bhat[idx[j],],sd=res$Shat[idx[j],],halfK=50 )

    est_prob[[j]] <-  s$prob[,1]

    fitted_trend[[j]] <- s$x_post
    if( j ==length(idx)){
      idx_var <- (1:length(idx)) [- (1)]
    }else{
      idx_var <- (1:length(idx))[- (j+1)]
    }


    temp_Y <- Y - Reduce("+", lapply( idx_var, function( j){
      X[,idx[j] ]%*%t(fitted_trend[[j]])
    }

    )
    )

  }

}

fitted_trend <- lapply(1:length(idx), function(l)
  fitted_trend[[l]]/susiF.obj$csd_X[idx[l]]
)


test_that("performance in low variance and sharp transition should be",{
  expect_gt( cor(f1, fitted_trend[[2]] ), 0.9999992  )
  expect_gt( cor(f2, fitted_trend[[1]] ), 0.9999  )
  expect_equal(which(est_prob[[2]]<0.05),70)
  expect_equal(which(est_prob[[1]]<0.05),13:103)

}

)





### test when constant obs
Y[,10:20] <- 0*Y[,10:20]# low variance




out <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale', filter.number =8  )


X <- colScale(X)
# centering input
Y <- colScale(Y, scale=FALSE)
susiF.obj <- out


idx <- do.call( c, lapply( 1:length(susiF.obj$cs),
                           function(l){
                             tp_id <-  which.max( susiF.obj$pip[susiF.obj$cs[[l]]])
                             susiF.obj$cs[[l]][tp_id]
                           }
)
)

temp_Y <- Y
res <- cal_Bhat_Shat(temp_Y,X )







est_prob <- list()
for ( k in 1:2){
  for (j in 1:length(idx)){
    res <- cal_Bhat_Shat(temp_Y,X )

    s =fit_hmm(x=res$Bhat[idx[j],],sd=res$Shat[idx[j],],halfK=50 )

    est_prob[[j]] <-  s$prob[,1]

    fitted_trend[[j]] <- s$x_post
    if( j ==length(idx)){
      idx_var <- (1:length(idx)) [- (1)]
    }else{
      idx_var <- (1:length(idx))[- (j+1)]
    }


    temp_Y <- Y - Reduce("+", lapply( idx_var, function( j){
      X[,idx[j] ]%*%t(fitted_trend[[j]])
    }

    )
    )

  }

}

fitted_trend <- lapply(1:length(idx), function(l)
  fitted_trend[[l]]/susiF.obj$csd_X[idx[l]]
)





f2_tilde <- f2
f2_tilde[10:20]<-0*f2_tilde[10:20]
test_that("performance in low variance and sharp transition should be",{
  expect_gt( cor(f1, fitted_trend[[2]] ), 0.9999992  )
  expect_gt( cor(f2_tilde, fitted_trend[[1]] ), 0.9999  )
  expect_equal(which(est_prob[[2]]<0.05),70)
  expect_equal(which(est_prob[[1]]<0.05),21:103)

}

)

