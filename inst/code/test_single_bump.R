rm(list=ls())
library(susiF.alpha)
library(ashr)
library(wavethresh)
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
sd_noise <- 1 #expected root signal noise ratio
N <- 100    #Number of individuals
P <- 10     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 5   #Position of the causal covariate for effect 2
lev_res <-7#length of the molecular phenotype (2^lev_res)
f1 <-  rep(0, 2^lev_res)
f1[ 20:25] <-2
f1[ 50:55] <-1
f1[ 20:25] <-2
f1[ 60:75] <- -1
#first effect)
f2 <-  rep(0, 2^lev_res)
f2[ 10:25] <-1
f2[ 50:65] <-1
f1[ 60:75] <- -1
plot( f1, type ="l", ylab="effect", col="blue")
lines(f2, col="red")
beta0       <- 0
beta1       <- 1
beta2       <- 1
noisy.data  <- list()
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

out2 <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale', TI=FALSE,filter.number =8  )


X <- colScale(X)
# centering input
Y <- colScale(Y, scale=FALSE)
susiF.obj <- out
L_points=20

idx <- do.call( c, lapply( 1:length(susiF.obj$cs),
                           function(l){
                             tp_id <-  which.max( susiF.obj$pip[susiF.obj$cs[[l]]])
                             susiF.obj$cs[[l]][tp_id]
                           }
)
)

temp_Y <- Y
fitted_trend <- list()
for ( k in 1:3){
  for (j in 1:length(idx)){
    res <- cal_Bhat_Shat(temp_Y,X )
    #plot(res$Bhat[idx[j],])
    s =susieR:: susie_trendfilter(res$Bhat[idx[j],],
                                  order=0,
                                  L=L_points)
    #plot(predict(s))
    fitted_trend[[j]] <- predict(s)
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

plot( f1, type="l", main="Estimated effect 1, sd=1",
      xlab="", ylim=c(-1.50,2.5),lwd=2   )


lines( out$fitted_func[[1]] ,col='blue',lwd=1.5   )

lines(fitted_trend[[1]],col='red' ,lwd=1.5 )
source("D:/Document/Serieux/Travail/Package/susiF.alpha/inst/code/fit_hmm.R" )
tt <- fit_hmm(res$Bhat[idx[1],],res$Shat[idx[1],],100)
lines(apply(tt$prob[,-1],1, sum), lwd=1.5, col="darkgreen")
lines( tt$prob%*%tt$mu/ susiF.obj$csd_X[idx[2]], lwd=1.5, col="green")
 legend(x= 80,
        y=2.5 ,
        lty= rep(1.5,2),
         legend = c("Truth",
                    "wave est ",
                    "Susie trend est",
                    'HMM grid est',
                    "prob neq 0 HMM"),
        col=c("black",
              'blue',
              "red",
              "green",
              "darkgreen" ),
        bty = "n"
 )


plot( f2, type="l", main="Estimated effect 2, sd=1", xlab="", ylim=c(-0.50,2.5),lwd=1.5  )


lines( out$fitted_func[[2]] ,col='blue' ,lwd=1.5)

lines(fitted_trend[[2]],col='red' ,lwd=1.5)

tt <- fit_hmm(res$Bhat[idx[2],],res$Shat[idx[2],],100)
lines(apply(tt$prob[,-1],1, sum), lwd=2, col="darkgreen")
lines( tt$prob%*%tt$mu/ susiF.obj$csd_X[idx[2]], lwd=2, col="green")


temp_Y <- Y
res <- cal_Bhat_Shat(temp_Y,X )


