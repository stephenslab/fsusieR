rm(list = ls()) 
library(fsusieR)
library(susieR)
library(ebnm)
set.seed(1)
'%!in%' <- function(x,y)!('%in%'(x,y))
data(N3finemapping)
X <- N3finemapping$X
mysd= 0.1
N =200
genotype <-X[1:N,1:100]
data(N3finemapping)
X <- N3finemapping$X
genotype <-X[sample(1:nrow(X), size=N),1:100]

idx <- which( apply( genotype,2, var ) <1e-15)
if( length(idx)==0){
  X <-genotype
  
  Rtrue <- cor (genotype )
}else{
  genotype <- genotype [, -idx]
  X <-genotype
  
}
G<- genotype
X <- (X -0.99*min(X))/(0.5*max(X ))

G <-  (G -0.99*min(G ))/(0.5*max(G ))



idx <- which( apply( genotype,2, var ) <1e-15)

if ( length(idx)>0){
  genotype <- genotype [, -idx]
} 
lev_res =6
count.data  <- list()
L <-2# sample(1:2, size =1)#actual number of effect

lf <-  list()
lf[[1]]<-  rep(0,2^6)# cos((1:2^6) /(0.5*2^6))# #rep(0.1, 2^6)
lf[[1]][20:30]=1.5
#lf[[1]][10:20] <-2

lf[[2]]<- sin((1:2^6) /(0.5*2^6)) # rep(0.1, 2^6)
#lf[[2]][50:60] <-2

true_pos <- sample(1:ncol(genotype), replace = FALSE,size=2)

plot (lf[[1]], type="l", main=paste ( "effect of SNP",true_pos[1]))

plot (lf[[2]], type="l", main=paste ( "effect of SNP",true_pos[2]))


pos1 <- true_pos[1]
pos2 <- true_pos[2]
if( length(which(apply(G,2,var)==0))>0){
  G <- G[,-which(apply(G,2,var)==0)]
}
# G <- matrix( rnorm(nrow(genotype)*300), nrow = nrow(genotype))


predictor <-rep (0, length(lf[[1]] ))
count.data  <- list()

G[ , true_pos[1]] <-G[ , true_pos[1]] -min(G[ , true_pos[1]] )
G[ , true_pos[2]] <-G[ , true_pos[2]] -min(G[ , true_pos[2]] )
for ( i in 1:N)
{
  
  predictor <-rep (0, length(lf[[1]] ))
  
  for ( l in 1:L){
    predictor <-predictor + G[i, true_pos[l]]*lf[[l]]+0.3
  }
  predictor <- exp( predictor+ rnorm(  length(lf[[1]]), sd=mysd))
  
  count.data [[i]] <-   rpois(n= length(lf[[1]]) ,
                              lambda =predictor  )
  
}
count.data <- do.call(rbind, count.data)


Y <- count.data


res0= susiF(Y=log1p(Y),X=X, L=3)

res= Pois_fSuSiE(Y=Y,X=X, L=3, max.iter =   10, post_processing = "smash")

plot(res0$ind_fitted_func,Y, pch=19)
points(res$Mu_pm,Y, pch=19, col="blue")


cor(c(res$Mu_pm),log1p(c(Y)))
cor(c(res0$ind_fitted_func),log1p(c(Y)))

plot(res$susiF.obj$fitted_func[[1]], type="l")
lines(res0$fitted_func[[1]], col="green")
lines(lf[[1]], lwd=2) 
plot(res$susiF.obj$fitted_func[[2]], type="l")
lines(res0$fitted_func[[1]], col="green")
lines(lf[[2]], lwd=2) 
