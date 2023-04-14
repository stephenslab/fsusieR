set.seed(1)
library(susiF.alpha)
library(haarfisz)
func <- wavethresh::DJ.EX(n = 512, noisy = FALSE)
baseline <-  func[[1]]-min(func[[1]])
plot( baseline, type="l")
effect <-   func[[2]]-min(func[[2]])
plot(effect, type="l")


N=50
nullweight = 10/sqrt(N)
P=10
G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
beta2       <- 1
noisy.data  <- list()

for ( i in 1:N)
{

  noisy.data [[i]] <- rpois(n=length(baseline), lambda = (baseline+G[i,1]*effect))

}
noisy.data <- do.call(rbind, noisy.data)

pos1<-1


plot( noisy.data[1,], type = "l", col=(G[1, pos1]*3+1),
      main="Observed curves \n colored by the causal effect", ylim= c(0,100), xlab="")
for ( i in 2:N)
{
  lines( noisy.data[i,], type = "l", col=(G[i, pos1]*3+1))

}
legend(x=0.3,
       y=-10,
       lty = rep(1,3),
       legend= c("0", "1","2"),
       col=c("black","blue","yellow"))



Y <- noisy.data
Y_0 <- noisy.data
X <- G
library(wavethresh)
library(haarfisz)
Y_f <-DWT2_pois(Y)
Y <- Y_f




plot(Y_0[2,],hft.inv(Y_f[2,]))


#Y <- colScale(Y, scale=FALSE)
#X <- colScale(X)


#Using a column like phenotype

v1 <- rep(1, nrow(X))
tt <-  cal_Bhat_Shat2(Y_f,X,v1,lowc_wc = NULL)


indx_lst <- gen_wavelet_indx(9)
B0hat <- tt$B0hat
Bhat <- tt$Bhat
Shat <- tt$Shat
init_pi0_w=1
control_mixsqp =  list(verbose=FALSE)

tt <-  ash(Y_f[1,], rep( 0.1, 512))

k=0
plot( hft.inv( Bhat[1, ]/0.66 )  , ylim=c(0,70), type="l", col="green")
abline(a=0,b=0)
lines(  effect)

plot((   hft.inv( Bhat[1, ]/var(X[,1]) )  ),effect)
abline(a=0,b=1)


tt2 <- ash(Bhat[1, ],( Shat[1, ]))
plot(hft.inv(tt2$result$PosteriorMean/var(X[,1]))  , ylim=c(0,70), type="l", col="blue")

lines(  effect)
abline(a=1,b=0)








plot(hft.inv( B0hat[1, ]) , ylim=c(0,70), type="l", col="blue")
abline(a=1,b=0)
lines(  baseline)
