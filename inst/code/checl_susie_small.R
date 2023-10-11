
library(susiF.alpha)
library(susieR)
set.seed(10)
data(N3finemapping)
attach(N3finemapping)


N =30
if(file.exists("D:/Document/Serieux/Travail/Package/susiF.alpha/pb_susie.RData")){
  load("D:/Document/Serieux/Travail/Package/susiF.alpha/pb_susie.RData")
}else{
  res <- list()
}


for (o in 1:10000){

  G <- N3finemapping$X[sample(1:nrow(N3finemapping$X), size=N, replace= FALSE),]



  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }

  L <- sample (1:10, size=1)
  # G <- matrix( rnorm(N*300), nrow = N)
  true_pos <- sample( 1:ncol(G), L)

  beta <- 5
  Y <- rep(NA, N)
  for( i in 1:N){
    Y[i] <- sum( beta*G[,true_pos])+ rnorm(1,sd=0.51)
  }

  m1 <- susie( X,Y, L=10)
   if (length(m1$sets$cs)==0){
     out <- NULL
   }else{
     out <-   c( length(m1$sets$cs), #number of CS
                 length(which(true_pos%in% do.call(c, m1$sets$cs))), #number of effect found
                 Reduce("+",sapply(1:length(m1$sets$cs), function(k)
                   ifelse( length(which(true_pos%in%m1$sets$cs[[k]] ))==0, 1,0)
                 )
                 )
     )
     res[[o]] <- out
   }


    save( res, file="D:/Document/Serieux/Travail/Package/susiF.alpha/pb_susie.RData")
    print(res)

}




