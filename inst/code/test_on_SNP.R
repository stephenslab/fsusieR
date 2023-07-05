library(susiF.alpha)
library(susieR)
library(wavethresh)
 data(N3finemapping)
attach(N3finemapping)
res <- list()
lev_res=3
N <- 100
for ( o in (length(res)+1):1000){

  L <- sample( 1:10, size=1)
  eff <- list()
  for (l in 1:L) {
    eff[[l]] <- simu_IBSS_per_level(lev_res )$sim_func
  }

  noisy.data  <- list()
  X <- N3finemapping$X[1:N,1:100]
  true_pos <- sample( 1:ncol(X), L)
  for ( i in 1:nrow(X))
  {
    noise <- rnorm(length( eff[[1]] ), sd= 1)
    noisy.data [[i]] <-    noise
    for ( l in 1:L){
      noisy.data [[i]] <-    noisy.data [[i]] +X[i,true_pos[l]]*eff[[l]]
    }


  }
  noisy.data <- do.call(rbind, noisy.data)

  Y <- noisy.data
  m1 <- susiF(Y,X, L=10)
  res[[ o]] <- c(length(m1$cs), #number of CS
                 length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
                 Reduce("+",sapply(1:length(m1$cs), function(k)
                   ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
                 )
                 ))
  print( o)
}


df_simu <- do.call(rbind, res)
sum(df_simu[,3])/sum(df_simu[,2]+df_simu[,3])
