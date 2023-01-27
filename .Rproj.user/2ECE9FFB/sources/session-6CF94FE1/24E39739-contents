
rm(list=ls())
library(susiF.alpha)
N <- 40
lev_res <- 7
temp <-  wavethresh::DJ.EX(n=2^lev_res,noisy = FALSE)

idx <- sample(1:4,
              replace = TRUE,
              size=N)
Batch <- paste("Batch",
               idx ,
               sep="_")

#length of the molecular phenotype (2^lev_res)
f1 <-  simu_IBSS_per_level(lev_res )$sim_func
sim_func <- list()
for ( i in 1:N){

  sim_func[[i]] <-  unlist((temp[idx[i]]))+ 1*f1+ rnorm(2^lev_res,sd=1)
}
df <- do.call( rbind, sim_func)
Batch <- factor(Batch)
temp <- data.frame(Batch =Batch)
 ff <- ~-1+  Batch  #ensure that each Batch effect is estimated individually
utils::str(m <- model.frame(ff, temp))
mat <- model.matrix(ff, m)

mat
Cov <-   mat
adjust_data <- adjust_FM_covariate(Y=df ,X=Cov )


temp <-  wavethresh::DJ.EX(n=2^lev_res,noisy = FALSE)
plot(adjust_data$fitted_coef[,1], type="l")
lines(unlist(temp[[1]]))

plot(adjust_data$fitted_coef[,2], type="l")
lines(unlist(temp[[2]]))
plot(adjust_data$fitted_coef[,3], type="l")
lines(unlist(temp[[3]]))

plot(adjust_data$fitted_coef[,4], type="l")
lines(unlist(temp[[4]]))




par(mfrow=c(2,1))
plot(apply(df[which( Batch=="Batch_2"),],2,mean) - apply(df[which( Batch=="Batch_1"),],2,mean)  ,
     lwd=2, type="l",col="orange"
  )
lines(c(apply(df[which( Batch=="Batch_3"),],2,mean) - apply(df[which( Batch=="Batch_1"),],2,mean)),
      lwd=2,
      col="green")

lines(c(apply(df[which( Batch=="Batch_4"),],2,mean) - apply(df[which( Batch=="Batch_1"),],2,mean))  ,
      lwd=2,
      col="blue"  )






plot(apply(adjust_data$Y_adjusted[which( Batch=="Batch_2"),],2,mean) - apply(adjust_data$Y_adjusted[which( Batch=="Batch_1"),],2,mean)  ,
     lwd=2, type="l",col="orange", ylim=c(-20,20)
)
lines(c(apply(adjust_data$Y_adjusted[which( Batch=="Batch_3"),],2,mean) - apply(adjust_data$Y_adjusted[which( Batch=="Batch_1"),],2,mean)),
      lwd=2,
      col="green")

lines(c(apply(adjust_data$Y_adjusted[which( Batch=="Batch_4"),],2,mean) - apply(adjust_data$Y_adjusted[which( Batch=="Batch_1"),],2,mean))  ,
      lwd=2,
      col="blue"  )




