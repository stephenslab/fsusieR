
i=1

est=list()

for ( i in 1:n.shifts){
  shifted_x <- c(x[(i*k + 1):n], x[1:(i*k)])
  wd_shifted <- wavethresh::wd(c(shifted_x,rev(shifted_x)), filter.number = filter.number, family = family, type = "station")
  D=ashr::ash( c( wd_shifted$D) , rep( 0.2,length( wd_shifted$D)  ))$result$PosteriorMean 
  
  wd_shifted$D =   D
  tt <-      wavethresh::av.basis(
    wavethresh::convert(wd_shifted),
    level = wd_shifted$nlevels - 1,  # Reconstruct at the finest level
    ix1 = 0,                         # Start index
    ix2 = 1,                         # End index
    filter = wd_shifted$filter       # Wavelet filter
  ) 
 
 
  
  est[[i]] = 0.5 *(tt[1:1024]+rev(tt[-(1:1024)]))
}

recover_x <- function(shifted_x, i) {
  n <- length(shifted_x) 
  shift_back <- (n  - i*k  ) %% n  # Ensure non-negative index
  recovered_x <- c(shifted_x[(shift_back + 1):n], shifted_x[1:shift_back])
  return(recovered_x)
}


for ( i in 1: length(est)){
  est[[i]]=recover_x(est[[i]],i)
}
plot(colMeans(do.call(rbind, est)))
