################################## WAVELET UTILITY FUNCTIONS ############################

# @title Lifting scheme for non decimated and unevenly spaced data on matrix
#
# @description  Interpolation procedure from Kovac and Silvermann 2000 for a matrix of functions
#
#
# @details We suppose that we observe n functions/curves at T different time points. We remap these functions into a grid of length equal to a power of two using the lifting scheme of Kovac and Silvermann 2000.
#
# @param Y matrix of curves of size nXT, where the curves are stored by row
#
# @param pos vector of positive values, each the vector component corresponds to the measurement time. (the values of pos should increase along with the index)
#
# @importFrom  wavethresh makegrid
#
# @export
#
interpol_mat <- function(Y, pos, max_scale=10)
{

  bp    <- (pos- min(pos))/(max(pos)-min(pos))
  Y_new <- t(apply(Y,
                   1, 
                   interpolKS2, 
                   bp=bp))
  grid  <- wavethresh::makegrid(t=bp,
                                y = 1:dim(Y)[2],
                                gridn = min( 2^max_scale, 
                                             2^(floor(log(length(pos)-1,2)) + 1)   )
                                )$gridt
    grid <-  (grid  - min(grid) )*length(grid)/(max(grid)- min(grid))  #* (max(bp)-  min(bp))/ (max(grid)-  min(grid))
  out <- list(Y    = Y_new,
              grid = grid,
              bp   = bp)
  return(out)
}




interpolKS2 <-  function (y, bp, max_scale=10)
{
  wavethresh::makegrid(t=bp,y=y, gridn = min( 2^max_scale, 2^(floor(log(length(bp)-1,2)) + 1)   ))$gridy
}



# @title Interpolation procedure from Kovac and Silverman 2000
#
# @description  Interpolation procedure from Kovac and Silveramnn 2000
#
# @param y vector of observed curves
#
# @param bp scaled observation time points on 0,1
#
# @importFrom wavethresh makegrid
#
# @export
interpolKS <-  function (y, bp)
{
  wavethresh::makegrid(t=bp,y=y, gridn = 2^(floor(log(length(bp)-1,2)) + 1)   )$gridy
}


# @title  wavelet transform on a matrix of functions
#
# @description function adapted from grove R package from Ma and Soriano. Perform wavelet transform of each row of a matrix
#
# @param data matrix of size NxJ where J is a power of two
#
# @param filter.number This selects the smoothness of the wavelet that you want to use in the decomposition. By default, this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments. Description from Guy Nason (wavethresh package)
#
# @param family specifies the family of wavelets that you want to use. Two popular options are "DaubExPhase" and "DaubLeAsymm" but see the help for filter.select for more possibilities.. Description from Guy Nason (wavethresh package)
#
# @return A list with the following components
#
#\item{C}{ vector length n containing the wavelet C coefficient}
#\item{D}{ matrix of size nx 2^(J -1) where each row contains the wavelet D coefficients, ordered in the same way as in the wavethresh package}
#\item{family}{ used for the wavelet transform}
#\item{filter.number}{ }
#
# @export
#
#' @importFrom stats complete.cases
#' @importFrom wavethresh accessC
#' @importFrom wavethresh wd
#'
DWT2 <- function (data, filter.number = 10, family = "DaubLeAsymm" )
{

  NA_pos <- which(!complete.cases(data))
  if (length(NA_pos )>0){
    data [is.na(data)]<-0
  }
  J <- ncol(data)
  n <- nrow(data)
  D <- matrix(NA, nrow = n, ncol = J - 1)
  C <- rep(NA, n)
  for (i in 1:n) { ## Speed Gain

      temp <- wavethresh::wd(data[i, ],
                             filter.number = filter.number,
                             family = family,
                             min.scale=10 )

    D[i, ] <- temp$D
    C[i] <- wavethresh::accessC(temp, level = 0)
  }
  if (length(NA_pos )>0){
    D [NA_pos,] <- NA
    C [NA_pos]  <- NA
  }

  output <- list(C = C, D = D, J = log2(J), filter.number = filter.number,
                 family = family)
  class(output) <- "DWT"
  return(output)
}


#' @title Generate list of index wavelet coefficients
#
#' @description Generate list of index corresponding to where the wavelet coefficient of scale s are stored in wd$D.
#First element of the list corresponds to the index for the wavelet coefficient at scale 0, the second component corresponds to the index for the wavelet coefficient at scale 1 (and so on). The last componnet of the list correspond to where the C coefficients is stored in init_prior for class mixture_normal_per_scale
#
#' @param lev_res integer, corresponding to log2 of the signal length. WARNING the ordering change for different values of s.
#
#' @export
#' @keywords internal
#' @examples
#' library(wavethresh)
#' tem_func <- rnorm( 2^8)
#' twav <- wd(tem_func)
#' indx_lst <- gen_wavelet_indx(8)
#' plot(accessD(twav,level=6), ( twav$D[unlist(indx_lst[(6+1)])]) )
#'  #should a straightline
#'
gen_wavelet_indx <- function(lev_res)
{

  if( lev_res>15){
    stop("lev_res is larger than 15 please use a smaller level of resolution")
  }

  indx_lst <- list()
  indx_lst[[1]] <- 2^lev_res -1 #coefficient
  for ( s in 1:(lev_res-1))
  {

    indx  <- 2^(lev_res)- (2^((s+1))-1) :(2^s)


    indx_lst[[s+1]] <- indx
  }
  indx_lst[[length(indx_lst)+1]] <- 2^lev_res# C coefficient
  out <- indx_lst
  return(out)

}

#' @importFrom stats lm
wavelet_reg <-  function(Y, X, verbose=TRUE,  pos=NULL,   thresh_lowcount=0){



  if (is.null(pos))
  {
    pos <- 1:dim(Y)[2]
  }
  if(!is.wholenumber(log2(dim(Y)[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check whether dim(Y) not equal to 2^J or if the data are unevenly spaced
  {

    inter_pol.obj <- interpol_mat(Y, pos)
    Y             <- inter_pol.obj$Y
    bp            <- inter_pol.obj$bp

    start_pos <- min( pos)
    end_pos <-max(pos)
    outing_grid   <- start_pos + (end_pos-start_pos)/(length(pos))*inter_pol.obj$grid
    if(verbose)
    {
      message( "Response matrix dimensions not equal to nx 2^J \n or unevenly spaced data \n interpolation procedure used")
    }
  }  else{

    outing_grid   <- pos
  }


  W <- DWT2(Y)
  Y_f      <-  cbind( W$D,W$C)


  v1 <- rep(1, dim(X)[2])
  lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)


  Est <- lapply( 1:ncol(Y_f), function(j){
    fit <-  lm(Y_f[,j]~ X)
    out <-  list(coef=summary(fit)$coefficients[,1] ,
                 sd  =summary(fit)$coefficients[,2])
    return(out)
  }




  )
  Coef <- t(do.call(rbind,lapply(1:length(Est), function(j) Est[[j]]$coef)))
  Sd   <- t(do.call(rbind,lapply(1:length(Est), function(j) Est[[j]]$sd)))



  indx_lst <- gen_wavelet_indx(log2(dim(Y_f)[2]))
  if (is.null(nrow(Coef))){



    coef_p <- do.call(c,lapply( 1:length(indx_lst),
                      function (j) ashr::ash(Coef[indx_lst[[j]]],
                                             Sd  [indx_lst[[j]]]
                                              )$result$PosteriorMean
                      )
                    )
    temp <- wavethresh::wd(rep(0, length(Coef)[2]))

    temp$D                     <-  coef_p[ -indx_lst[[length(indx_lst)]]]
    temp$C[length(temp$C)]     <-  coef_p[ indx_lst[[length(indx_lst)]]]
    fitted_coef   <-  wavethresh::wr(temp)
  }else{
    fitted_coef <-list()
    for ( i in 1:nrow(Coef ))
    {
      coef_p <-  do.call(c,lapply( 1:length(indx_lst),
                                   function (j) ashr::ash(Coef[i,indx_lst[[j]]],
                                                          Sd  [i,indx_lst[[j]]]
                                   )$result$PosteriorMean
                                  )
                            )



      temp <- wavethresh::wd(rep(0, dim(Coef)[2]))

      temp$D                     <-  coef_p[  -indx_lst[[length(indx_lst)]]]
      temp$C[length(temp$C)]     <-   coef_p[  indx_lst[[length(indx_lst)]]]
      fitted_coef[[i]]   <-  wavethresh::wr(temp)
    }
    fitted_coef <- do.call(cbind,   fitted_coef )
  }


  Y_adjusted  = Y-X%*%t(fitted_coef[,-1])#removing estimated baseline
  out <- list( Y_adjusted  = Y_adjusted,

               fitted_coef = fitted_coef[,-1],#removing estimated baseline
               pos         = outing_grid
  )


  return(out)
}



#'@title  Ramp matrix of unevenly space data (or non power of two)
#
#' @description  Ramp matrix of unevenly space data (or non power of two)
#' @param Y  Matrix of observed curves
#' @param pos sampling position
#' @param verbose logical
#'@param max_scale numeric, define the maximum of wavelet coefficients used in the analysis (2^max_scale).
#'        Set 10 true by default.
#' @importFrom stats complete.cases
#' @export
remap_data <- function(Y,pos, verbose=TRUE, max_scale=10){

  NA_pos <- which(!complete.cases(Y))
  if (length(NA_pos )>0){
    Y [is.na(Y)]<-0
  }

  if(!is.wholenumber(log2(dim(Y)[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check whether dim(Y) not equal to 2^J or if the data are unevenly spaced
  {

    inter_pol.obj <-interpol_mat(Y, pos, max_scale =  max_scale)
    Y             <- inter_pol.obj$Y
    bp            <- inter_pol.obj$bp

    start_pos <- min( pos)
    end_pos <-max(pos)
    outing_grid   <- start_pos + (end_pos-start_pos)/(length(inter_pol.obj$grid))*inter_pol.obj$grid
    if(verbose)
    {
      message( "Response matrix dimensions not equal to nx 2^J \n or unevenly spaced data \n interpolation procedure used")
    }
  }  else{

    outing_grid   <- pos
  }
  if (length(NA_pos )>0){
    Y [NA_pos,]<-NA
  }
  return(list(Y           = Y,
              outing_grid = outing_grid))
}



# This function reconstructs variance of the mean estimate for a
# given wavelet basis given a wst object.
# code from the smashr package
#
#' @importFrom wavethresh av.basis
#' @importFrom wavethresh nlevelsWT
AvBasis.var <- function (wst, Ccode = TRUE, ...) {
  nlevels <- nlevelsWT(wst)
  
  if (is.null(wst$filter$G)) {
    if (Ccode == FALSE) {
      answer <- av.basis(wst, level = nlevels - 1, ix1 = 0,
                         ix2 = 1, filter = wst$filter)
    }
    else {
      error <- 0
      answer <- rep(0, 2^nlevels)
      H <- wst$filter$H
      aobj <-
        .C("av_basisWRAP", wstR = as.double(wst$wp),
           wstC = as.double(wst$Carray),
           LengthData = as.integer(length(answer)),
           level = as.integer(nlevels - 1), H = as.double(H),
           LengthH = as.integer(length(H)), answer = as.double(answer),
           error = as.integer(error), PACKAGE = "wavethresh")
      if (aobj$error != 0)
        stop(paste("av_basisWRAP returned error code",
                   aobj$error))
      answer <- aobj$answer
    }
  }
  else {
    error <- 0
    answerR <- answerI <- rep(0, 2^nlevels)
    H <- wst$filter$H
    G <- wst$filter$G
    aobj <- .C("comAB_WRAP", wstR = as.double(Re(wst$wp)),
               wstI = as.double(Im(wst$wp)),
               wstCR = as.double(Re(wst$Carray)),
               wstCI = as.double(Im(wst$Carray)),
               LengthData = as.integer(length(answerR)),
               level = as.integer(nlevels - 1), HR = as.double(Re(H)),
               HI = as.double(Im(H)), GR = as.double(Re(G)),
               GI = as.double(Im(G)), LengthH = as.integer(length(H)),
               answerR = as.double(answerR), answerI = as.double(answerI),
               error = as.integer(error),
               PACKAGE = "wavethresh")
    if (aobj$error != 0)
      stop(paste("av_basisWRAP returned error code", aobj$error))
    answer <- aobj$answerR
  }
  answer
}



# This function performs "decomposition" of variances of detail
# coefficients for a given wavelet basis in wavelet transformation.
#code from the smashr package
#' @importFrom stats tsp
#' @importFrom stats tsp<-
#' @importFrom wavethresh filter.select
#' @importFrom wavethresh first.last
#' @importFrom wavethresh wd.int
#' @importFrom wavethresh IsPowerOfTwo
wd.var <- function (data, filter.number = 10, family = "DaubLeAsymm",
                    type = "wavelet", bc = "periodic", verbose = FALSE,
                    min.scale = 0, precond = TRUE) {
  if (verbose == TRUE)
    cat("wd: Argument checking...")
  if (!is.atomic(data))
    stop("Data is not atomic")
  DataLength <- length(data)
  nlevels <- nlevelsWT(data)
  filter <- wavethresh::filter.select(filter.number, family)
  filter$H <- wavethresh::filter.select(filter.number,family)$H^2
  filter$G <- wavethresh::filter.select(filter.number,family)$H^2
  if (is.na(nlevels))
    stop("Data length is not power of two")
  if (type != "wavelet" && type != "station")
    stop("Unknown type of wavelet decomposition")
  if (type == "station" && bc != "periodic")
    stop("Can only do periodic boundary conditions with station")
  if (verbose == TRUE)
    cat("...done\nFilter...")
  if (verbose == TRUE)
    cat("...selected\nFirst/last database...")
  fl.dbase <- wavethresh::first.last(LengthH = length(filter$H),
                                     DataLength = DataLength,
                                     type = type, bc = bc)
 
  if (bc == "interval") {
    ans <- wavethresh::wd.int(data = data,
                              preferred.filter.number = filter.number,
                              min.scale = min.scale,
                              precond = precond)
    fl.dbase <- wavethresh::first.last(LengthH = length(filter$H),
                                       DataLength = DataLength,
                                       type = type, bc = bc, current.scale = min.scale)
    filter <- list(name = paste("CDV", filter.number, sep = ""),
                   family = "CDV", filter.number = filter.number)
    l <-
      list(transformed.vector = ans$transformed.vector,
           current.scale = ans$current.scale,
           filters.used = ans$filters.used,
           preconditioned = ans$preconditioned, date = ans$date,
           nlevels =
             wavethresh::IsPowerOfTwo(length(ans$transformed.vector)),
           fl.dbase = fl.dbase, type = type, bc = bc, filter = filter)
    class(l) <- "wd"
    return(l)
  }
  dtsp <- tsp(data)
  C <- rep(0, fl.dbase$ntotal)
  C[1:DataLength] <- data
  if (verbose == TRUE)
    error <- 1
  else error <- 0
  if (verbose == TRUE)
    cat("built\n")
  if (verbose == TRUE)
    cat("Decomposing...\n")
  nbc <- switch(bc, periodic = 1, symmetric = 2)
  if (is.null(nbc))
    stop("Unknown boundary condition")
  ntype <- switch(type, wavelet = 1, station = 2)
  if (is.null(filter$G)) {
    wavelet.decomposition <-
      .C("wavedecomp", C = as.double(C),
         D = as.double(rep(0, fl.dbase$ntotal.d)),
         H = as.double(filter$H),
         LengthH = as.integer(length(filter$H)),
         nlevels = as.integer(nlevels),
         firstC = as.integer(fl.dbase$first.last.c[, 1]),
         lastC = as.integer(fl.dbase$first.last.c[, 2]),
         offsetC = as.integer(fl.dbase$first.last.c[,3]),
         firstD = as.integer(fl.dbase$first.last.d[,
                                                   1]), lastD = as.integer(fl.dbase$first.last.d[,
                                                   2]), offsetD = as.integer(fl.dbase$first.last.d[,
                                                   3]), ntype = as.integer(ntype), nbc = as.integer(nbc),
         error = as.integer(error), PACKAGE = "wavethresh")
  }
  else {
    wavelet.decomposition <-
      .C("comwd", CR = as.double(Re(C)),
         CI = as.double(Im(C)), LengthC = as.integer(fl.dbase$ntotal),
         DR = as.double(rep(0, fl.dbase$ntotal.d)),
         DI = as.double(rep(0, fl.dbase$ntotal.d)),
         LengthD = as.integer(fl.dbase$ntotal.d),
         HR = as.double(Re(filter$H)), HI = as.double(-Im(filter$H)),
         GR = as.double(Re(filter$G)), GI = as.double(-Im(filter$G)),
         LengthH = as.integer(length(filter$H)),
         nlevels = as.integer(nlevels),
         firstC = as.integer(fl.dbase$first.last.c[, 1]),
         lastC = as.integer(fl.dbase$first.last.c[, 2]),
         offsetC = as.integer(fl.dbase$first.last.c[,3]),
         firstD = as.integer(fl.dbase$first.last.d[,1]),
         lastD = as.integer(fl.dbase$first.last.d[,2]),
         offsetD = as.integer(fl.dbase$first.last.d[,3]),
         ntype = as.integer(ntype), nbc = as.integer(nbc),
         error = as.integer(error), PACKAGE = "wavethresh")
  }
  if (verbose == TRUE)
    cat("done\n")
  error <- wavelet.decomposition$error
  if (error != 0) {
    cat("Error ", error, " occured in wavedecomp\n")
    stop("Error")
  }
  if (is.null(filter$G)) {
    l <- list(C = wavelet.decomposition$C, D = wavelet.decomposition$D,
              nlevels = nlevelsWT(wavelet.decomposition), fl.dbase = fl.dbase,
              filter = filter, type = type, bc = bc, date = date())
  }
  else {
    l <- list(C = complex(real = wavelet.decomposition$CR,
                          imaginary = wavelet.decomposition$CI),
              D = complex(real = wavelet.decomposition$DR,
                          imaginary = wavelet.decomposition$DI),
              nlevels = nlevelsWT(wavelet.decomposition),
              fl.dbase = fl.dbase, filter = filter, type = type,
              bc = bc, date = date())
  }
  class(l) <- "wd"
  if (!is.null(dtsp))
    tsp(l) <- dtsp
  return(l)
}


# This function converts variances of detail coefficients for a given
# wavelet basis from wd objects to wst objects.
#code from the smashr package
#' @importFrom wavethresh wst
#' @importFrom wavethresh getarrvec
#' @importFrom wavethresh accessD
#' @importFrom wavethresh accessC
#' @importFrom wavethresh putD
#' @importFrom wavethresh putC
convert.var <- function (wd, ...) {
  if (wd$type != "station")
    stop("Object to convert must be of type \"station\" ")
  n <- 2^nlevelsWT(wd)
  dummy <- rep(0, n)
  tmpwst <- wavethresh::wst(dummy, filter.number = wd$filter$filter.number,
                            family = wd$filter$family)
  tmpwst$filter <- wd$filter
  tmpwst$date <- wd$date
  arrvec <- wavethresh::getarrvec(nlevelsWT(wd), sort = FALSE)
  for (lev in (nlevelsWT(wd) - 1):1) {
    ds <- wavethresh::accessD.wd(wd, level = lev)
    cs <- wavethresh::accessC.wd(wd, level = lev)
    ds <- ds[arrvec[, nlevelsWT(wd) - lev]]
    cs <- cs[arrvec[, nlevelsWT(wd) - lev]]
    tmpwst <- wavethresh::putD(tmpwst, level = lev, v = ds)
    tmpwst <- wavethresh::putC(tmpwst, level = lev, v = cs)
  }
  tmpwst <- wavethresh::putC(tmpwst, level = nlevelsWT(wd),
                             v = accessC(wd, level = wd$nlevels))
  tmpwst <- wavethresh::putD(tmpwst, level = nlevelsWT(wd),
                             v = accessC(wd, level = wd$nlevels))
  tmpwst <- wavethresh::putC(tmpwst, level = 0,
                             v = wavethresh::accessC(wd, level = 0))
  arrvec <- sort.list(wavethresh::levarr(1:n, levstodo = nlevelsWT(wd)))
  tmpwst <- wavethresh::putD(tmpwst, level = 0,
                             v = wavethresh::accessD(wd, level = 0)[arrvec])
  return(tmpwst)
}






TI_ash_smooth=  function( betahat,sds ,n.shifts=10 ,family =  "DaubExPhase",
                          filter.number=10){
  
  
  x=betahat
  
  n <- length(x)
   
  
  n <- length(x)
  n.shifts=10  
  family =  "DaubExPhase" 
  filter.number=10
  # Initialize a matrix to store smoothed signals for each shift
  smoothed_signals <- matrix(0, nrow = n.shifts, ncol = n)
  smoothed_var <- matrix(0, nrow = n.shifts, ncol = n)
  
  
  k= floor(n/n.shifts)
  
  i=1
  shifted_x <- c(x[(i*k + 1):n], x[1:(i*k)])
  wd_shifted <- wavethresh::wd(shifted_x, filter.number = filter.number, family = family, type = "station")
  
  mat_Coef_D=  matrix( 0, nrow = n.shifts, ncol=  length( wd_shifted$D))
  mat_Coef_D_var = matrix( 0, nrow = n.shifts, ncol=  length( wd_shifted$D))
  
  shrunk_wc   <- matrix( 0,ncol=length( wd_shifted$D) , nrow=n.shifts,
                         byrow = FALSE)
  shrunk_var  <- matrix( 0,ncol=length( wd_shifted$D) , nrow=n.shifts,
                         byrow = FALSE)
  for (i in 1:n.shifts) {
    
    # Shift the signal
    shifted_x <- c(x[(i*k + 1):n], x[1:(i*k)])
    
    # Perform wavelet decomposition
    wd_shifted <- wavethresh::wd(shifted_x, filter.number = filter.number, family = family, type = "station")
    mat_Coef_D[i,]=wd_shifted$D
    
    res_ash= ashr::ash( c( mat_Coef_D[i,]) , rep( sds,length( wd_shifted$D)  ))
    
    
    shrunk_wc[i, ]   <-  res_ash$result$PosteriorMean 
    shrunk_var[i, ]  <-   res_ash$result$PosteriorSD^2 
  }
  
  
  recover_x <- function(shifted_x, i) {
    n <- length(shifted_x)
    n_tilt=n#+sample(c(-1,0,1),size=1)
    shift_back <- (n_tilt - i*k) %% n  # Ensure non-negative index
    recovered_x <- c(shifted_x[(shift_back + 1):n], shifted_x[1:shift_back])
    return(recovered_x)
  }
  
  for ( i in 1:n.shifts){
    
    # Reconstruct the smoothed signal using av.basis
    wd_shifted$D =   shrunk_wc[i, ] 
    smoothed_signals[i, ] <-      wavethresh::av.basis(
      wavethresh::convert(wd_shifted),
      level = wd_shifted$nlevels - 1,  # Reconstruct at the finest level
      ix1 = 0,                         # Start index
      ix2 = 1,                         # End index
      filter = wd_shifted$filter       # Wavelet filter
    ) 
    
    wd_shifted$D=(shrunk_var[i,])
    ab_var=  AvBasis.var(convert.var(wd_shifted))
    smoothed_var[i, ] =  recover_x(ab_var,i)
    
  }
  
  recover_x <- function(shifted_x, i) {
    n <- length(shifted_x)
    tilt= sample(c(-2:2),size=1)
    shift_back <- (n  - i*k +tilt ) %% n  # Ensure non-negative index
    recovered_x <- c(shifted_x[(shift_back + 1):n], shifted_x[1:shift_back])
    return(recovered_x)
  }
  
  for ( i in 1:n.shifts){
    smoothed_signals[i, ]= recover_x(smoothed_signals[i,],i)
  }
  for ( i in 1:n.shifts){
    smoothed_var[i, ]= recover_x(smoothed_var[i,],i)
  }
  
  return(list(smoothed_signal= colMeans(smoothed_signals),
              smoothed_var=  colMeans(smoothed_var)))
  
}





