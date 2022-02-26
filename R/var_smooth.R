

wd.D = function (data, filter.number = 10, family = "DaubLeAsymm",
                 type = "wavelet", bc = "periodic", verbose = FALSE,
                 min.scale = 0, precond = TRUE) {
  l = wavethresh::wd(data = data, filter.number = filter.number,
                     family = family, type = type, bc = bc,
                     verbose = verbose, min.scale = min.scale,
                     precond = precond)
  return(l$D)
}

# @description Computes the non-decimated wavelet transform matrix for
#   a given basis.
# @param n The sample size. Must be a power of 2.
# @param filter.number Specifies the type of wavelet basis used.
# @param family Specifies the type of wavelet basis used.
# @return The NDWT matrix for the specified basis, with the entries squared.
ndwt.mat = function (n, filter.number, family) {
  J = log2(n)
  X = diag(rep(1, n))
  W = matrix(0, n * J, n)
  W = apply(X, 1, wd.D, filter.number = filter.number, family = family )
  return(list(W2 = W^2))
}

x <- Bhat2[1,]

x.var <- (Shat2[1,])^2

x.w.d = wavethresh::wd(x, filter.number = filter.number,
                       family = family )



n <- length(x)

data <- x
data.var <- x.var
J = log2(n)
X = diag(rep(1, n))
W = matrix(0, n * J, n)

family="DaubLeAsymm"

filter.number = 10

Wl = ndwt.mat(n, filter.number = filter.number, family = family)


var.smooth = function (data, data.var, x.var.ini, basis, v.basis, Wl,
                       filter.number, family, post.var, ashparam, jash,
                       weight, J, n, SGD) {
  wmean = matrix(0, J, n)
  wvar = matrix(0, J, n)

    x.w = wavethresh::wd(data, filter.number = filter.number,
                         family = family, type = "station")

    # Diagonal of W*V*W'.
    x.w.v = apply((rep(1, n * J) %o% data.var) * Wl$W2, 1, sum)
    x.pm = rep(0, n)
    x.w.v.s = rep(0, n * J)
    for (j in 0:(J - 1)) {
      index = (((J - 1) - j) * n + 1):((J - j) * n)
      x.w.j = accessD(x.w, j)
      x.w.v.j = x.w.v[index]
      ind.nnull = (x.w.v.j != 0)
      zdat.ash = shrink.wc(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),
                           ashparam, jash = jash,
                           df = min(50, 2^(j + 1)), SGD = SGD)
      x.pm[ind.nnull] = get_pm(zdat.ash)
      x.pm[!ind.nnull] = 0
      if ((sum(is.na(x.pm)) > 0) & (SGD == TRUE)) {
        zdat.ash =
          shrink.wc(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),
                    ashparam, jash = jash,
                    df = min(50, 2^(j + 1)), SGD = FALSE)
        x.pm[ind.nnull] = get_pm(zdat.ash)
        x.pm[!ind.nnull] = 0
      }
      x.w = putD(x.w, j, x.pm)
      if (post.var == TRUE) {
        x.w.v.s[index[ind.nnull]] = get_psd(zdat.ash)^2
        x.w.v.s[index[!ind.nnull]] = 0
      }
    }
    var.est = AvBasis(convert(x.w))
    if (post.var == TRUE) {
      mv.wd = wd.var(rep(0, n), filter.number = basis$filter.number,
                     family = basis$family, type = "station")
      mv.wd$D = x.w.v.s
      var.est.var = AvBasis.var(convert.var(mv.wd))
    }

  if (post.var == TRUE) {
    return(list(var.est = var.est, var.est.var = var.est.var))
  } else {
    return(var.est)
  }
}




















n=2^10
 t=1:n/n
 spike.f = function(x) (0.75*exp(-500*(x-0.23)^2) +
   1.5*exp(-2000*(x-0.33)^2) +
   3*exp(-8000*(x-0.47)^2) +
   2.25*exp(-16000*(x-0.69)^2) +
   0.5*exp(-32000*(x-0.83)^2))
 mu.s = spike.f(t)

 # Gaussian case
 # -------------
 mu.t=(1+mu.s)/5
 plot(mu.t,type='l')
 var.fn = (0.0001+4*(exp(-550*(t-0.2)^2) +
   exp(-200*(t-0.5)^2) + exp(-950*(t-0.8)^2)))/1.35
 plot(var.fn,type='l')
 rsnr=sqrt(5)
 sigma.t=sqrt(var.fn)/mean(sqrt(var.fn))*sd(mu.t)/rsnr^2
 X.s=rnorm(n,mu.t,sigma.t)
 mu.est.rmad<-ti.thresh(X.s,method='rmad')
 mu.est.smash<-ti.thresh(X.s,method='smash')
 plot(mu.t,type='l')
 lines(mu.est.rmad,col=2)
 lines(mu.est.smash,col=4)





 ti.thresh = function (x, sigma = NULL, method = "smash", filter.number = 1,
                       family = "DaubExPhase", min.level = 3,
                       ashparam = list()) {
   n = length(x)
   J = log2(n)
   if (length(sigma) == 1)
     sigma = rep(sigma, n)

   if (is.null(sigma) & method == "rmad") {
     lambda.thresh = 1
   } else {
     lambda.thresh = 2
   }

   if (is.null(sigma)) {
     if (method == "smash") {
       sigma = sqrt(smash.gaus(x, v.est = TRUE, v.basis = TRUE,
                               filter.number = filter.number,
                               family = family, ashparam = ashparam,
                               weight = 1))
     } else if (method == "rmad") {
       x.w = wavethresh::wd(x, filter.number = filter.number,
                            family = family, type = "station")
       win.size = round(n/10)
       odd.boo = (win.size%%2 == 1)
       win.size = win.size + (1 - odd.boo)
       sigma = runmad(accessD(x.w, J - 1), win.size,
                      endrule = "func")
     } else {
       stop("Error: Method not recognized")
     }
   }

   if (filter.number == 1 & family == "DaubExPhase") {
     tsum = sum(x)
     vdtable = cxxtitable(x)$difftable
     vtable = cxxtitable(sigma^2)$sumtable
     wmean = threshold.haar(vdtable, vtable, lambda.thresh, levels = 0:(J - 1 - min.level), type = "hard")
     wwmean = -wmean
     mu.est = cxxreverse_gwave(tsum, wmean, wwmean)
   } else {
     x.w = wavethresh::wd(x, filter.number = filter.number,
                          family = family, type = "station")
     Wl = ndwt.mat(n, filter.number = filter.number, family = family)

     # Diagonal of W*V*W'.
     x.w.v = apply((rep(1, n * J) %o% (sigma^2)) * Wl$W2, 1, sum)
     x.w.t = threshold.var(x.w, x.w.v, lambda.thresh,
                           levels = (min.level):(J - 1),
                           type = "hard")
     mu.est = AvBasis(convert(x.w.t))
   }
   return(mu.est)
 }





 # This function performs "decomposition" of variances of detail
 # coefficients for a given wavelet basis in wavelet transformation.
 #
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






 lol0 <- c(rWShat2)
 lol1 <- c(Shat)

 plot(Shat[10,])
 plot(rWShat2[10,])


 which(Shat == max(Shat), arr.ind = TRUE)
 which(rWShat2 == max(rWShat2), arr.ind = TRUE)

 idx_perm <- c()
 for ( i in 1:length(lol0))
 {
   print(which(lol1== lol0[i]))
   idx_perm <- c( which(lol1== lol0[i]),idx_perm)
 }


 tt <- wd.var(rep(0, 1024))
 plot(tt$C)

 W1 <- (GenW(n=  ncol(Shat2)  , filter.number = 10, family = "DaubLeAsymm"))
 W1[, c(1,ncol(W1 ))]<- W1[, c(ncol(W1),1)]

 Wl<- W1
 # Diagonal of W*V*W'.


 x.w.v = apply((rep(1, n * J) %o% data.var) * W1^2, 1, sum)
 x.pm = rep(0, n)
 x.w.v.s = rep(0, n * J)
