rm(list=ls())
devtools::load_all(".")
load("D:/Document/Serieux/Travail/Package/mvf.susie.alpha/inst/pb_input.RData")

source("D:/Document/Serieux/Travail/Package/susiF.alpha/R/operation_on_susiF_obj.R")
X <- pb_input[[3]]$X
Y <- pb_input[[3]]$Y$Y_f[[1]]
true_pos  <- pb_input[[3]]$true_pos
is.pois=FALSE
check=10

L=11

pos = NULL
prior = "mixture_normal_per_scale"
verbose = TRUE
maxit = 100
tol = 1e-3
cov_lev = 0.95
min.purity=0.5
filter.cs =TRUE
init_pi0_w= 1
nullweight =10
control_mixsqp =  list(verbose=FALSE,
                       eps = 1e-6,
                       numiter.em = 4
)
thresh_lowcount=0
cal_obj=FALSE
L_start = 3
quantile_trans=FALSE
greedy =TRUE
backfit =TRUE
gridmult= sqrt(2)
max_scale=10
max_SNP_EM=1000
max_step_EM=1
cor_small=FALSE






if(max_SNP_EM<10){
  stop("Argument max_SNP_EM has to be larger than 10")
}

####Cleaning input -----
pt <- proc.time()
if( prior %!in% c("normal", "mixture_normal", "mixture_normal_per_scale"))
{
  stop("Error: provide valid prior input")
}

if(!cal_obj){
  tol <-10^-3
}
if(L_start >L)
{
  L_start <- L
}
## Input error messages

if (is.null(pos))
{
  pos <- 1:dim(Y)[2]
}

#reshaping of the data
if ( !(length(pos)==dim(Y)[2])) #miss matching positions and number of observations
{
  stop("Error: number of position provided different from the number of column of Y")
}



map_data <- remap_data(Y=Y,
                       pos=pos,
                       verbose=verbose,
                       max_scale=max_scale)

outing_grid <- map_data$outing_grid
Y           <- map_data$Y
rm( map_data)
# centering and scaling covariate
tidx <- which(apply(X,2,var)==0)
if( length(tidx)>0){
  warning(paste("Some of the columns of X are constants, we removed" ,length(tidx), "columns"))
  X <- X[,-tidx]
}
X <- colScale(X)
# centering input
Y <- colScale(Y, scale=FALSE)
W <- DWT2(Y)
Y_f      <-  cbind( W$D,W$C)

if(verbose){
  print("Starting initialization")
}
if(verbose){
  print("Data transform")
}

#X <- matrix(X)
### Definition of some static parameters ---

indx_lst <-  gen_wavelet_indx(log2(length( outing_grid)))
#removing wc with variance 0 or below a certain level

lowc_wc <-   which_lowcount(Y_f,thresh_lowcount)
if(verbose){
  print( paste("Discarding ", length(lowc_wc), "wavelet coefficients out of ", ncol(Y_f)))
}
if(length(lowc_wc)> (ncol(Y_f )-3)){
  stop("almost all the wavelet coefficients are null/low variance, consider using univariate fine mapping")
}


if(quantile_trans)# important to do it after testing for lowcount
{
  W$C <- Quantile_transform(W$C )

  W$D <- apply( W$D,2,  Quantile_transform )
  Y_f      <-  cbind( W$D,W$C)
}

v1       <-  rep(1, dim(X)[1])### used in fit_lm to add a column of 1 in the design matrix
# Wavelet transform of the inputs


update_D <- W
if(verbose){
  print("Data transform done")
}
### Definition of some dynamic parameters ------
if(verbose){
  print("Initializing prior")
}
update_Y    <- cbind( W$D,W$C) #Using a column like phenotype, temporary matrix that will be regularly updated
temp        <- init_prior(Y              = update_Y,
                          X              = X,
                          prior          = prior ,
                          v1             = v1,
                          indx_lst       = indx_lst,
                          lowc_wc        = lowc_wc,
                          control_mixsqp = control_mixsqp,
                          nullweight     = nullweight,
                          gridmult       = gridmult,
                          max_SNP_EM     = max_SNP_EM,
                          max_step_EM    = max_step_EM,
                          cor_small      = cor_small)
G_prior     <- temp$G_prior
tt          <- temp$tt

#Recycled for the first step of the while loop
susiF.obj   <-  init_susiF_obj(L_max=L,
                               G_prior=G_prior,
                               Y=Y,
                               X=X,
                               L_start=L_start,
                               greedy=greedy,
                               backfit=backfit)

if(verbose){
  print("Initialization done")
}


#### start while ------


iter <- 1
init=TRUE
df=NULL
  for( l in 1:susiF.obj$L)
  {

    #print(susiF.obj$alpha[[l]])
    update_Y  <-  cal_partial_resid(
      susiF.obj = susiF.obj,
      l         = (l-1)  ,
      X         = X,
      D         = W$D,
      C         = W$C,
      indx_lst  = indx_lst
    )

    if(verbose){
      print(paste("Fitting effect ", l,", iter" ,  iter ))
    }

    if(init){#recycle operation used to fit the prior


      EM_out <- gen_EM_out (tpi_k= get_pi_G_prior(G_prior),
                            lBF  = log_BF  (G_prior,
                                            Bhat     = tt$Bhat,
                                            Shat     = tt$Shat,
                                            lowc_wc  = lowc_wc,
                                            indx_lst = indx_lst,
                                            df       = df
                            )
      )

      init <- FALSE
    }else{


      tt   <- cal_Bhat_Shat(update_Y,X,v1 ,
                            lowc_wc =lowc_wc ,

                            cor_small2     = cor_small2)

      tpi <-  get_pi(susiF.obj,l)
      G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)
      if(is.pois){
        tt$Shat <- update_Shat_pois(Shat=tt$Shat,
                                    indx_lst=indx_lst,
                                    lowc_wc=lowc_wc
        )
        #   print( tt$Shat)
      }

      EM_out  <- EM_pi(G_prior        = G_prior,
                       Bhat           = tt$Bhat,
                       Shat           = tt$Shat,
                       indx_lst       = indx_lst,
                       init_pi0_w     = init_pi0_w,
                       control_mixsqp = control_mixsqp,
                       lowc_wc        = lowc_wc,
                       nullweight     = nullweight,
                       max_SNP_EM     = max_SNP_EM,
                       max_step       = max_step_EM,
                       df             = df
      )
      #plot(EM_out$lBF/(sum( EM_out$lBFlBF)))
    }

    #print(h)
    # print(EM_out$lBF[1:10])
    susiF.obj <-  update_susiF_obj(susiF.obj   = susiF.obj ,
                                   l           = l,
                                   EM_pi       = EM_out,
                                   Bhat        = tt$Bhat,
                                   Shat        = tt$Shat,
                                   indx_lst    = indx_lst,
                                   lowc_wc     = lowc_wc,
                                   df          = df
    )


  }#end for l in 1:L  -----

temp <- wavethresh::wd(rep(0, susiF.obj$n_wac))


temp$D     <- susiF.obj$fitted_wc[[3]][410,-32]
temp$C[length(temp$C)]     <- susiF.obj$fitted_wc[[3]][410,32]
plot(wavethresh::wr(temp),type="l")




temp$D     <- susiF.obj$fitted_wc[[3]][388,-32]
temp$C[length(temp$C)]     <- susiF.obj$fitted_wc[[3]][388,32]
lines(wavethresh::wr(temp),type="l")
temp$D     <- susiF.obj$fitted_wc[[5]][388,-32]
temp$C[length(temp$C)]     <- susiF.obj$fitted_wc[[5]][388,32]
lines(wavethresh::wr(temp),type="l")

temp$D     <- susiF.obj$fitted_wc[[4]][388,-32]
temp$C[length(temp$C)]     <- susiF.obj$fitted_wc[[4]][388,32]
lines(wavethresh::wr(temp),type="l")

susiF.obj <- update_alpha_hist(susiF.obj)

susiF.obj <- update_cal_cs(susiF.obj,
                           cov_lev=cov_lev)



susiF.obj$cs



dummy.cs <-  which_dummy_cs(susiF.obj ,
                            min.purity = min.purity,
                            X=X)


dummy.cs


  ####Check greedy/backfit and stopping condition -----
  susiF.obj <- greedy_backfit (susiF.obj,
                               verbose    = verbose,
                               cov_lev    = cov_lev,
                               X          = X,
                               min.purity = min.purity
  )
  sigma2    <- estimate_residual_variance(susiF.obj,Y=Y_f,X)
  #print(sigma2)
  susiF.obj <- update_residual_variance(susiF.obj, sigma2 = sigma2 )
  susiF.obj <- test_stop_cond(susiF.obj = susiF.obj,
                              check     = check,
                              cal_obj   = cal_obj,
                              Y         = Y_f,
                              X         = X,
                              D         = W$D,
                              C         = W$C,
                              indx_lst  = indx_lst)
  #  print(susiF.obj$alpha)
  #print(susiF.obj$ELBO)
  check <- susiF.obj$check



  iter <- iter +1


 #end while
