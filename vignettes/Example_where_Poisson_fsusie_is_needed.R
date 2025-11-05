## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.height = 3,
  fig.align = "center",
  fig.cap = "&nbsp;",
  dpi = 175
  )

## ----setup, include=FALSE-----------------------------------------------------
rm(list = ls()) 
library(fsusieR)
library(susieR)
library(ebnm)
load("C:/Document/Serieux/Travail/Package/susiF.alpha/data/another_example.RData")

## -----------------------------------------------------------------------------
true_pos=example$true_pos
true_pos

## -----------------------------------------------------------------------------
plot(example$Y[1,], main="RNAseq count for individual 1 over the region of interest")

## -----------------------------------------------------------------------------
 
plot(example$lf , type="l", main = "effect of the causal SNP") 
X=example$X
Y=example$Y
 

## -----------------------------------------------------------------------------

res0 = susiF(X=X, Y=log1p(Y),L=3)
res1= susiF(X=X,Y=HFT(Y),L=3)

## -----------------------------------------------------------------------------
res0$cs
res1$cs
true_pos

## -----------------------------------------------------------------------------

res_poisF = Pois_fSuSiE(Y=Y,X=X ,L=3 , max.iter=2)


## -----------------------------------------------------------------------------
res_poisF$susiF.obj$cs

## -----------------------------------------------------------------------------
plot(res_poisF$susiF.obj$fitted_func[[1]])
lines(example$lf )

