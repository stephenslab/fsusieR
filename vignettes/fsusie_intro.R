## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,comment = "#>",fig.width = 5,
                      fig.height = 3,fig.align = "center",
					  fig.cap = "&nbsp;",dpi = 120)

## ----load-packages------------------------------------------------------------
library(susieR)
library(fsusieR)
library(ggplot2)
library(cowplot)

## ----sim-effects, fig.height=2, fig.width=5-----------------------------------
set.seed(1)
f1 <- simu_IBSS_per_level(7)$sim_func
f2 <- simu_IBSS_per_level(7)$sim_func
par(mar = c(2,2,2,2))
plot(f1,type = "l",col = "royalblue",lwd = 2)
abline(a = 0,b = 0,lty = "dotted")
lines(f2,type = "l",col = "limegreen",lwd = 2)

## ----sim-trait----------------------------------------------------------------
data(N3finemapping)
rsnr  <- 0.5
pos1  <- 25
pos2  <- 75
X     <- N3finemapping$X[,1:100]
n     <- nrow(X)
nloc  <- length(f1)
noise <- matrix(rnorm(n*nloc,sd = var(f1)/rsnr),n,nloc)
Y     <- tcrossprod(X[,pos1],f1) + tcrossprod(X[,pos2],f2) + noise

## ----run-susiF, results="hide"------------------------------------------------
fit <- susiF(Y,X,L = 10)

## ----credible-sets------------------------------------------------------------
fit$cs

## ----pips---------------------------------------------------------------------
cs1 <- fit$cs[[1]]
cs2 <- fit$cs[[2]]
fit$pip[cs1]
fit$pip[cs2]

## ----cor-cs1------------------------------------------------------------------
cor(X[,cs1])

## ----plot-effects, fig.height=2.5, fig.width=5--------------------------------
plot_susiF_effect(fit)

## ----plot-pips, fig.height=1.75, fig.width=5----------------------------------
plot_susiF_pip(fit)

## ----plot-all, fig.height=3.5, fig.width=5------------------------------------
res <- plot(fit)
plot_grid(res$pip,res$effect,
          nrow = 2,ncol = 1,
		  rel_heights = c(3,5))

## ----sim-effects-uneven, fig.height=2, fig.width=5----------------------------
set.seed(2)
data(N3finemapping)
f1 <- simu_IBSS_per_level(9)$sim_func
f2 <- simu_IBSS_per_level(9)$sim_func
pos <- sort(sample(512,130))
par(mar = c(2,2,2,2))
plot(f2,type = "l",lwd = 1.25)
points(pos,f2[pos],col = "black",pch = 20,cex = 0.75)
abline(a = 0,b = 0,col = "black",lty = "dotted")

## ----sim-data-uneven----------------------------------------------------------
rsnr   <- 0.03
f1_obs <- f1[pos]
f2_obs <- f2[pos]
X      <- N3finemapping$X[,1:100]
n      <- nrow(X)
nloc   <- length(pos)
noise  <- matrix(rnorm(n*nloc,sd = var(f1)/rsnr),n,nloc)
Y      <- tcrossprod(X[,pos1],f1_obs) + tcrossprod(X[,pos2],f2_obs) + noise

## ----run-susiF-uneven, results="hide"-----------------------------------------
fit <- susiF(Y,X,L = 10,pos=pos)

## ----plot-effects-uneven, fig.height=2.5, fig.width=5-------------------------
plot_susiF_effect(fit)

## ----fsusie-priors, message=FALSE---------------------------------------------
fit_mnps <- susiF(Y,X,L = 10,prior = "mixture_normal_per_scale",verbose = FALSE)
fit_mn <- susiF(Y,X,L = 10,prior = "mixture_normal",verbose = FALSE)
fit_mnps$runtime
fit_mn$runtime

