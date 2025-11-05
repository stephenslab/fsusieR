## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment = "#",collapse = TRUE,results = "hold",
                      fig.align = "center",dpi = 120)

## ----load-pkgs----------------------------------------------------------------
library(fsusieR)
library(reshape2)
library(ggplot2)
library(cowplot)

## ----set-seed-----------------------------------------------------------------
set.seed(1)

## ----sim-params---------------------------------------------------------------
n <- 100
m <- 32
p <- 12

## ----sim-mafs-----------------------------------------------------------------
maf <- 0.05 + 0.45*runif(p)

## ----sim-ids------------------------------------------------------------------
snpids <- paste0("SNP-",1:p)
cpgids <- paste0("CpG-",1:m)

## ----sim-geno-----------------------------------------------------------------
X <- (runif(n*p) < maf) +
     (runif(n*p) < maf)
X <- matrix(X,n,p,byrow = TRUE)
storage.mode(X) <- "double"
X[,4] <- X[,3] + 0.03*rnorm(n)
colnames(X) <- snpids

## ----sim-effect-matrix--------------------------------------------------------
F <- matrix(0,p,m)
F[1,9:16] <- 2.3
F[9,9:16] <- (-2.3)
F[3,25:32] <- 2
rownames(F) <- snpids
colnames(F) <- cpgids

## ----sim-Y--------------------------------------------------------------------
E <- matrix(3*rnorm(n*m),n,m)
Y <- X %*% F
Y <- Y + E
baseline <- min(Y)
Y <- Y - baseline

## ----plot-data, fig.height=2, fig.width=5-------------------------------------
pdat <- melt(Y)
x    <- X[,3]
pdat <- data.frame(cpg        = rep(1:m,each = n) +
                                runif(m*n,min = -0.2,max = 0.2),
                   meth_level = as.vector(Y),
				   geno       = factor(x))
pdat_lines <- data.frame(cpg = rep(1:m,times = 3),
                         geno = factor(rep(0:2,each = m)),
	 			 		 meth_level = rep(c(rep(0,8),rep(-0.75,8),rep(0,16)),
 						                  times = 3))
pdat_lines$meth_level <- pdat_lines$meth_level - baseline
rows <- which(with(pdat_lines,geno == 1 & cpg > 24))
pdat_lines[rows,"meth_level"] <- pdat_lines[rows,"meth_level"] + 2
rows <- which(with(pdat_lines,geno == 2 & cpg > 24))
pdat_lines[rows,"meth_level"] <- pdat_lines[rows,"meth_level"] + 4
p1 <- ggplot(pdat,aes(x = cpg,y = meth_level,color = geno)) +
  geom_point(shape = 20,size = 1.25) +
  scale_x_continuous(breaks = c(0,seq(4,32,4))) +
  scale_color_manual(values = c("darkblue","limegreen","darkorange")) +
  geom_line(data = pdat_lines,size = 0.75) +
  labs(x = "CpG",y = "methylation level",color = "genotype") +
  theme_cowplot(font_size = 11)
print(p1)

## ----plot-data-ggsave, echo=FALSE, eval=FALSE, message=FALSE------------------
# ggsave("demo_data.pdf",
#        p1 + scale_x_continuous(breaks = 1:m),
# 	   height = 2,width = 5)

## ----map-qtls-----------------------------------------------------------------
assoc <- matrix(0,m,p)
rownames(assoc) <- cpgids
colnames(assoc) <- snpids
for (i in 1:m) {
  for (j in 1:p) {
    dat <- data.frame(x = X[,j],y = Y[,i])
    fit <- lm(y ~ x,dat)
    assoc[i,j] <- summary(fit)$coefficients["x","Pr(>|t|)"]
  }
}

## ----gwas-by-snp, fig.height=2.5, fig.width=3---------------------------------
pdat <- data.frame(cpg    = rep(1:m,times = p),
                   snp    = rep(1:p,each = m),
				   effect = as.vector(t(F != 0)),
                   pval   = as.vector(assoc))
pdat <- transform(pdat,pval = -log10(pval))
threshold <- -log10(0.05/(m*p))
p2 <- ggplot(pdat,aes(x = snp,y = pval,shape = effect,color = effect)) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = threshold,color = "black",linetype = "dotted") +
  scale_x_continuous(breaks = 1:p) +
  scale_shape_manual(values = c(1,4)) +
  scale_color_manual(values = c("darkblue","darkorange")) +
  labs(x = "SNP",y = "-log10 p-value") +
  theme_cowplot(font_size = 11)
print(p2)

## ----gwas-by-snp-ggsave, echo=FALSE, eval=FALSE-------------------------------
# ggsave("demo_gwas_by_snp.pdf",p2,height = 3,width = 3)

## ----gwas-by-cpg, fig.height=2.5, fig.width=3.5-------------------------------
p3 <- ggplot(pdat,aes(x = cpg,y = pval,shape = effect,color = effect)) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = threshold,linetype = "dotted") +
  scale_x_continuous(breaks = c(0,seq(4,32,4))) +
  scale_shape_manual(values = c(1,4)) +
  scale_color_manual(values = c("darkblue","darkorange")) +
  labs(x = "CpG",y = "-log10 p-value") +
  theme_cowplot(font_size = 11)
print(p3)

## ----gwas-by-cpg-ggsave, echo=FALSE, eval=FALSE, message=FALSE----------------
# ggsave("demo_gwas_by_cpg.pdf",
#        p3 + scale_x_continuous(breaks = seq(1,32)),
# 	   height = 3,width = 3.75)

## ----run-fsusie, message=FALSE, results="hide"--------------------------------
fit <- susiF(Y,X,L = 3,filter_cs = FALSE,prior = "mixture_normal",
             post_processing = "HMM")

## ----fsusie-cs----------------------------------------------------------------
fit$cs

## ----fsusie-pip-plot, fig.height=1.5, fig.width=3-----------------------------
cs_colors <- c("dodgerblue","darkorange","red")
pdat <- data.frame(SNP = 1:p,
                   PIP = fit$pip,
				   CS  = rep("none",p))
pdat$CS[fit$cs[[1]]] <- "L1"
pdat$CS[fit$cs[[2]]] <- "L2"
pdat$CS[fit$cs[[3]]] <- "L3"
pdat <- transform(pdat,CS = factor(CS))
p4 <- ggplot(pdat,aes(x = SNP,y = PIP,fill = CS)) +
  geom_point(shape = 21,size = 2,color = "white") +
  scale_x_continuous(breaks = 1:p) +
  scale_fill_manual(values = c(cs_colors,"gray")) +
  theme_cowplot(font_size = 10)
print(p4)

## ----fsusie-pip-plot-ggsave, echo=FALSE, eval=FALSE---------------------------
# ggsave("demo_pip_plot.pdf",p4,height = 1.5,width = 3)

## ----fsusie-affected-cpgs, fig.height=2, fig.width=6, eval=FALSE--------------
# i <- sapply(fit$cs,function (x) intersect(x,c(1,3,9)))
# pdat <- data.frame(CS   = rep(c("L1","L2","L3"),each = m),
#                    CpG  = rep(1:m,times = 3),
#                    lfsr = unlist(fit$lfsr_func),
# 				   affected = as.vector(t(F[i,]) != 0),
# 				   stringsAsFactors = FALSE)
# pdat <- transform(pdat,
#                   CS   = factor(CS,c("L3","L2","L1")),
# 				  lfsr = -log10(lfsr))
# ggplot(pdat,aes(x = CpG,y = CS,size = lfsr,color = affected)) +
#   geom_point(shape = 1) +
#   scale_x_continuous(breaks = c(0,seq(4,32,4))) +
#   scale_color_manual(values = c("darkblue","darkorange")) +
#   scale_size(range = c(0.5,10),breaks = c(1.3,5,10)) +
#   labs(size = "-log10(lfsr)") +
#   theme_cowplot(font_size = 10)

## ----run-fsusie-TI, message=FALSE, results="hide"-----------------------------
fit_TI <- susiF(Y,X,L = 3,filter_cs = FALSE,prior = "mixture_normal",
                post_processing = "TI")

## ----fsusie-cs-TI-------------------------------------------------------------
range(fit$pip - fit_TI$pip)
fit_TI$cs

## ----fsusie-affected-cpgs-TI, fig.height=4.5, fig.width=3.5-------------------
effect_plot <- function (i) {
  pdat <- data.frame(cpg      = 1:m,
                     estimate = fit_TI$fitted_func[[i]],
					 lower    = fit_TI$cred_band[[i]]["low",],
					 upper    = fit_TI$cred_band[[i]]["up",])
  rows <- with(pdat,which(lower > 0 | upper < 0))
  pdat2 <- data.frame(cpg = rows,estimate = 0)
  return(ggplot(pdat,aes(x = cpg,y = estimate,ymin = lower,ymax = upper)) +
         geom_point(color = cs_colors[i],size = 1) +
         geom_errorbar(color = cs_colors[i],linewidth = 0.4) +
         geom_hline(yintercept = 0,color = "black",linewidth = 0.4) +
		 geom_point(data = pdat2,mapping = aes(x = cpg,y = estimate),
		            shape = 20,color = "black",size = 1.5,
					inherit.aes = FALSE) +
         scale_x_continuous(breaks = c(0,seq(4,32,4))) +
         labs(x = "CpG",y = "change",title = paste0("CS",i)) +
		 theme_cowplot(font_size = 11))
}
plot_grid(effect_plot(1),
          effect_plot(2),
   	      effect_plot(3),
		  nrow = 3,ncol = 1)

## ----fsusie-affected-cpgs-TI-ggsave, echo=FALSE, message=FALSE, eval=FALSE----
# ggsave("demo_affected_cpgs.pdf",
#        plot_grid(effect_plot(1) + scale_x_continuous(breaks = 1:m),
#                  effect_plot(2) + scale_x_continuous(breaks = 1:m),
#    	             effect_plot(3) + scale_x_continuous(breaks = 1:m),
# 		  nrow = 3,ncol = 1),
# 	   height = 4.5,width = 4)

## ----session-info-------------------------------------------------------------
sessionInfo()

