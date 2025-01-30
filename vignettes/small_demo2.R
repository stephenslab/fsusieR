# TO DO: Try doubling the number of CpG sites to increase power for
# fsusie.
library(fsusieR)
library(ggplot2)
library(cowplot)
set.seed(1)

# Simulate the toy methylation data.
n <- 100 # Number of samples.
m <- 32 # Number of CpGs.
p <- 12 # Number of SNPs.

# Generate the MAFs.
maf <- 0.05 + 0.45*runif(p)

# Generate ids for the SNPs and CpGs.
snp_ids <- paste0("rs",sample(100:1000,p))
cpg_ids <- paste0("CpG",1:m)

# Simulate the genotypes.
X <- (runif(n*p) < maf) +
     (runif(n*p) < maf)
X <- matrix(X,n,p,byrow = TRUE)
colnames(X) <- snp_ids
storage.mode(X) <- "double"

# This is the effect matrix.
F <- matrix(0,p,m)
F[1,9:16] <- 2.5
F[3,9:16] <- (-1.75)
F[9,25:32] <- 2
rownames(F) <- snp_ids
colnames(F) <- cpg_ids

# Simulate the methylation levels at the CpGs.
E <- matrix(4 * rnorm(n*m),n,m)
Y <- X %*% F + E
Y <- Y - min(Y)

# Run fsusie.
fit <- susiF(Y,X,L = 3,filter_cs = FALSE,post_processing = "TI",
             prior = "mixture_normal")
print(fit$cs)
plot_susiF_pip(fit)
plot_susiF_effect(fit)

# Compute association statistics for all CpG-SNP pairs.
pvals <- matrix(0,m,p)
rownames(pvals) <- cpg_ids
colnames(pvals) <- snp_ids
for (i in 1:m) {
  for (j in 1:p) {
    dat <- data.frame(x = X[,j],y = Y[,i])
    out <- lm(y ~ x,dat)
    pvals[i,j] <- summary(out)$coefficients["x","Pr(>|t|)"]
  }
}
pdat <- data.frame(cpg  = rep(cpg_ids,times = p),
                   snp  = rep(snp_ids,each = m),
                   pval = as.vector(pvals),
                   stringsAsFactors = FALSE)
pdat <- transform(pdat,
                  snp = factor(snp,snp_ids),
                  cpg = factor(cpg,rev(cpg_ids)))
pdat <- transform(pdat,pval = p.adjust(pval,method = "bonferroni"))
pdat <- transform(pdat,pval = -log10(pval))
i <- which(pdat$pval < 1.3)
pdat[i,"pval"] <- NA
p1 <- ggplot(pdat,aes(x = snp,y = cpg,size = pval)) +
  geom_point() +
  theme_cowplot(font_size = 12) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
print(p1)
