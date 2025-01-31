# Run fsusie.
fit <- susiF(Y,X,L = 3,filter_cs = FALSE,post_processing = "TI",
             prior = "mixture_normal")
print(fit$cs)
plot_susiF_pip(fit)
plot_susiF_effect(fit)

stop()

# Compute association statistics for all CpG-SNP pairs, and create
# various plots of the CpG-SNP association statistics.
pdat <- data.frame(cpg  = rep(cpg_ids,times = p),
                   snp  = rep(snp_ids,each = m),
                   pval = as.vector(pvals),
                   stringsAsFactors = FALSE)
pdat <- transform(pdat,
                  snp = factor(snp,snp_ids),
                  cpg = factor(cpg,cpg_ids))
# pdat <- transform(pdat,pval = p.adjust(pval,method = "bonferroni"))
pdat <- transform(pdat,pval = -log10(pval))
p2 <- ggplot(pdat,aes(x = snp,y = pval)) +
  geom_point(shape = 21,size = 2,color = "white",fill = "royalblue") +
  labs(x = "",y = "-log10 p-value") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
p3 <- ggplot(pdat,aes(x = cpg,y = pval)) +
  geom_point(shape = 21,size = 2,color = "white",fill = "darkorange") +
  labs(x = "",y = "-log10 p-value") +
  theme_cowplot(font_size = 10) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
stop()
i <- which(pdat$pval < 1.3)
pdat[i,"pval"] <- NA
p1 <- ggplot(pdat,aes(x = snp,y = cpg,size = pval)) +
  geom_point() +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
