pkgload::load_all('.', quiet = TRUE)
results <- list()
for (seed in 1:3) {
  set.seed(seed + 100)
  X <- matrix(rnorm(48 * 9), 48, 9)
  X[, 4] <- .9 * X[, 2] + sqrt(.19) * X[, 4]
  Y <- outer(X[, 2], c(rep(0, 4), rep(2, 8), rep(0, 4))) +
    outer(X[, 7], sin(seq(0, 2*pi, length.out = 16))) +
    matrix(rnorm(48 * 16), 48, 16)
  X <- cbind(1, X[, 1:4], 0, X[, 5:9])
  colnames(X) <- paste0('variant', seq_len(ncol(X)))
  for (prior in c('mixture_normal', 'mixture_normal_per_scale')) {
    for (penalty in c(10, 100, 1000)) {
      for (L in c(1L, 4L)) {
        fit <- suppressWarnings(susiF(Y, X, L = L, L_start = 1,
          prior = prior, nullweight = penalty, standardize_Y = seed == 2,
          maxit = 12, cov_lev = c(.8, .95, .99)[seed],
          post_processing = 'none', verbose = FALSE, e = 0))
        null <- fsusieR:::.fs_null_effects(fit)
        a <- do.call(rbind, fit$alpha[!null])
        expected <- if (is.null(a)) rep(0, ncol(X)) else 1-apply(1-a,2,prod)
        stopifnot(max(abs(fit$pip - expected)) < 1e-12,
          all(is.finite(fit$lbf)), all(fit$pip[c(1,6)] == 0),
          all(abs(vapply(fit$alpha,sum,numeric(1))-1) < 1e-12),
          fit$n_cs == sum(lengths(fit$cs)>0), fit$niter <= 12)
        for (l in seq_along(fit$cs)) {
          j <- fit$cs[[l]]
          if (!length(j)) next
          stopifnot(sum(fit$alpha[[l]][j]) >= fit$cov_lev - 1e-12,
            all(fit$pip[j] >= fit$alpha[[l]][j] - 1e-12),
            identical(names(j),colnames(X)[j]))
        }
        results[[length(results)+1L]] <- data.frame(seed,prior,penalty,L,
          n_cs=fit$n_cs,niter=fit$niter,converged=fit$converged)
      }
    }
  }
}
write.csv(do.call(rbind,results),'tmp/package-review/simulation-checks.csv',row.names=FALSE)
cat(length(results),'simulated fits passed all output invariants.\n')

obj <- structure(list(L=2L,P=3L,alpha=list(c(.94,.06,0),c(.94,0,.06)),
                      cs=vector('list',2)),class='susiF')
obj <- update_cal_pip(update_cal_cs(obj,.95))
before <- ggplot2::ggplot(data.frame(snp=1:3,pip=obj$pip,CS=factor(c(2,1,2))),
    ggplot2::aes(snp,pip,color=CS)) + ggplot2::geom_point(size=3) +
  ggplot2::scale_color_manual(values=c('dodgerblue2','green4')) +
  cowplot::theme_cowplot() + ggplot2::labs(title='Before: shared SNP assigned only to CS2',x='SNP',y='PIP') +
  ggplot2::ylim(0,1.08)
after <- plot_susiF_pip(obj,point_size=3,title='After: both CS memberships are visible') +
  ggplot2::ylim(0,1.08)
ggplot2::ggsave('tmp/package-review/cs-overlap-before-after.png',
  cowplot::plot_grid(before,after,nrow=1),width=11,height=4,dpi=160,bg='white')
