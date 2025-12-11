
# check_lfsr_calibration_functional.R
library(mvf.susie.alpha)
set.seed(2025)

# -------------------------
# Parameters (tune as needed)
# -------------------------
N <- 100            # samples
P <- 10         # SNPs
Tlen <- 1024          # functional domain length (2^6)
nsim <- 50         # number of replicates (300 is a good compromise)
L_fit <- 3          # number of effects to fit
lfsr_thresholds <- c(0.01, 0.02, 0.05, 0.1)  # thresholds to evaluate
postproc <- "HMM"                # use HMM postprocessing
verbose_fit <- FALSE
lev_res=9

# bump specification (location and amplitude)
bump_center <- round(Tlen * 0.5)
bump_width <- round(Tlen * 0.12)    # half-width of bump
bump_amp <- 1.5                     # amplitude of bump

# utility: build bump vector length Tlen
build_bump <- function(Tlen, center, width, amp) {
  x <- seq_len(Tlen)
  w <- width
  bump <- amp * exp(-((x - center)^2) / (2 * (w/2)^2))
  return(bump)
}
true_bump <- fsusieR::simu_IBSS_ash_vanilla(lev_res = lev_res)$sim_func
true_bump[sample(1: (2^lev_res),
                 size=floor(.7* (2^lev_res)),
                 replace = FALSE
                 )]=0
# Option: you can set some domain points to exactly zero to test null calibration outside bump
is_nonzero_point <- which(true_bump != 0)

# storage
# For pointwise calibration we'll collect, across replicates, the lfsr vector when
# there is a causal SNP. We treat points outside the bump as "null" (true effect ~ 0).
pointwise_calls <- matrix(0L, nrow=nsim, ncol=Tlen)  # 1 if lfsr < thresh (we'll evaluate per threshold)
trait_min_lfsr <- numeric(nsim)  # min lfsr across domain per replicate
cs_contains_causal <- logical(nsim)

result=list()

# Optionally allow null-only replicates (no causal SNP) to assess behavior under pure null:
# Here we simulate the causal SNP present always, but the bump amplitude can be set to 0 to test null.
for (rep in 1:nsim) {
  # 1) Simulate genotype
  G <- matrix(sample(0:2, size=N*P, replace=TRUE), nrow=N, ncol=P)
  causal_snp <- sample(1:P, 1)

  # 2) Simulate functional phenotype Y_f (N x Tlen): Gaussian noise + genotype*bump at causal SNP
  Yf <- matrix(rnorm(N * Tlen, sd=1), nrow=N, ncol=Tlen)
  # add genotype-mediated bump
  for (i in 1:N) {
    Yf[i, ] <- Yf[i, ] + G[i, causal_snp] * true_bump
  }
  Y <- list(Y_f = list(Yf), Y_u = NULL)   # single functional trait

  # pos argument: domain positions (needed by multfsusie)
  pos <- list(pos1 = 1:Tlen)

  # 3) Fit multfsusie
  mfit <-   multfsusie(Y = Y, X = G , L = L_fit, post_processing = postproc, verbose = verbose_fit)

  if (is.null(mfit)) {
    # record NA for this replicate
    pointwise_calls[rep, ] <- NA_integer_
    trait_min_lfsr[rep] <- NA_real_
    cs_contains_causal[rep] <- NA
    next
  }

  # 4) Extract pointwise functional LFSR.
  #    mfit$lfsr is a list of effects; each mfit$lfsr[[e]]$est_lfsr_functional[[1]] is a vector length Tlen
  # Aggregate across effects by taking minimum (conservative): for each domain point take min LFSR across effects
  Llen <- length(mfit$lfsr)
  lfsr_pointwise_by_effect <- matrix(NA_real_, nrow=Llen, ncol=Tlen)
  for (e in 1:Llen) {
    lf_e <- mfit$lfsr[[e]]$est_lfsr_functional
    if (!is.null(lf_e) && length(lf_e) >= 1) {
      # lf_e[[1]] should be vector length Tlen
      vec <- lf_e[[1]]
      if (length(vec) == Tlen) {
        lfsr_pointwise_by_effect[e, ] <- vec
      } else {
        # unexpected shape - try to coerce
        lfsr_pointwise_by_effect[e, ] <- rep(NA_real_, Tlen)
      }
    } else {
      lfsr_pointwise_by_effect[e, ] <- rep(NA_real_, Tlen)
    }
  }
  # collapse to per-point LFSR (min across effects)
  per_point_lfsr <- apply(lfsr_pointwise_by_effect, 2, min, na.rm=TRUE)
  # If all NA -> mark replicate NA
  if (all(is.na(per_point_lfsr))) {
    pointwise_calls[rep, ] <- NA_integer_
    trait_min_lfsr[rep] <- NA_real_
    cs_contains_causal[rep] <- NA
    next
  }

  result[[rep]]= list( lfsr=per_point_lfsr,
                       is_nonzero_point=is_nonzero_point)
  print(result)
}
##############################################
# Updated Calibration Evaluation
##############################################

# Identify null and signal points
null_points <- which(true_bump == 0)
nonnull_points <- which(true_bump != 0)

# Extract all valid LFSR arrays
valid_idx <- which(!sapply(result, is.null))
cat("Successful replicates:", length(valid_idx), "out of", nsim, "\n")

# Build matrix of dimension nsim_valid × Tlen
lfsr_mat <- do.call(rbind, lapply(result[valid_idx], function(x) x$lfsr))

# -------------------------------------------------
# 1. POINTWISE CALIBRATION ACROSS MULTIPLE THRESHOLDS
# -------------------------------------------------

calib_table <- data.frame(
  threshold = lfsr_thresholds,
  empirical_fp = NA_real_,
  empirical_tp = NA_real_   # optional: gives power at signal points
)

for (i in seq_along(lfsr_thresholds)) {
  thr <- lfsr_thresholds[i]

  # FP rate: fraction of null points incorrectly called
  emp_fp <- mean(lfsr_mat[, null_points, drop=FALSE] < thr)

  # Power (TP): fraction of nonzero points correctly called
  emp_tp <- mean(lfsr_mat[, nonnull_points, drop=FALSE] < thr)

  calib_table$empirical_fp[i] <- emp_fp
  calib_table$empirical_tp[i] <- emp_tp
}

print(calib_table)

# -------------------------------------------------
# 2. TRAIT-LEVEL DECISION (min LFSR across domain)
# -------------------------------------------------

trait_min_lfsr <- apply(lfsr_mat, 1, min)

trait_level_call <- trait_min_lfsr < 0.05

cat(sprintf("Trait-level detection rate (min LFSR<0.05): %.3f\n",
            mean(trait_level_call)))

# -------------------------------------------------
# 3. PLOT POINTWISE CALIBRATION FOR THRESHOLD 0.05
# -------------------------------------------------

emp_fp_pointwise <- colMeans(lfsr_mat[, null_points, drop=FALSE] < 0.05)
emp_tp_pointwise <- colMeans(lfsr_mat[, nonnull_points, drop=FALSE] < 0.05)

emp_freq <- rep(NA, Tlen)
emp_freq[null_points]    <- emp_fp_pointwise
emp_freq[nonnull_points] <- emp_tp_pointwise

plot(1:Tlen, emp_freq, type='h', lwd=2,
     xlab="Domain index",
     ylab="freq(lfsr<0.05)",
     main="Pointwise empirical freq of LFSR<0.05")
abline(h=0.05, col="red", lty=2)
points(nonnull_points, emp_freq[nonnull_points], pch=19, col="blue")
#legend("topright", legend=c("null points", "signal points", "nominal 0.05"),
#       pch=c(1,19,NA), col=c("black","blue","red"), lty=c(NA,NA,2))

cat("Calibration evaluation complete.\n")
