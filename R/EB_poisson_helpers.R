# Helpers that were referenced from R/EB_poisson_mean_routines.R but
# never defined. Adding them so the package loads and the bisection
# fallback in vga_pois_solver() works.

#' @title Vectorised bisection for monotonic functions
#'
#' @description Element-wise bisection root-finder used as the fallback
#' inside [vga_pois_solver()]. Assumes `f(lower, ...) < 0` and
#' `f(upper, ...) > 0` (the Poisson-VGA variance equation `h_v` is
#' monotone increasing in `v`). When `auto_adjust_interval = TRUE`, the
#' `upper` boundary is doubled (capped at `1e10`) for any element whose
#' bracket is not yet valid.
#'
#' @param f Function returning a numeric vector of the same length as
#'   `lower`. The first argument is the value of `v`.
#' @param lower,upper Numeric vectors giving the per-element bracket.
#' @param ... Extra arguments forwarded to `f`.
#' @param auto_adjust_interval Logical. If `TRUE`, try to expand `upper`
#'   when the root is not initially bracketed. Default `FALSE`.
#' @param maxiter,tol Iteration cap and tolerance on |f(mid)|.
#'
#' @return Numeric vector of roots.
#' @keywords internal
#' @export
bisection <- function(f, lower, upper, ...,
                      auto_adjust_interval = FALSE,
                      maxiter = 1000, tol = 1e-5) {

  fl <- f(lower, ...)
  fu <- f(upper, ...)

  if (auto_adjust_interval) {
    expand <- 0L
    while (any(sign(fl) == sign(fu)) && expand < 50L) {
      bad <- sign(fl) == sign(fu)
      upper[bad] <- pmin(upper[bad] * 2, 1e10)
      fu <- f(upper, ...)
      expand <- expand + 1L
    }
  }

  for (i in seq_len(maxiter)) {
    mid <- (lower + upper) / 2
    fm  <- f(mid, ...)
    same_as_lower <- sign(fm) == sign(fl)
    lower <- ifelse(same_as_lower, mid, lower)
    upper <- ifelse(same_as_lower, upper, mid)
    fl    <- ifelse(same_as_lower, fm,  fl)
    if (max(abs(fm), na.rm = TRUE) < tol) break
    if (max(upper - lower, na.rm = TRUE) < tol) break
  }
  (lower + upper) / 2
}


#' @title VGA objective with log1p link (not yet implemented)
#'
#' @description
#' Placeholder for the Poisson VGA objective under the alternative
#' `link = "log1p"` mode of [pois_smooth_split()]. The Poisson likelihood
#' under that link involves `E_q[log(exp(mu) - 1)]` (or `softplus`),
#' which has no closed-form Gaussian expectation. Implementing this
#' properly requires either a quadrature scheme or a different
#' variational family; we have not done that here.
#'
#' Calling this function raises an informative error rather than
#' failing silently with `object not found`. Use `link = "log"` for now.
#'
#' @keywords internal
#' @export
pois_mean_GP_log1p_opt_obj <- function(theta, x, s, beta, sigma2, n) {
  stop("link = \"log1p\" is not implemented: ",
       "pois_mean_GP_log1p_opt_obj() has no closed-form VGA objective. ",
       "Use link = \"log\" or implement a quadrature-based version.",
       call. = FALSE)
}

#' @rdname pois_mean_GP_log1p_opt_obj
#' @keywords internal
#' @export
pois_mean_GP_log1p_opt_obj_gradient <- function(theta, x, s, beta, sigma2, n) {
  stop("link = \"log1p\" is not implemented: ",
       "pois_mean_GP_log1p_opt_obj_gradient() has no closed-form VGA gradient. ",
       "Use link = \"log\" or implement a quadrature-based version.",
       call. = FALSE)
}
