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
#' @param maxiter,tol Iteration cap and relative tolerance on the bracket width.
#'
#' @return Numeric vector of roots.
#' @keywords internal
#' @export
bisection <- function(f, lower, upper, ...,
                      auto_adjust_interval = FALSE,
                      maxiter = 1000, tol = 1e-5) {

  fl <- f(lower, ...)
  fu <- f(upper, ...)

  bracketed <- function(fl, fu) {
    (fl <= 0 & fu >= 0) | (fl >= 0 & fu <= 0)
  }

  if (auto_adjust_interval) {
    expand <- 0L
    while (any(!bracketed(fl, fu)) && expand < 50L) {
      bad <- !bracketed(fl, fu)
      upper[bad] <- pmin(upper[bad] * 2, 1e10)
      fu <- f(upper, ...)
      expand <- expand + 1L
    }
  }

  if (any(is.na(fl)) || any(is.na(fu)) ||
      any(!bracketed(fl, fu))) {
    stop("bisection: the root is not bracketed", call. = FALSE)
  }

  root <- rep(NA_real_, length(lower))
  root[fl == 0] <- lower[fl == 0]
  root[fu == 0] <- upper[fu == 0]
  active <- is.na(root)

  for (i in seq_len(maxiter)) {
    if (!any(active)) break
    mid <- (lower + upper) / 2
    fm  <- f(mid, ...)
    if (any(is.na(fm[active]))) {
      stop("bisection: function returned NA inside the bracket", call. = FALSE)
    }

    width_scale <- pmax(abs(mid), sqrt(.Machine$double.eps))
    done <- active & (fm == 0 | (upper - lower) <= tol * width_scale)
    root[done] <- mid[done]
    active[done] <- FALSE

    same_as_lower <- active &
      ((fm >= 0 & fl >= 0) | (fm <= 0 & fl <= 0))
    opposite_lower <- active & !same_as_lower
    lower[same_as_lower] <- mid[same_as_lower]
    fl[same_as_lower] <- fm[same_as_lower]
    upper[opposite_lower] <- mid[opposite_lower]
  }

  if (any(active)) {
    warning("bisection did not converge within maxiter", call. = FALSE)
    root[active] <- (lower[active] + upper[active]) / 2
  }
  root
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
