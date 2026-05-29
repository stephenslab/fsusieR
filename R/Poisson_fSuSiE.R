# This file used to define an alternate Pois_fSuSiE0 implementation.
# It has been superseded by R/Pois_fSuSiE.R. The function name is kept
# as a thin alias for backward compatibility; it will be removed in a
# future release.
#
# To get the new behaviour explicitly, call Pois_fSuSiE().

#' @rdname Pois_fSuSiE
#' @export
Pois_fSuSiE0 <- function(...) {
  .Deprecated("Pois_fSuSiE",
              msg = "Pois_fSuSiE0() is deprecated; use Pois_fSuSiE() instead.")
  Pois_fSuSiE(...)
}
