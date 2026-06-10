#' @title toolCorrectOpenEcosystemMismatch
#' @description Corrects mismatches in the classification of open ecosystems
#'   into pasture and other between land-use initialisation (LUH) and ESA CCI.
#'   Excess "other" land relative to LUH is shifted to "past", but only up to
#'   the available LUH pasture surplus, to avoid spurious pasture expansion.
#'
#' @param x magpie object with "other" and "past" land-use classes.
#'   Cells and years already aligned with \code{luIni}
#' @param luIni magpie object from \code{calcOutput("LanduseInitialisation")} with aggregated
#'   "past" and "other" classes. Cells and years already aligned with \code{x}
#'
#' @return \code{x} with the other/pasture mismatch corrected
#'
#' @author Patrick v. Jeetze
#' @importFrom magclass setItems setNames
#' @export
toolCorrectOpenEcosystemMismatch <- function(x, luIni) {
  # calculate mismatch that is shifted to pasture
  otherMismatch <- x[, , "other"] - luIni[, , "other"]
  otherMismatch <- toolConditionalReplace(otherMismatch, c("<0", "is.na()"), 0)
  # but shift cannot be bigger than current LUH pasture
  # to avoid additional pasture expansion due to conservation measure
  pastMaxShift <- setNames(luIni[, , "past"], NULL) - x[, , "past"]
  pastMaxShift <- toolConditionalReplace(pastMaxShift, c("<0", "is.na()"), 0)
  otherMismatch <- pmin(otherMismatch, setItems(pastMaxShift, dim = 3.2, "other"))
  # subtract other land mismatch
  x[, , "other"] <- x[, , "other"] - otherMismatch
  # add to ESA CCI "pasture & rangeland" class instead
  x[, , "past"] <- x[, , "past"] + setItems(otherMismatch, dim = 3.2, "past")
  return(x)
}
