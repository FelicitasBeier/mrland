#' @title toolMagpieLandClasses
#' @description Harmonises cellular land area data in terms of "crop", "past", "forest" and
#'   "other" (dim 3.2) with the 2020 LanduseInitialisation and returns it in MAgPIE land classes:
#'   \enumerate{
#'   \item each set in dim 3.1 is scaled down where its total exceeds land area minus urban land,
#'   \item open ecosystems are shifted from "other" to "past" via
#'     \code{\link{toolCorrectOpenEcosystemMismatch}},
#'   \item land types are split into MAgPIE land classes via \code{\link{toolSplitLandClasses}}.
#'   }
#'
#' @param x magpie object on lpjcell resolution without years; dim 3.2 holds
#'   "crop", "past", "forest", "other"
#' @param nclasses "seven" or "nine", see \code{\link{toolSplitLandClasses}}
#'
#' @return magpie object with the same dim 3.1 as \code{x} and the MAgPIE land classes in dim 3.2
#'
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolSplitLandClasses}}
#' @export
toolMagpieLandClasses <- function(x, nclasses = "seven") {
  .landUse2020 <- function(n) {
    lu <- calcOutput("LanduseInitialisation", nclasses = n, cellular = TRUE,
                     input_magpie = TRUE, aggregate = FALSE, years = 2020)
    getYears(lu) <- NULL
    getCells(lu) <- getCells(x)
    return(lu)
  }
  landUse5 <- .landUse2020("five")
  landUse9 <- .landUse2020("nine")

  landNoUrban <- setNames(dimSums(landUse5, dim = 3) - landUse9[, , "urban"], NULL)
  getYears(landNoUrban) <- getYears(x)
  landMismatch <- toolConditionalReplace(landNoUrban / dimSums(x, dim = 3.2), c(">1", "is.na()"), 1)
  x <- x * landMismatch

  x <- toolCorrectOpenEcosystemMismatch(x, landUse5)
  return(toolSplitLandClasses(x, landUse9, nclasses))
}
