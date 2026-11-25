#' @title toolSplitLandClasses
#' @description Disaggregates land types "crop", "past", "forest", "other" (dim 3.2)
#'   into MAgPIE land classes. Forest (and other land for "nine") is split into primary
#'   and secondary by shares from \code{landUse9}; where \code{landUse9} has no such land,
#'   all of it goes to secondary. Forestry and urban are set to 0.
#'
#' @param x magpie object with 2 data dims; dim 3.2 holds "crop", "past", "forest", "other"
#' @param landUse9 nine-class LanduseInitialisation, cells and years aligned with \code{x}
#' @param nclasses "seven" or "nine"
#' \itemize{
#' \item "seven" returns "crop", "past", "forestry", "primforest", "secdforest", "urban" and "other"
#' \item "nine" returns "crop", "past" (0), "range", "forestry", "primforest", "secdforest",
#' "urban", "primother" and "secdother"
#' }
#'
#' @return magpie object with the same dim 3.1 as \code{x} and the MAgPIE land classes in dim 3.2
#'
#' @author Patrick v. Jeetze
#' @importFrom magclass dimSums getCells getItems getYears mbind new.magpie setItems setNames
#' @export
toolSplitLandClasses <- function(x, landUse9, nclasses = "seven") {
  if (!nclasses %in% c("seven", "nine")) {
    stop("Option specified for argument 'nclasses' does not exist.")
  }
  .split <- function(pool, prim, secd) {
    total <- setNames(dimSums(landUse9[, , c(prim, secd)], dim = 3) + 1e-10, NULL)
    primShr <- setNames(landUse9[, , prim] / total, NULL)
    secdShr <- setNames(landUse9[, , secd] / total, NULL)
    secdShr[secdShr == 0 & primShr == 0] <- 1
    mbind(setItems(primShr * x[, , pool], dim = 3.2, prim),
          setItems(secdShr * x[, , pool], dim = 3.2, secd))
  }
  # built from scratch rather than x * 0, so that they stay 0 where x is NA
  .zero <- function(name) {
    new.magpie(getCells(x), getYears(x), paste(getItems(x, dim = 3.1), name, sep = "."), fill = 0)
  }

  if (nclasses == "seven") {
    out <- mbind(x[, , c("crop", "past")], .zero("forestry"),
                 .split("forest", "primforest", "secdforest"),
                 .zero("urban"), x[, , "other"])
  } else {
    out <- mbind(x[, , "crop"], .zero("past"), setItems(x[, , "past"], dim = 3.2, "range"),
                 .zero("forestry"), .split("forest", "primforest", "secdforest"),
                 .zero("urban"), .split("other", "primother", "secdother"))
  }
  return(out)
}
