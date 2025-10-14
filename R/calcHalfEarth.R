#' @title calcHalfEarth
#' @description Function calculates land area in 'Half Earth' conservation priority area
#'
#' @param cells (deprecated) always lpjcell (67420 cells)
#' @param nclasses Options are either "seven" or "nine".
#' \itemize{
#' \item "seven" separates primary and secondary forest and includes "crop", "past", "forestry",
#' "primforest", "secdforest", "urban" and "other"
#' \item "nine" adds the separation of pasture and rangelands, as well as a differentiation of
#' primary and secondary non-forest vegetation and therefore returns "crop", "past", "range",
#' "forestry", "primforest", "secdforest", "urban", "primother" and "secdother"
#' }
#'
#' @return magpie object in cellular resolution with different protection options in conservation priority areas
#' @author Patrick v. Jeetze, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("HalfEarth", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset addLocation
#' @importFrom magclass collapseDim
#' @importFrom mstools toolCoord2Isocell
#'

calcHalfEarth <- function(cells = "lpjcell", nclasses = "seven") {
  # Land area (in Mha):
  iniLU <- calcOutput("LanduseInitialisation",
                      cellular = TRUE,
                      nclasses = nclasses, input_magpie = TRUE,
                      years = "y1995", aggregate = FALSE)
  landArea <- dimSums(iniLU, dim = 3)

  # Half Earth Protection Share
  halfEarthShr <- readSource("HalfEarth", convert = "onlycorrect")
  getNames(halfEarthShr) <- "HalfEarth"
  getSets(halfEarthShr)  <- c("x", "y", "iso", "year", "data")

  # Land area to be protected by 2050 (in Mha)
  x <- halfEarthShr * landArea

  # get WDPA baseline data (1995 - 2020)
  wdpaBase <- calcOutput("ProtectedAreaBaseline",
                         aggregate = FALSE,
                         nclasses = nclasses, magpie_input = TRUE)

  # make sure that current WDPA protected areas
  # are not part of conservation priority targets
  x <- x - setYears(dimSums(wdpaBase[, "y2020", ], dim = 3), NULL)
  x <- toolConditionalReplace(x, "<0", 0)

  urbanLand <- calcOutput("UrbanLandFuture",
                          subtype = "LUH3", aggregate = FALSE,
                          timestep = "5year")

  # Conservation potential after 2020
  consvPot <- landArea - dimSums(wdpaBase[, "y2020", ], dim = 3) -
    setCells(urbanLand[, "y2020", "SSP2"], getCells(wdpaBase))
  consvPot <- toolConditionalReplace(consvPot, "<0", 0)

  # Where conservation priority area is larger than conservation potential
  # replace by conservation potential
  for (c in getNames(x)) {
    tmp <- x[, , c]
    tmp[tmp > consvPot] <- consvPot[tmp > consvPot]
    x[, , c] <- tmp
  }

  if (nclasses == "seven") {
    # calulate share of respective natveg classes
    natvegShr <- iniLU[, , c("primforest", "secdforest", "other")] /
      dimSums(iniLU[, , c("primforest", "secdforest", "other")], dim = 3)
    natvegShr <- toolConditionalReplace(natvegShr, "is.na()", 0)

    # magpie object containing all land classes and their conservation land shares
    landShr <- mbind(
      new.magpie(getCells(natvegShr), NULL, c("crop", "past", "forestry"), fill = 0),
      natvegShr[, , c("primforest", "secdforest")],
      new.magpie(getCells(natvegShr), NULL, "urban", fill = 0),
      natvegShr[, , "other"]
    )
  } else if (nclasses == "nine") {
    # calulate share of respective natveg classes
    natvegShr <- iniLU[, , c("primforest", "secdforest", "primother", "secdother")] /
      dimSums(iniLU[, , c("primforest", "secdforest", "primother", "secdother")], dim = 3)
    natvegShr <- toolConditionalReplace(natvegShr, "is.na()", 0)

    # magpie object containing all land classes and their conservation land shares
    landShr <- mbind(
      new.magpie(getCells(natvegShr), NULL, c("crop", "past", "forestry"), fill = 0),
      natvegShr[, , c("primforest", "secdforest")],
      new.magpie(getCells(natvegShr), NULL, "urban", fill = 0),
      natvegShr[, , "primother"],
      natvegShr[, , "secdother"]
    )
  }

  # compute land area reserved for conservation
  x <- x * setCells(landShr, getCells(x))

  return(list(x = x,
              weight = NULL,
              unit = "Mha",
              description = "Half earth conservation priority area in each land type",
              isocountries = FALSE))
}
