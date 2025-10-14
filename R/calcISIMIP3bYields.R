#' @title calcISIMIP3bYields
#'
#' @description reads and cleans up ISIMIP3b crop yield data
#'
#' @param subtype subtype of yield based on readISIMIPoutputs, for crop yields
#' @param cells   magpie or lpjcell
#' @param smooth  smooth cells via spline
#'
#' @return magpie object in cellular resolution
#' @author David Meng-Chuen Chen, Edna Molina Bacca
#'
#' @examples
#' \dontrun{
#' calcOutput("ISIMIP3bYields", aggregate = FALSE)
#' }
#'
#' @importFrom mstools toolHoldConstant
#' @importFrom madrat toolSplitSubtype
#' @importFrom magclass dimOrder magpply dimSums getNames mbind time_interpolate
#' @importFrom mstools toolCoord2Isocell toolSmooth

calcISIMIP3bYields <- function(subtype = "yields:EPIC-IIASA:ukesm1-0-ll:ssp585:default:3b",
                               smooth = TRUE, cells = "lpjcell") {

  if (grepl("historical", subtype)) {
    stop("Can only read full future scenarios for now, with historical already added")
  }

  st <- toolSplitSubtype(subtype, list(dataset = "yields",
                                       model   = c("LPJmL", "EPIC-IIASA", "pDSSAT",
                                                   "CYGMA1p74", "PROMET", "CROVER",
                                                   "ISAM", "LDNDC", "PEPIC"),
                                       gcm     = c("gfdl-esm4", "ipsl-cm6a-lr",
                                                   "mpi-esm1-2-hr", "mri-esm2-0",
                                                   "ukesm1-0-ll"),
                                       scen    = c("historical", "ssp126", "ssp370", "ssp585"),
                                       co2     = c("default", "2015co2"),
                                       version = c("2a", "2b", "3a", "3b")))

  pastSubtype <- paste(st$dataset, st$model, st$gcm, "historical", st$co2, st$version, sep = ":")
  past <- readSource("ISIMIP", subtype = pastSubtype, convert = FALSE)
  scen <- readSource("ISIMIP", subtype = subtype, convert = FALSE)

  x <- mbind(past[, 2:length(getYears(past)), ], scen)

  # Interpolation of y2015 after mbind of past and scen (NA due to harvest year correction)
  nameClean <- function(x, subtype) {
    getNames(x, dim = 1)[getNames(x, dim = 1) == "mai"] <- "maiz"
    getNames(x, dim = 1)[getNames(x, dim = 1) == "soy"] <- "soybean"
    getNames(x, dim = 1)[getNames(x, dim = 1) == "ri1"] <- "ricea"
    getNames(x, dim = 1)[getNames(x, dim = 1) == "ri2"] <- "riceb"
    getNames(x, dim = 1)[getNames(x, dim = 1) == "swh"] <- "springwheat"
    getNames(x, dim = 1)[getNames(x, dim = 1) == "wwh"] <- "winterwheat"
    getNames(x, dim = 2)[getNames(x, dim = 2) == "ir"] <- "irrigated"
    getNames(x, dim = 2)[getNames(x, dim = 2) == "rf"] <- "rainfed"
    return(x)
  }

  plantDay    <- readSource("GGCMICropCalendar", subtype = "planting_day")
  plantDay    <- collapseNames(plantDay[, , c("ri1", "ri2", "wwh", "swh", "soy", "mai")][, , c("rf", "ir")])
  maturityDay <- readSource("GGCMICropCalendar", subtype = "maturity_day")
  maturityDay <- collapseNames(maturityDay[, , c("ri1", "ri2", "wwh", "swh", "soy", "mai")][, , c("rf", "ir")])

  diff <- maturityDay - plantDay
  diff <- nameClean(diff, subtype)

  for (n in getNames(x)) {
    cellsCorr <- where(diff[, , n] < 0)$true$regions
    x[cellsCorr, "y2015", n] <- setYears((x[cellsCorr, "y2016", n] +
                                            x[cellsCorr, "y2014", n]) / 2,
                                         "y2015")
  }


  # read in mask
  harvArea <- collapseNames(readSource("GGCMICropCalendar",
                                       subtype = "fraction_of_harvested_area",
                                       convert = FALSE))
  harvArea <- nameClean(harvArea)

  # for wheat take higher yielding variety  based on highest mean yield between 1981 and 2011
  if (st$model == "CYGMA1p74") { # CYGMA has no winter wheat
    getNames(x, dim = 1)[getNames(x, dim = 1) == "springwheat"] <- "tece"
  } else {
    # use mask to select between spring and winter wheat yields
    tece <- collapseNames(x[, , "springwheat"] * harvArea[, , "springwheat"] +
                            x[, , "winterwheat"] * harvArea[, , "winterwheat"])
    # tece mask does not cover all cells, only current harv area.
    # Fill in other areas with higher yielding variety, based on historical 30 year averages
    higherw <- magpply(x[, 1981:2011, "springwheat"],
                       FUN = mean, MARGIN = c(1, 3)) > magpply(x[, 1981:2011, "winterwheat"],
                                                               FUN = mean, MARGIN = c(1, 3))
    higherw <- time_interpolate(setYears(higherw, 1961),
                                interpolated_year = getYears(x),
                                integrate_interpolated_years = TRUE)
    higherw <- collapseNames(ifelse(higherw == 1, x[, , "springwheat"], x[, , "winterwheat"]))

    tece <- ifelse(tece == 0, higherw, tece)

    tece <- add_dimension(collapseNames(tece), dim = 3.1, nm = "tece")

    tece <- add_dimension(collapseNames(tece), dim = 3.1, nm = "tece")
    x <- x[, , c("springwheat", "winterwheat"), invert = TRUE]
    x <- mbind(x, tece)
  }

  if (st$model == "CROVER") {
    # CROVER doesn't have ri2 data
    getNames(x, dim = 1)[getNames(x, dim = 1) == "ricea"] <- "rice_pro"
  } else {
    # take weighted average of rice yields by crop area,
    # multiply by two in cells where there is both
    # (to get multicropped yield per year instead of yield per harvest)
    multiMask  <- collapseNames(harvArea[, , "riceb"])
    rice <- multiMask * collapseNames(x[, , "ricea"] * harvArea[, , "ricea"] + x[, , "riceb"] * harvArea[, , "riceb"])
    rice <- add_dimension(collapseNames(rice), dim = 3.1, nm = "rice_pro")




    x <- x[, , c("ricea", "riceb"), invert = TRUE]
    x <- mbind(x, rice)
  }

  if (smooth == TRUE) {
    # smooth with spline
    x <- toolSmooth(x)
  }
  # set very small yields from smoothing to 0
  x[x < 0.001] <- 0

  # weight for aggregation
  cropAreaWeight <- dimSums(calcOutput("Croparea", sectoral = "kcr", physical = TRUE, irrigation = FALSE,
                                       cellular = TRUE, aggregate = FALSE, years = "y1995",
                                       round = 6)[, , getNames(x, dim = 1)],
                            dim = 3)

  # reduce number of cells if 59199 cells required
  if (cells == "magpiecell") {
    x <- toolCoord2Isocell(x, cells = cells)
  }

  return(list(x            = x,
              weight       = cropAreaWeight,
              unit         = "t/ha",
              description  = "ISIMIP3b GGCMI yields for soy rice wheat maize",
              isocountries = FALSE))
}
