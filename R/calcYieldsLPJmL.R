#' @title calcYieldsLPJmL
#'
#' @description This function extracts yields from LPJmL
#'              for all years
#'
#' @param lpjml         Defines LPJmL version for main crop inputs
#' @param climatetype   Switch between different climate scenarios
#' @param cells         if cellular is TRUE: "magpiecell" for 59199 cells
#'                      or "lpjcell" for 67420 cells
#'
#' @return magpie object in cellular resolution
#'
#' @author Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("YieldsLPJmL", aggregate = FALSE)
#' }
#'
#' @importFrom madrat toolGetMapping
#' @importFrom withr local_options

calcYieldsLPJmL <- function(lpjml = "ggcmi_phase3_nchecks_bft_e511ac58",
                            climatetype = "GSWP3-W5E5:historical", cells = "lpjcell") {

  # Extract argument information
  if (grepl("historical", climatetype)) {
    stage <- "smoothed"
  } else {
    stage <- "harmonized2020"
  }

  # Increase object size limit
  local_options(magclass_sizeLimit = 1e+12)

  lpj2mag     <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")
  cropsLPJmL  <- unique(lpj2mag$LPJmL5)
  irrigTypes  <- c("irrigated", "rainfed")
  yields      <- list()

  for (crop in cropsLPJmL) {
    subdata        <- as.vector(outer(crop, irrigTypes, paste, sep = "."))
    yields[[crop]] <- calcOutput("LPJmL_new", version = lpjml, climatetype = climatetype,
                                 subtype = "harvest", subdata = subdata, stage = stage, aggregate = FALSE)
  }
  yields  <- mbind(yields)

  # Select cells to be returned
  if (cells == "magpiecell") {
    yields <- toolCoord2Isocell(yields)
  }

  # Check for NAs
  if (any(is.na(yields))) {
    stop("produced NA yields")
  }

  return(list(x            = yields,
              weight       = NULL,
              unit         = "t per ha",
              description  = "Yields in tons per hectar for LPJmL crop types.",
              isocountries = FALSE))
}
