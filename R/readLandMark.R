#' @title readLandMark
#' @description Reads Indigenous Peoples’ and Local Community Lands from LandMark
#' "LandMark. 2025. Indicative Areas of Indigenous and Community Land Rights. Data
#' file from LandMark: The Global Platform of Indigenous and Community Lands.
#' Available at www.landmarkmap.org."
#'
#' Protected areas (WDPA+China) were masked at a spatial resolution of
#' 10 arcseconds before aggregating the data to 0.5°.
#'
#' @param subtype Whether \code{"delineated"} or \code{"indicative"} lands and
#' territories are returned.
#'
#' @return Returns magpie objects with the land area covered by the IPLC lands
#' that was unprotected in 2020.
#' @author Patrick v. Jeetze
#'
#' @examples
#' \dontrun{
#' readSource("LandMark", convert = "onlycorrect")
#' }
#'
#' @importFrom terra rast classify cellSize aggregate terraOptions
#' @importFrom withr local_tempdir defer
#' @importFrom mstools toolGetMappingCoord2Country
#'

readLandMark <- function(subtype = "delineated") {
  # Set up 'terra' options
  terraOptions(tempdir = local_tempdir(tmpdir = getConfig("tmpfolder")), todisk = TRUE, memfrac = 0.5)
  defer(terraOptions(tempdir = tempdir()))

  # read data
  # the data includes unprotected land not covered
  # by WDPA and protected areas in China
  if (subtype == "delineated") {
    unprotectedLandMark <- rast(paste0(
      "./LandMark_IPLC_land_unprotected_0.5.tif"
    ))
  } else if (subtype == "indicative") {
    unprotectedLandMark <- rast(paste0(
      "./LandMark_IPLC_land_indicative_unprotected_0.5.tif"
    ))
  } else {
    stop("Please select an existing subtype")
  }

  # get spatial mapping
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  # transform raster to magpie object
  out <- mbind(
    as.magpie(extract(unprotectedLandMark[["crop"]], map[c("lon", "lat")])[, "crop"], spatial = 1),
    as.magpie(extract(unprotectedLandMark[["past"]], map[c("lon", "lat")])[, "past"], spatial = 1),
    as.magpie(extract(unprotectedLandMark[["forest"]], map[c("lon", "lat")])[, "forest"], spatial = 1),
    as.magpie(extract(unprotectedLandMark[["other"]], map[c("lon", "lat")])[, "other"], spatial = 1)
  )

  # set dimension names
  dimnames(out) <- list(
    "x.y.iso" = paste(map$coords, map$iso, sep = "."),
    "t" = NULL,
    "data" = paste(paste0("LandMark_IPLC_", subtype), c("crop", "past", "forest", "other"), sep = ".")
  )

  return(out)
}
