#' @title readKeyBiodiversityAreas
#' @description Reads land area covered by for Key Biodiversity Areas
#' (https://www.keybiodiversityareas.org/) that was unprotected in 2020.
#' Protected areas were masked at a spatial resolution of 10 arc seconds
#' before aggregating the data to 0.5°.

#' @param subtype "unprotected" or "all"
#'
#' @return Returns magpie objects with the area covered by unprotected Key Biodiversity Areas per grid cell
#' @author Patrick v. Jeetze
#'
#' @examples
#' \dontrun{
#' readSource("KeyBiodiversityAreas", convert = "onlycorrect")
#' }
#'
#' @importFrom terra rast classify cellSize aggregate terraOptions
#' @importFrom withr local_tempdir defer
#' @importFrom mstools toolGetMappingCoord2Country
#'

readKeyBiodiversityAreas <- function(subtype = "unprotected") {

  # Set up 'terra' options
  terraOptions(tempdir = local_tempdir(tmpdir = getConfig("tmpfolder")), todisk = TRUE, memfrac = 0.5)
  defer(terraOptions(tempdir = tempdir()))

  # read data
  if (subtype == "unprotected") {
    # the data includes unprotected land not covered
    # by WDPA and protected areas in China
    kbaLand <- rast("./kba_land_unprotected_0.5.tif")
  } else if (subtype == "all") {
    kbaLand <- rast("./kba_land_0.5.tif")
  }

  # get spatial mapping
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  # transform raster to magpie object
  out <- mbind(
    as.magpie(extract(kbaLand[["crop"]], map[c("lon", "lat")])[, "crop"], spatial = 1),
    as.magpie(extract(kbaLand[["past"]], map[c("lon", "lat")])[, "past"], spatial = 1),
    as.magpie(extract(kbaLand[["forest"]], map[c("lon", "lat")])[, "forest"], spatial = 1),
    as.magpie(extract(kbaLand[["other"]], map[c("lon", "lat")])[, "other"], spatial = 1)
  )

  # set dimension names
  dimnames(out) <- list(
    "x.y.iso" = paste(map$coords, map$iso, sep = "."),
    "t" = NULL,
    "data" = paste("KBA", c("crop", "past", "forest", "other"), sep = ".")
  )

  return(out)
}
