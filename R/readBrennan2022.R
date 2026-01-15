#' @title readBrennan2022
#' @description Reads Critical Connectivity Areas as defined in Brennan, A.,
#' Naidoo, R., Greenstreet, L., Mehrabi, Z., Ramankutty, N., & Kremen, C. (2022).
#' Functional connectivity of the world’s protected areas. Science, 376(6597),
#' 1101–1104. https://doi.org/10.1126/science.abl8974
#' Protected areas (2020) and Key Biodiversity Areas/Global Safet Net areas were masked at
#' a spatial resolution of 10 arc seconds before aggregating the data to 0.5°.
#'
#' @param subtype Defines whether land area covered by Critical Connectivity Areas has
#' been masked by other conservation priority data. If Key Biodiversity Areas have
#' only been masked the option is \code{"KBA_masked"}. With \code{"KBA_GSN_masked"}, land
#' area covered by the Global Safety Net (distinct species assemblages cluster) is also
#' masked. This is useful for complementary scenario building.
#'
#' @return Returns magpie objects with the land area covered by Critical Connectivity areas
#'  that is NOT already covered by Key Biodiversity Areas or the Global Safety Net
#'  (distinct species assemblages cluster) and was unprotected in 2020.
#' @author Patrick v. Jeetze
#'
#' @examples
#' \dontrun{
#' readSource("Brennan2022", convert = "onlycorrect")
#' }
#'
#' @importFrom terra rast classify cellSize aggregate terraOptions
#' @importFrom withr local_tempdir defer
#' @importFrom mstools toolGetMappingCoord2Country
#'

readBrennan2022 <- function(subtype = "KBA_GSN_masked") {

  # Set up 'terra' options
  terraOptions(tempdir = local_tempdir(tmpdir = getConfig("tmpfolder")), todisk = TRUE, memfrac = 0.5)
  defer(terraOptions(tempdir = tempdir()))

  # read data
  if (subtype == "KBA_masked") {
    # WDPA and China protected areas are masked
    unprotectedCCALand <- rast(paste0(
      "./cca_land_unprotected_0.5.tif"
    ))
  } else if (subtype == "KBA_GSN_masked") {
    # WDPA and China protected areas, as well as areas
    # under the "global safety net" are masked
    unprotectedCCALand <- rast(paste0(
      "./cca_land_unprotected_gsn_masked_0.5.tif"
    ))
  } else {
    stop("Please select an existing subtype")
  }

  # get spatial mapping
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  # transform raster to magpie object
  out <- mbind(
    as.magpie(extract(unprotectedCCALand[["crop"]], map[c("lon", "lat")])[, "crop"], spatial = 1),
    as.magpie(extract(unprotectedCCALand[["past"]], map[c("lon", "lat")])[, "past"], spatial = 1),
    as.magpie(extract(unprotectedCCALand[["forest"]], map[c("lon", "lat")])[, "forest"], spatial = 1),
    as.magpie(extract(unprotectedCCALand[["other"]], map[c("lon", "lat")])[, "other"], spatial = 1)
  )

  # set dimension names
  dimnames(out) <- list(
    "x.y.iso" = paste(map$coords, map$iso, sep = "."),
    "t" = NULL,
    "data" = paste("CCA", c("crop", "past", "forest", "other"), sep = ".")
  )

  return(out)
}
