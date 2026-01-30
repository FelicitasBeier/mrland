#' @title readDinerstein2020
#' @description Reads Global Safety Net data set published by Dinerstein, E., Joshi,
#' A. R., Vynne, C., Lee, A. T. L., Pharand-Deschênes, F., França, M., Fernando, S.,
#' Birch, T., Burkart, K., Asner, G. P., & Olson, D. (2020). A “Global Safety Net”
#' to reverse biodiversity loss and stabilize Earth’s climate. Science Advances, 6(36),
#' eabb2824. https://doi.org/10.1126/sciadv.abb2824
#'
#' Protected areas and Key Biodiversity Areas were masked at
#' a spatial resolution of 10 arc seconds before aggregating the data to 0.5°.
#'
#' @param subtype Defines which cluster (see Dinerstein et al. 2020) of the
#' Global Safety Net is returned.The different subtypes for land are:
#' \code{"GSN:distinct_species_assemblages"}, \code{"GSN:rare_phenomena"},
#' \code{"GSN:areas_of_intactness"}, \code{"GSN:climate_stabilisation_tier1"} and
#' \code{"GSN:climate_stabilisation_tier2"}.
#'
#' @return Returns magpie objects with the land area covered by the Global Safety Net that
#' is NOT already covered by Key Biodiversity Areas and was unprotected in 2020.
#' @author Patrick v. Jeetze
#'
#' @examples
#' \dontrun{
#' readSource("Dinerstein2020", convert = "onlycorrect")
#' }
#'
#' @importFrom terra rast classify cellSize aggregate terraOptions
#' @importFrom withr local_tempdir defer
#' @importFrom mstools toolGetMappingCoord2Country
#'

readDinerstein2020 <- function(subtype = "GSN:distinct_species_assemblages") {

  # extract subtype
  subtype <- unlist(strsplit(subtype, split = ":"))[2]

  # Set up 'terra' options
  terraOptions(tempdir = local_tempdir(tmpdir = getConfig("tmpfolder")), todisk = TRUE, memfrac = 0.5)
  defer(terraOptions(tempdir = tempdir()))


  # read data
  # the data includes unprotected land not covered
  # by WDPA and protected areas in China
  if (subtype == "distinct_species_assemblages") {
    unprotectedGSNLand <- rast(paste0(
      "./GSN_distinct_species_assemblages_unprotected_land_0.5.tif"
    ))
    id <- "DSA"
  } else if (subtype == "rare_phenomena") {
    unprotectedGSNLand <- rast(paste0(
      "./GSN_rare_phenomena_unprotected_land_0.5.tif"
    ))
    id <- "RarePhen"
  } else if (subtype == "areas_of_intactness") {
    unprotectedGSNLand <- rast(paste0(
      "./GSN_areas_of_intactness_unprotected_land_0.5.tif"
    ))
    id <- "AreaIntct"
  } else if (subtype == "climate_stabilisation_tier1") {
    unprotectedGSNLand <- rast(paste0(
      "./GSN_climate_stabilisation_tier1_unprotected_land_0.5.tif"
    ))
    id <- "ClimTier1"
  } else if (subtype == "climate_stabilisation_tier2") {
    unprotectedGSNLand <- rast(paste0(
      "./GSN_climate_stabilisation_tier2_unprotected_land_0.5.tif"
    ))
    id <- "ClimTier2"
  } else {
    stop("Please select an existing subtype")
  }

  # get spatial mapping
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  # transform raster to magpie object
  out <- mbind(
    as.magpie(extract(unprotectedGSNLand[["crop"]], map[c("lon", "lat")])[, "crop"], spatial = 1),
    as.magpie(extract(unprotectedGSNLand[["past"]], map[c("lon", "lat")])[, "past"], spatial = 1),
    as.magpie(extract(unprotectedGSNLand[["forest"]], map[c("lon", "lat")])[, "forest"], spatial = 1),
    as.magpie(extract(unprotectedGSNLand[["other"]], map[c("lon", "lat")])[, "other"], spatial = 1)
  )

  # set dimension names
  dimnames(out) <- list(
    "x.y.iso" = paste(map$coords, map$iso, sep = "."),
    "t" = NULL,
    "data" = paste(paste0("GSN_", id), c("crop", "past", "forest", "other"), sep = ".")
  )

  return(out)
}
