#' @title readNoon2022
#' @description Reads irrecoverble carbon data set published by Noon, M. L., Goldstein, A.,
#' Ledezma, J. C., Roehrdanz, P. R., Cook-Patton, S. C., Spawn-Lee, S. A., Wright, T. M.,
#' Gonzalez-Roglich, M., Hole, D. G., Rockström, J., & Turner, W. R. (2022). Mapping the
#' irrecoverable carbon in Earth’s ecosystems. Nature Sustainability, 5(1),
#' Article 1. https://doi.org/10.1038/s41893-021-00803-6
#' Protected areas were masked at a spatial resolution of 10 arc seconds
#' before aggregating the data to 0.5°.

#' @param subtype Defines whether carbon data or land area and related subtypes should be
#' returned (see options below).
#' Carbon or land subtypes need to be specified via ":"
#' The different subtypes for land are: "IrrC_30pc", "IrrC_40pc", "IrrC_50pc",
#' "IrrC_60pc", "IrrC_70pc", "IrrC_80pc", "IrrC_90pc", "IrrC_100pc"
#' which corresponds to the land area that was unprotected in 2020 and is covered by the respective
#' percentile of all irrecoverable carbon. \code{IrrC_50pc} e.g. returns all unprotected land that
#' contains the top 50\,\% of global irrecoverable carbon.

#' @return Returns magpie objects with the area of unprotected irrecoverable carbon land per grid cell
#' @author Patrick v. Jeetze

#' @examples
#' \dontrun{
#' readSource("Noon2022", convert = "onlycorrect")
#' }
#'
#' @importFrom terra rast classify cellSize aggregate terraOptions
#' @importFrom withr local_tempdir defer
#' @importFrom mstools toolGetMappingCoord2Country
#'

readNoon2022 <- function(subtype = "land:IrrC_50pc") {

  # extract subtype
  datatype <- unlist(strsplit(subtype, split = ":"))[1]
  subtype <- unlist(strsplit(subtype, split = ":"))[2]

  # Set up 'terra' options
  terraOptions(tempdir = local_tempdir(tmpdir = getConfig("tmpfolder")), todisk = TRUE, memfrac = 0.5)
  defer(terraOptions(tempdir = tempdir()))

  if (datatype == "land") {

    # read data
    # the data includes unprotected land not covered
    # by WDPA and protected areas in China
    if (subtype == "IrrC_30pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_30pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_40pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_40pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_50pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_50pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_60pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_60pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_70pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_70pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_75pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_75pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_80pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_80pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_90pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_90pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_95pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_95pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_99pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_99pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_100pc") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/Irrecoverable_C_100pc_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_75pc_30by30") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/IrrC_75pc_30by30_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_95pc_30by30") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/IrrC_95pc_30by30_unprotected_land_area_0.5deg.tif"
      ))
    } else if (subtype == "IrrC_99pc_30by30") {
      unprotectedICLand <- rast(paste0(
        "./unprotected_irrecoverable_C_land",
        "/IrrC_99pc_30by30_unprotected_land_area_0.5deg.tif"
      ))
    } else {
      stop("Please select an existing subtype")
    }

    # get spatial mapping
    map <- toolGetMappingCoord2Country(pretty = TRUE)
    # transform raster to magpie object
    out <- mbind(
      as.magpie(extract(unprotectedICLand[["crop"]], map[c("lon", "lat")])[, "crop"], spatial = 1),
      as.magpie(extract(unprotectedICLand[["past"]], map[c("lon", "lat")])[, "past"], spatial = 1),
      as.magpie(extract(unprotectedICLand[["forest"]], map[c("lon", "lat")])[, "forest"], spatial = 1),
      as.magpie(extract(unprotectedICLand[["other"]], map[c("lon", "lat")])[, "other"], spatial = 1)
    )

    # set dimension names
    dimnames(out) <- list(
      "x.y.iso" = paste(map$coords, map$iso, sep = "."),
      "t" = NULL,
      "data" = paste(subtype, c("crop", "past", "forest", "other"), sep = ".")
    )

    return(out)
  }
}
