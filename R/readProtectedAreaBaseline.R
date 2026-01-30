#' @title readProtectedAreaBaseline
#' @description Reads spatial land cover information within protected areas.
#'              Land cover information for protected areas has been extracted from ESA CCI
#'              land use/land cover data (https://www.esa-landcover-cci.org/) and data from
#'              the WDPA data base (https://www.protectedplanet.net).
#' @return Returns magpie object with the protected area separated for each land type
#'         (cropland, pasture, forest, other land) per grid cell from 1995 to 2020.
#' @author Patrick v. Jeetze
#'
#' @examples
#' \dontrun{
#' readSource("ProtectedAreaBaseline", convert = "onlycorrect")
#' }
#' @importFrom terra aggregate cellSize classify extract rast segregate terraOptions tmpFiles
#' @importFrom mstools toolGetMappingCoord2Country
#' @importFrom withr local_tempdir defer
#'

readProtectedAreaBaseline <- function() {
  # Set up 'terra' options
  terraOptions(tempdir = local_tempdir(tmpdir = getConfig("tmpfolder")), todisk = TRUE, memfrac = 0.5)
  defer(terraOptions(tempdir = tempdir()))

  # Historic reference years
  refyears <- c(1995, 2000, 2005, 2010, 2015, 2020)

  # create output object
  out <- NULL

  for (yr in refyears) {
    message(paste("Beginning year", yr))

    # create output object for selected year
    # to iterate over IUCN protected area categories
    yrOut <- NULL

    for (cat in c("all", "I-II-III")) {
      # Read raster file: The file contains land cover information for legally
      # designated protected areas in the reference year, estimated based on
      # ESA CCI land use/land cover (LULC)
      if (cat == "all") {
        wdpaRast <- rast(paste0("./wdpa+china_esacci_land_", yr, "_0.5deg.tif"))
      } else {
        wdpaRast <- rast(paste0("./wdpa+china_", cat, "_esacci_land_", yr, "_0.5deg.tif"))
      }

      # get spatial mapping
      map <- toolGetMappingCoord2Country(pretty = TRUE)
      # transform raster to magpie object
      wdpaMag <- mbind(
        as.magpie(extract(wdpaRast[["crop"]], map[c("lon", "lat")])[, "crop"], spatial = 1),
        as.magpie(extract(wdpaRast[["past"]], map[c("lon", "lat")])[, "past"], spatial = 1),
        as.magpie(extract(wdpaRast[["forest"]], map[c("lon", "lat")])[, "forest"], spatial = 1),
        as.magpie(extract(wdpaRast[["other"]], map[c("lon", "lat")])[, "other"], spatial = 1)
      )
      # set dimension names
      if (cat == "all") {
        dimnames(wdpaMag) <- list(
          "x.y.iso" = paste(map$coords, map$iso, sep = "."),
          "t" = paste0("y", yr),
          "data" = c("WDPA.crop", "WDPA.past", "WDPA.forest", "WDPA.other")
        )
      } else {
        dimnames(wdpaMag) <- list(
          "x.y.iso" = paste(map$coords, map$iso, sep = "."),
          "t" = paste0("y", yr),
          "data" = paste0("WDPA_", cat, c(".crop", ".past", ".forest", ".other"))
        )
      }

      # bind to output
      yrOut <- mbind(yrOut, wdpaMag)
    }

    # bind to output
    out <- mbind(out, yrOut)

    message(paste("Finished year", yr))
  }

  return(out)
}
