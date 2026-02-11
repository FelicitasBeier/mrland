#' @title readLUH2UrbanFuture
#' @description read in gridded future urban land use datasets, from LUH2 Hurtt data
#' @return magpie object of gridded future urban land use in Mha, 2015-2100
#' @author David Chen, Patrick v. Jeetze
#' @seealso \code{\link[madrat]{readSource}}
#' @importFrom magclass as.magpie
#' @importFrom ncdf4 nc_open
#' @importFrom terra ext rast aggregate extract ext<-
#' @importFrom magclass as.magpie mbind
readLUH2UrbanFuture <- function() {
  files <- c("states_ssp1_19.nc", "states_ssp2_45.nc", "states_ssp3_70.nc",
             "states_ssp4_60.nc", "states_ssp5_85.nc")

  # select years
  startYear <- 2015 # min 850
  endYear <- 2100
  timesteps <- 1
  offset <- 2015 # year 2016=1, y2100=86

  map <- toolGetMappingCoord2Country(pretty = TRUE)

  # land area
  carea <- rast("staticData_quarterdeg.nc")[["carea"]]
  ext(carea) <- c(-180, 180, -90, 90)

  dataSel <- "urban" # can change if ever other future LUH2v2 land use types are needed
  timeSel <- seq(startYear - offset, endYear - offset, by = timesteps)

  yrs <- paste0("y", seq(startYear, endYear, timesteps))
  ssps <- paste0("SSP", 1:5)

  urbanMag <- NULL
  for (i in seq_along(files)) {
    tmpMag <- NULL
    for (t in seq_along(timeSel)) {
      shr <- rast(files[i], subds = dataSel)
      shr <- shr[[t]]
      x <- shr * carea
      x <- aggregate(x, fact = 2, fun = sum, na.rm = TRUE)
      x <- as.magpie(extract(x, map[c("lon", "lat")])[, 2], spatial = 1)
      dimnames(x) <- list("x.y.iso" = paste(map$coords, map$iso, sep = "."),
                          "t" = yrs[t],
                          "data" = ssps[i])
      tmpMag <- mbind(tmpMag, x)
    }
    urbanMag <- mbind(urbanMag, tmpMag)
  }

  # convert from km^2 to Mha
  urbanMag <- urbanMag / 10000

  out <- urbanMag

  return(out)
}
