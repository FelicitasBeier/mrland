#' @title readH08evapotranspiration
#' @description Read evapotranspiration data
#' @param subtype Switch between different inputs
#' @return MAgPIE objects with results on cellular level.
#' @author Marcos Alves
#' @seealso
#' \code{\link{readH08evapotranspiration}}
#' @examples
#' \dontrun{
#' readSource("H08evapotranspiration", subtype, convert = "onlycorrect")
#' }
#'
#' @importFrom terra rast
#' @importFrom raster brick
#' @importFrom magclass as.magpie
#' @export

readH08evapotranspiration <-
  function(subtype = "H08:mri-esm2-0:historical") {
    x <- toolSplitSubtype(subtype, list(water_model = NULL, climate_model = NULL, scenario = NULL))
    files <- Sys.glob("*.nc")
    fileName <- grep(x$scenario, files, value = TRUE)
    r <-  suppressWarnings(brick(fileName))

    y <- as.magpie(r, temporal = 1)
    y <- toolCoord2Isocell(y)

    if (x$scenario == "historical") {
      getItems(y, dim = 2) <- paste0("y", seq(1850, 2014, 1))
    } else {
      getItems(y, dim = 2) <- paste0("y", seq(2015, 2100, 1))
    }
    y <- setNames(y, x$scenario)

    return(y)
  }
