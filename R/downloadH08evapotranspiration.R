#' @title downloadH08evapotranspiration
#' @description Download water models evapotranspiration data
#' @param subtype Switch between different inputs
#' @author  Marcos Alves
#' @examples
#'
#' \dontrun{readSource("H08evapotranspiration",  convert="onlycorrect")}

downloadH08vapotranspiration <- function(subtype = "H08:mri-esm2-0:historical") {

  x        <- toolSplitSubtype(subtype,
                               list(waterModel = NULL,
                                    climateModel = NULL,
                                    scenario = NULL))
  path     <- "/p/projects/rd3mod/inputdata/sources/Evapotranspiration_raw" #nolint

  listFiles  <- list.files(path)
  file        <- grep("_evap", listFiles, value = TRUE)
  if (x$scenario == "historical") {
    file  <- grep("historical_histsoc", file, value = TRUE)
  } else {
    file  <- grep(paste0(x$scenario, "_2015soc"), file, value = TRUE)
  }

  filePath   <- file.path(path, file)

  if (file.exists(filePath)) {
    file.copy(filePath, file)
  } else {
    stop(paste("Data is not available so far:", filePath))
  }

  # Compose meta data
  return(list(url           = paste0(filePath),
              doi           = NULL,
              title         = x$scenario,
              author        = NULL,
              version       = x$waterModel,
              release_date  = NULL,
              description   = NULL,
              license       = NULL,
              reference     = NULL,
              unit          = "kg m-2 s-1")
  )
}
