#' @title correctOzoneYieldShock
#' @description correct Ozone Yield shock data
#' @return x corrected magpie object containing all ISO countries
#' @param x magpie object provided by the read function
#' @author Jake Tommey
#' @examples
#' \dontrun{
#' readSource("OzoneShock", convert="onlycorrect")
#' }
#' @importFrom madrat toolGetMapping


correctOzoneYieldShock <- function(x) {
  x <- toolConditionalReplace(x, conditions = "is.na()", replaceby = 0)
  isocountries <- unique(toolGetMapping("mapCoords2Country.rds", where = "mstools")$iso)
  missingImport <- c("QAT", "HTI", "IRL", "CHN", "IND",
    "ZAF", "MDG", "CRI", "AUS", "HTI",
    "MYS", "ZAF", "IDN", "VEN", "CHN"
  )
  names(missingImport) <- c("BHR", "DOM", "GBR", "HKG", "LKA",
    "LSO", "MUS", "NIC", "NZL", "PRI",
    "SGP", "SWZ", "TLS", "TTO", "TWN"
  )
  x <- toolCountryFill(x, fill = 0.0, countrylist = isocountries, verbosity = 0, map = missingImport)
  return(x)
}
