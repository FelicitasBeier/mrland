#' Tool: spatial_header
#'
#' Given a regionmapping (mapping between ISO countries and regions) the
#' function calculates a 0.5 degree spatial header for 0.5 degree magclass
#' objects
#'
#' @param mapping Either a path to a mapping or an already read-in mapping as
#' data.frame.
#' @return A vector with 59199 elements
#' @author Jan Philipp Dietrich
#' @export
#' @seealso \code{\link[madrat]{regionscode}}
#' @examples
#' \dontrun{
#' spatial_header("regionmappingMAgPIE.csv")
#' }
#'
spatial_header <- function(mapping) { # nolint
  if (is.character(mapping)) {
    map <- read.csv(mapping, sep = ";")
  } else if (is.data.frame(mapping)) {
    map <- mapping
  } else {
    stop("Mapping is provided in an unsupported format. It should be either a character or a data.frame!")
  }

  .getColumn <- function(map, what) {
    col <- which(names(map) %in% what)
    if (length(col) == 0) stop("No fitting column found")
    else if (length(col) > 1) stop("Ambiguous column selection")
    return(as.character(map[[col]]))
  }

  regionscode <- regionscode(map)
  reg <- .getColumn(map, c("region", "RegionCode"))
  names(reg) <- .getColumn(map, c("country", "CountryCode"))
  iso <- toolGetMapping("CountryToCellMapping.rds", where = "mstools")$iso
  spatialHeader <- paste(reg[iso], seq_along(iso),
                         sep = ".")
  return(spatialHeader)
}
