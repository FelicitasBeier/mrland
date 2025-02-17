convertREMMAG <- function(x, subtype) {
  map <- toolGetMapping(type = "regional", name = "regionmappingMAgPIE.csv", where = "mappingfolder")

  if (subtype == "ghgprices") {
    y <- toolAggregate(x, map)
  } else if (subtype == "biodem") {
    pop <- calcOutput("Population", scenario = "SSP2", aggregate = FALSE)
    y <- toolAggregate(x, map, weight = pop[, 2010, ])
  } else {
    stop("Unknown subtype ", subtype)
  }
  return(y)
}
