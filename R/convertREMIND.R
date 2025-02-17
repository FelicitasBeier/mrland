convertREMIND <- function(x, subtype) {

  # remove global dimension
  x <- x["GLO", , , invert = TRUE]

  map <- toolGetMapping(type = "regional", name = "regionmappingH12.csv", where = "madrat")

  if (grepl("intensive", subtype)) {
    # No weight for disaggregation because it's prices
    y <- toolAggregate(x, map)

  } else if (grepl("extensive", subtype)) {
    # Use population of 2010 as weight for disaggregation
    pop <- calcOutput("Population", scenario = "SSP2", aggregate = FALSE)
    y <- toolAggregate(x, map, weight = pop[, 2010, ])

  } else {
    stop("Unknown subtype ", subtype)
  }

  return(y)
}
