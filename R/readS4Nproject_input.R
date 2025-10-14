#' readS4Nproject_input
#'
#' reads in total bioenergy (1st gen, 2nd gen and residues)
#' demand and co2 prices from IMAGE model for Sim4Nexus project
#'
#' @param subtype IMAGE input to be read in: co2prices or bioenergy
#' @return magpie object at country-level resolution
#'
#' @author Felicitas Beier
#' @seealso \code{\link[madrat]{readSource}}
#' @examples
#' \dontrun{
#' a <- readSource("S4Nproject_input", convert = "onlycorrect", aggregate = FALSE)
#' }
readS4Nproject_input <- function(subtype = "co2prices") { # nolint: object_name_linter.
  # IMAGE-iso mapping file
  regionmapping <- read.csv("image_isocodes.csv")
  regionmapping <- data.frame(iso = regionmapping$ISO.code, Region = regionmapping$IMAGE.region.name)
  regionmapping <- as.data.frame(as.matrix(regionmapping), stringsAsFactors = FALSE)

  # match IMAGE countries with all iso countries
  iso           <- getISOlist(type = "all")
  mapping       <- data.frame(iso = iso, stringsAsFactors = FALSE)
  regionmapping <- merge(mapping, regionmapping, "iso", all.x = TRUE)
  # missing matches:
  regionmapping$Region[regionmapping$iso == "SRB"] <- "Eastern Europe"
  regionmapping$Region[regionmapping$iso == "SSD"] <- "Eastern Africa"
  regionmapping$Region[regionmapping$iso == "PSE"] <- "Middle East"
  regionmapping$Region[regionmapping$iso == "SXM"] <- "Rest Centr America"
  regionmapping$Region[regionmapping$iso == "ALA"] <- "OECD Europe"
  regionmapping$Region[regionmapping$iso == "ATA"] <- "Oceania"
  regionmapping$Region[regionmapping$iso == "BES"] <- "OECD Europe"
  regionmapping$Region[regionmapping$iso == "BLM"] <- "OECD Europe"
  regionmapping$Region[regionmapping$iso == "CUW"] <- "Rest Centr America"
  regionmapping$Region[regionmapping$iso == "GGY"] <- "OECD Europe"
  regionmapping$Region[regionmapping$iso == "GRL"] <- "OECD Europe"
  regionmapping$Region[regionmapping$iso == "IMN"] <- "OECD Europe"
  regionmapping$Region[regionmapping$iso == "JEY"] <- "OECD Europe"
  regionmapping$Region[regionmapping$iso == "MAF"] <- "Rest Centr America"
  regionmapping$Region[regionmapping$iso == "MNE"] <- "Eastern Europe"

  # read in data
  if (subtype == "co2prices") {
    data <- read.csv("IMAGE_data_26112020_co2price.csv", stringsAsFactors = FALSE)
  } else if (subtype == "bioenergy") {
    data <- read.csv("IMAGE_data_26112020_bioenergy.csv", stringsAsFactors = FALSE)
  } else {
    stop("Specify subtype to be read in: co2prices or bioenergy (total bioenergy demand)")
  }

  # region name correction
  data$Region[data$Region == "Russia+"]         <- "Russia +"
  data$Region[data$Region == "China+"]          <- "China +"
  data$Region[data$Region == "Rest C. America"] <- "Rest Centr America"
  data$Region[data$Region == "Rest S. America"] <- "Rest South America"
  data$Region[data$Region == "Rest S. Asia"]    <- "Rest South Asia"
  data$Region[data$Region == "Rest S. Africa"]  <- "Rest Southern Africa"

  data      <- merge(data, regionmapping, "Region")
  data      <- data.frame(scenario = data$Scenario,
                          year = data$Year,
                          value = data$Value,
                          iso = data$iso,
                          stringsAsFactors = FALSE)
  data$year <- paste0("y", data$year)
  x         <- collapseNames(as.magpie(data))

  if (subtype == "bioenergy") {
    # Disaggregate regional bioenergy to country values using population weight
    # country population in SSP2
    pop     <- collapseNames(calcOutput("Population",
                                        scenario = "SSP2",
                                        extension2150 = "none",
                                        aggregate = FALSE,
                                        years = paste0("y", c(2005, seq(2010, 2100, by = 10))),
                                        round = 6))

    # regional population
    popReg <- toolAggregate(pop, rel = regionmapping, weight = NULL)
    popReg <- as.data.frame(popReg)
    popReg <- data.frame(Region = popReg$Region,
                         year = as.character(popReg$Year),
                         Value = popReg$Value,
                         stringsAsFactors = FALSE)
    popReg <- merge(popReg, regionmapping, "Region")
    popReg <- collapseNames(as.magpie(data.frame(iso = popReg$iso, year = popReg$year, Value = popReg$Value)))

    # bioenergy demand weighted by population
    x <- x * pop / popReg
  }

  return(x)
}
