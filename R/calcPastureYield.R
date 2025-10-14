#' @title calcPastureYield
#' @description Provides pasture yields defined as
#'              ratio of grazed biomass to grazed area
#' @param range_pastr Boolean value indicating if the grass yields
#'                    should be split between rangelands and pastures.
#' @return Pasture yields and corresponding weights as a list of two MAgPIE objects
#' @author Isabelle Weindl, Marcos Alves
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link[mrcommons]{calcFAOmassbalance}},
#' \code{\link[madrat]{readSource}}
#' @examples
#' \dontrun{
#' calcOutput("PastureYield")
#' }
#' @importFrom stats quantile

calcPastureYield <- function(range_pastr = FALSE) { # nolint

  if (range_pastr) { # nolint

    # mapping country to coordinates
    mapping <- toolGetMappingCoord2Country()
    mapping$coordiso <- paste(mapping$coords, mapping$iso, sep = ".")

    # get years
    magYearsPast <- findset("past_til2020")[7:12]

    # read in country-level biomass
    biomass <- calcOutput("FAOmassbalance", aggregate = FALSE)[, , "production.dm"][, magYearsPast, "pasture"]
    biomass <- collapseNames(biomass)

    # read in cellular land data
    land         <- calcOutput("LanduseInitialisation", cellular = TRUE, cells = "lpjcell",
                               nclasses = "nine", aggregate = FALSE)

    # use only common years for the following calculations
    commonYears <- intersect(getYears(land), magYearsPast)
    biomass <- biomass[, commonYears, ]
    land <- land[, commonYears, ]

    grasslLand   <- land[, , c("past", "range")]
    grasslLand   <- setNames(grasslLand, c("pastr", "range"))
    grasslShares <- setNames(grasslLand[, , "pastr"] / dimSums(grasslLand, dim = 3), "pastr")
    grasslShares <- add_columns(grasslShares, addnm = "range", dim = 3.1)
    grasslShares[, , "range"] <- 1 - grasslShares[, , "pastr"]
    grasslShares[is.nan(grasslShares) | is.infinite(grasslShares)] <- 0

    # cellular livestock data
    livestock      <- setNames(readSource("GLW3"), "liv_numb")
    livstSplit     <- livestock * grasslShares
    livstSplit     <- collapseNames(livstSplit)

    # aggregate to country-level
    livstSplitCtry <- dimSums(livstSplit, dim = c("x", "y"))
    livstShareCtry <- livstSplitCtry[, , "pastr"] / dimSums(livstSplitCtry, dim = 3)
    livstShareCtry[is.nan(livstShareCtry) | is.infinite(livstShareCtry)] <- 0
    livstShareCtry <- add_columns(livstShareCtry, addnm = "range", dim = 3.1)
    livstShareCtry[, , "range"] <- 1 - livstShareCtry[, , "pastr"]

    # I am splitting biomass consumption assuming the share
    # between animals reared on rangelands and pastures correlates linearly
    # with the production of grass in pastures and rangelands in a country. That can be
    # derived by the fact that the feedbaskets assume the same feed ingredients shares
    # within a country.

    # countrylist
    cntry <- intersect(getItems(livstShareCtry, dim = 1),
                       getItems(biomass, dim = 1))

    biomassSplit <- biomass[cntry, , ] * livstShareCtry[cntry, , ]
    grasslLandCtry <- dimSums(grasslLand, dim = c("x", "y"))
    pstrYield <- biomassSplit / grasslLandCtry

    pstrYield[pstrYield > 100] <- 100
    pstrYield <- toolCountryFill(pstrYield)
    pstrYield[is.nan(pstrYield) | is.na(pstrYield)] <- 1

    grasslLandCtry <- toolCountryFill(grasslLandCtry)
    grasslLandCtry[is.na(grasslLandCtry)] <- 0

    return(list(x = pstrYield,
                weight = grasslLandCtry,
                isocountries = FALSE,
                unit = "ton DM per ha",
                description = "Pasture yields"))

  } else {

    magYearsPast <- findset("past_til2020")
    biomass   <- calcOutput("FAOmassbalance", aggregate = FALSE)[, , "production.dm"][, magYearsPast, "pasture"]
    biomass   <- collapseNames(biomass)
    pastLand  <- calcOutput("LanduseInitialisation", aggregate = FALSE)[, , "past"]

    commonYears <- intersect(getYears(pastLand), magYearsPast)
    pstrYield <- biomass[, commonYears, ] / pastLand[, commonYears, ]
    pstrYield[is.nan(pstrYield)] <- 1
    pstrYield[pstrYield > 100]   <- 100
    getNames(pstrYield) <- NULL

    weight <- pastLand[, commonYears, ]

    return(list(x = pstrYield,
                weight = weight,
                unit = "ton DM per ha",
                description = "Pasture yields"))
  }
}
