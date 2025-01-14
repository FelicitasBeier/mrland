#' @title calcTradeTariff
#' @description calculate tarde tariffs from GTAP dataset
#' @param gtap_version type of GTAP data version
#' \itemize{
#' \item \code{GTAP7}
#' \item \code{GTAP8}
#' }
#'
#' @param type_tariff which producer price should be used
#' \itemize{
#' \item \code{type_tariff}
#' }
#'
#'@param bilateral calculates whether tariffs should be bilateral
#'
#' @return Trade tariffs as an MAgPIE object
#' @author Xiaoxi Wang, David M Chen
#' @examples
#'     \dontrun{
#'     x <- calcTradeTariff("GTAP7")
#'     }
#' @importFrom magclass as.data.frame add_columns
#' @importFrom magpiesets findset
#' @importFrom GDPuc toolConvertGDP
#' @importFrom mstools toolCountryFillBilateral
#'

calcTradeTariff<- function(gtap_version = "GTAP9", type_tariff = "total", bilateral = FALSE) { #nolint

  stopifnot(gtap_version %in% c("GTAP7", "GTAP8", "GTAP9"))
  vom  <- calcOutput(type = "GTAPTrade", subtype = paste(gtap_version, "VOM", sep = "_"),
                     aggregate = FALSE)
  voa  <- calcOutput(type = "GTAPTrade", subtype = paste(gtap_version, "VOA", sep = "_"),
                     aggregate = FALSE)
  vxmd <- calcOutput(type = "GTAPTrade", subtype = paste(gtap_version, "VXMD", sep = "_"),
                     aggregate = FALSE, bilateral = bilateral)
  viws <- calcOutput(type = "GTAPTrade", subtype = paste(gtap_version, "VIWS", sep = "_"),
                     aggregate = FALSE, bilateral = bilateral)
  x <- vom[, , getNames(viws)] / voa[, , getNames(viws)] #output at dom mkt prices/prod payment(farm gate)
  fillMean <- function(x) {
    stopifnot(is.magpie(x))
    for (k in getNames(x)) {
      mean <- mean(x[, , k], na.rm = TRUE)
      x[, , k][is.na(x[, , k])] <- mean
    }
    return(x)
  }

  x <- fillMean(x)
  if (type_tariff == "export") {
    #export tax: positive values are taxes, negative values are export subsidies

    vxwd <- calcOutput(type = "GTAPTrade", subtype = paste(gtap_version, "VXWD", sep = "_"),
                       aggregate = FALSE, bilateral = bilateral)
    y <- (vxwd - vxmd) / vxmd #(exports World Price - export mkt price)/export mkt price
    y[is.nan(y)] <- NA
  }else if (type_tariff == "import") {
    #import tariff: positive values are tariffs

    vims <- calcOutput(type = "GTAPTrade", subtype = paste(gtap_version, "VIMS", sep = "_"),
                       aggregate = FALSE, bilateral = bilateral)
    y <- (vims - viws) / vxmd #(imports mkt Price - imports wld price)/export mkt price
    y[is.nan(y)] <- NA
    y[is.infinite(y)] <- NA
  } else if (type_tariff == "total") {
    vxwd <- calcOutput(type = "GTAPTrade", subtype = paste(gtap_version, "VXWD", sep = "_"),
                       aggregate = FALSE, bilateral = bilateral)
    vims <- calcOutput(type = "GTAPTrade", subtype = paste(gtap_version, "VIMS", sep = "_"),
                       aggregate = FALSE, bilateral = bilateral)
    y <- (vxwd - vxmd + vims - viws) / vxmd
    y[is.nan(y)] <- NA
    y[is.infinite(y)] <- NA
  } else {
    stop(paste("type_tariff", type_tariff, "is not supported!"))
  }

  if (bilateral) {
    x <- add_dimension(x, dim = 1.2, nm = getItems(x, dim = 1))
  }

  y <- fillMean(y) * x #multiply tariff ratio by exporting country's value ratio

  kTrade <- findset("k_trade")
  missing <- setdiff(kTrade, getNames(y))
  y <- add_columns(y, dim = 3.1, addnm = missing)
  y[, , missing] <- 0
  y <- y[, , setdiff(getNames(y), kTrade), invert = TRUE]

  if (gtap_version == "GTAP7") {
    py <- 2005
  } else if (gtap_version == "GTAP8") {
    py <- 2007
  } else if (gtap_version == "GTAP9") {
    py <- 2011
  }

  p <- collapseNames(calcOutput("PriceAgriculture",
                                datasource = "FAO", aggregate = FALSE))[, py, ]

  pGlo <- calcOutput("IniFoodPrice", aggregate = FALSE, products = "k_trade")

  missing2 <- setdiff(getNames(pGlo), getNames(p))
  p <- add_columns(p, dim = 3.1, addnm = missing2)

  for (i in missing2) {
    p[, , i] <- as.numeric(pGlo[, , i])
  }

  if (bilateral) {
    p <- add_dimension(p, dim = 1.2,   nm = getItems(p, dim = 1))
  }
  if (gtap_version == "GTAP9") {
    y <- GDPuc::toolConvertGDP(y, unit_in = "current US$MER",
                               unit_out = "constant 2017 US$MER",
                               replace_NAs = "no_conversion")

    y <- setYears(y[, 2011, ], NULL)
  }

  out <- setYears(y[, , kTrade] * p[, , kTrade], NULL) #multiply tariff by exporter producer price

  if (bilateral) {
    out <- toolCountryFillBilateral(out, 0)
  } else {
    out <- toolCountryFill(out, 0)
  }

  getSets(out)[1] <- "Regions"
  # set export trade tariffs of ruminant meant in Inida high to prevent exports; Japan too
  #exporting region is indicated by dim 1.1 with set name "Regions", importers are "REG"
  out[list(Regions = "IND"), , "livst_rum"]   <- 10^3
  out[list(Regions = "JPN"), , "livst_rum"]   <- 10^5
  out[list(Regions = "JPN"), , "livst_pig"]   <- 10^5
  out[list(Regions = "JPN"), , "livst_chick"] <- 10^5
  out[list(Regions = "JPN"), , "livst_egg"]   <- 10^5
  out[list(Regions = "JPN"), , "livst_milk"]  <- 10^5
  # set trade tariffs of alcohol in Japan high to prevent unrealistic exports
  out[list(Regions = "JPN"), , "alcohol"]  <- 10^5

  # use sugar tariffs as a surrogate for sugar crops
  out[, , "sugr_cane"] <- setNames(out[, , "sugar"], "sugr_cane")
  out[, , "sugr_beet"] <- setNames(out[, , "sugar"], "sugr_beet")


  # remove tariffs for primary products that are not produced in a certain region i.e. sugarcane in EUR,
  # and only assigned because they are part of a larger GTAP group (i.e. sugarcane as part of sugarcrops)

  prod <- calcOutput("Production", aggregate = FALSE)[, , "dm", drop = TRUE]
  prodL <- calcOutput("Production", products = "kli", aggregate = FALSE)[, , "dm", drop = TRUE]
  prod <- mbind(prod, prodL)
  rm(prodL)
  #take last decade
  prod <- dimSums(prod[, 2000:2013, ], dim = 2)

  sectorMapping <- toolGetMapping(type = "sectoral", name = "mappingGTAPMAgPIETrade.csv", where = "mrland")
  sectorMapping <- sectorMapping[which(sectorMapping$gtap != "xxx" & sectorMapping$magpie != "zzz"), ]

  citems <- intersect(getItems(prod, dim = 3), sectorMapping$magpie)

  #some items have 2 mappings, use the just one
  sectorMapping <- sectorMapping[-which(sectorMapping$gtap %in% c("ctl", "pcr", "rmk")), ]
  mag <- sectorMapping$magpie
  names(mag) <- sectorMapping$gtap

  #make production share
  prodShrs <- new.magpie(cells_and_regions = getItems(prod, dim = 1), years = NULL, names = citems, fill = 0)
  for (i in citems) {
    group <- sectorMapping[which(sectorMapping$gtap == sectorMapping[which(
                                                                           sectorMapping$magpie == i), "gtap"]),
    "magpie"]
    group <- group[which(group %in% citems)]
    prodShrs[, , i] <- prod[, , i] / dimSums(prod[, , group])
  }

  prodShrs[is.na(prodShrs)] <- 0
  #give threshold of 5%
  prodShrs <- ifelse(prodShrs < 0.05, 0, 1)
  c2items <- intersect(getItems(out, dim = 3), citems)

  #rename dim2 to multiply only by the exporters
  if (bilateral) {
    getItems(out, dim = 1.2) <- paste0(getItems(out, dim = 1.2), "2")
    out[, , c2items] <- prodShrs[, , c2items] * out[, , c2items]
    getItems(out, dim = 1.2) <- gsub("[0-9]+", "", getItems(out, dim = 1.2))
  } else {
    out[, , c2items] <-  prodShrs[, , c2items] * out[, , c2items]
  }

  if (gtap_version == "GTAP9") {
    yr <- 2011
  } else {
    yr <- NULL
  }

  if (bilateral) {
    weight <- setYears(toolCountryFillBilateral(vxmd[, yr, ] *
                                                  voa[, yr, getNames(viws)], 0), NULL)
  } else {
    weight <- setYears(toolCountryFill(vxmd[, yr, ] *
                                         voa[, yr, getNames(viws)], 0), NULL)
  }
  missing <- setdiff(kTrade, getNames(weight))
  weight <- add_columns(weight, dim = 3.1, addnm = missing)
  weight[, , missing] <- 1
  weight <- weight[, , setdiff(getNames(weight), kTrade), invert = TRUE] + 10^-10


  description <- paste0(type_tariff, "trade tariff")
  unit <- "US$2017/tDM"

  return(list(x = out,
              weight = weight,
              unit = unit,
              description = description))
}
