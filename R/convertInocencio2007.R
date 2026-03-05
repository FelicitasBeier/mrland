#' @title convertInocencio2007
#'
#' @description Convert irrigation cost data by Inocencio et al. (2007)
#'              from regional to ISO country resolution
#'
#' @param x MAgPIE object containing irrigation cost data at region level resolution
#'
#' @return MAgPIE object on ISO country level
#'
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ a <- convertInocencio2005(x) }
#'
convertInocencio2007 <- function(x) {

  # list of iso countries
  isoCountry <- toolGetMapping("iso_country.csv", where = "mrland")

  # mapping of iso countries to regions
  iso2region <- data.frame(
    iso = c(
      "AFG","AGO","AIA","ALA","ALB","AND","ARE","ARG","ARM","ASM","ATA","ATF","ATG","AUS","AUT","AZE",
      "BDI","BEL","BEN","BES","BFA","BGD","BGR","BHR","BHS","BIH","BLM","BLR","BLZ","BMU","BOL","BRA",
      "BRB","BRN","BTN","BVT","BWA","CAF","CAN","CCK","CHE","CHL","CHN","CIV","CMR","COD","COG","COK",
      "COL","COM","CPV","CRI","CUB","CUW","CXR","CYM","CYP","CZE","DEU","DJI","DMA","DNK","DOM","DZA",
      "ECU","EGY","ERI","ESH","ESP","EST","ETH","FIN","FJI","FLK","FRA","FRO","FSM","GAB","GBR","GEO",
      "GGY","GHA","GIB","GIN","GLP","GMB","GNB","GNQ","GRC","GRD","GRL","GTM","GUF","GUM","GUY","HKG",
      "HMD","HND","HRV","HTI","HUN","IDN","IMN","IND","IOT","IRL","IRN","IRQ","ISL","ISR","ITA","JAM",
      "JEY","JOR","JPN","KAZ","KEN","KGZ","KHM","KIR","KNA","KOR","KWT","LAO","LBN","LBR","LBY","LCA",
      "LIE","LKA","LSO","LTU","LUX","LVA","MAC","MAF","MAR","MCO","MDA","MDG","MDV","MEX","MHL","MKD",
      "MLI","MLT","MMR","MNE","MNG","MNP","MOZ","MRT","MSR","MTQ","MUS","MWI","MYS","MYT","NAM","NCL",
      "NER","NFK","NGA","NIC","NIU","NLD","NOR","NPL","NRU","NZL","OMN","PAK","PAN","PCN","PER","PHL",
      "PLW","PNG","POL","PRI","PRK","PRT","PRY","PSE","PYF","QAT","REU","ROU","RUS","RWA","SAU","SDN",
      "SEN","SGP","SGS","SHN","SJM","SLB","SLE","SLV","SMR","SOM","SPM","SRB","SSD","STP","SUR","SVK",
      "SVN","SWE","SWZ","SXM","SYC","SYR","TCA","TCD","TGO","THA","TJK","TKL","TKM","TLS","TON","TTO",
      "TUN","TUR","TUV","TWN","TZA","UGA","UKR","UMI","URY","USA","UZB","VAT","VCT","VEN","VGB","VIR",
      "VNM","VUT","WLF","WSM","YEM","ZAF","ZMB","ZWE","ABW"
    ),
    region = c(
      "SA","SSA","LAM","GLO","GLO","GLO","MEA","LAM","GLO","GLO","GLO","GLO","LAM","GLO","GLO","GLO",
      "SSA","GLO","SSA","LAM","SSA","SA","GLO","MEA","LAM","GLO","LAM","GLO","LAM","LAM","LAM","LAM",
      "LAM","SEA","SA","GLO","SSA","SSA","GLO","GLO","GLO","LAM","EA","SSA","SSA","SSA","SSA","GLO",
      "LAM","SSA","SSA","LAM","LAM","LAM","GLO","LAM","GLO","GLO","GLO","SSA","LAM","GLO","LAM","MEA",
      "LAM","MEA","SSA","MEA","GLO","GLO","SSA","GLO","GLO","GLO","GLO","GLO","GLO","SSA","GLO","GLO",
      "GLO","SSA","GLO","SSA","LAM","SSA","SSA","SSA","GLO","LAM","GLO","LAM","LAM","GLO","LAM","EA",
      "GLO","LAM","GLO","LAM","GLO","SEA","GLO","SA","GLO","GLO","MEA","MEA","GLO","MEA","GLO","LAM",
      "GLO","MEA","EA","GLO","SSA","GLO","SEA","GLO","LAM","EA","MEA","SEA","MEA","SSA","MEA","LAM",
      "GLO","SA","SSA","GLO","GLO","GLO","EA","LAM","MEA","GLO","GLO","SSA","SA","LAM","GLO","GLO",
      "SSA","GLO","SEA","GLO","EA","GLO","SSA","SSA","LAM","LAM","SSA","SSA","SEA","SSA","SSA","GLO",
      "SSA","GLO","SSA","LAM","GLO","GLO","GLO","SA","GLO","GLO","MEA","SA","LAM","GLO","LAM","SEA",
      "GLO","GLO","GLO","LAM","EA","GLO","LAM","MEA","GLO","MEA","SSA","GLO","GLO","SSA","MEA","SSA",
      "SSA","SEA","GLO","SSA","GLO","GLO","SSA","LAM","GLO","SSA","GLO","GLO","SSA","SSA","LAM","GLO",
      "GLO","GLO","SSA","LAM","SSA","MEA","LAM","SSA","SSA","SEA","GLO","GLO","GLO","SEA","GLO","LAM",
      "MEA","GLO","GLO","EA","SSA","SSA","GLO","GLO","LAM","LAM","GLO","GLO","LAM","LAM","LAM","LAM",
      "SEA","GLO","GLO","GLO","MEA","SSA","SSA","SSA","LAM"
    ), stringsAsFactors = FALSE)

  # Build region vector in ISO order
  regions <- iso2region$region

  # Check that all regions exist in x
  missingRegions <- setdiff(unique(regions), getItems(x, dim = "region"))
  if (length(missingRegions) > 0) {
    stop("x does not contain these regions: ", paste(missingRegions, collapse = ", "))
  }

  # Expand from object with ordered region vector and rename to countries
  out <- x[regions, , ]
  getItems(out, dim = "region") <- iso2region$iso
  out <- out[isoCountry$x, , ]

  return(out)
}
