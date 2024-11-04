#' @title calcYieldsLPJmL
#'
#' @description This function extracts yields from LPJmL in the first season (main growing period)
#'              and calculates the yields for the off season (second growing period)
#'
#' @param lpjml         Defines LPJmL version for main crop inputs
#' @param climatetype   Switch between different climate scenarios
#' @param multicropping Multicropping activated (TRUE) or not (FALSE) and
#'                      Multiple Cropping Suitability mask selected
#'                      (mask can be:
#'                      "none": no mask applied (only for development purposes)
#'                      "actual:total": currently multicropped areas calculated from total harvested areas
#'                                      and total physical areas per cell from readLandInG
#'                      "actual:crop" (crop-specific), "actual:irrigation" (irrigation-specific),
#'                      "actual:irrig_crop" (crop- and irrigation-specific),
#'                      "potential:endogenous": potentially multicropped areas given
#'                                              growing conditions from LPJmL
#'                      "potential:exogenous": potentially multicropped areas given
#'                                             GAEZ suitability classification)
#'                      (e.g. TRUE:actual:total; TRUE:none; FALSE)
#'
#' @return magpie object in cellular resolution
#'
#' @author Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("YieldsLPJmL", aggregate = FALSE)
#' }
#'
#' @importFrom madrat toolGetMapping
#' @importFrom withr local_options

calcYieldsLPJmL <- function(lpjml = "ggcmi_phase3_nchecks_bft_e511ac58",
                            climatetype = "GSWP3-W5E5:historical",
                            multicropping = FALSE) {
  # Extract multiple cropping argument information
  areaMask      <- paste(str_split(multicropping, ":")[[1]][2],
                         str_split(multicropping, ":")[[1]][3], sep = ":")
  multicropping <- as.logical(str_split(multicropping, ":")[[1]][1])

  # Increase object size limit
  local_options(magclass_sizeLimit = 1e+12)

  # LPJmL crop types
  lpj2mag     <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mrlandcore")
  cropsLPJmL  <- unique(lpj2mag$LPJmL5)

  # Combine irrigated and rainfed yields for all crops
  irYlds <- list()
  rfYlds <- list()
  yields <- list()
  # Read in by crop because of memory issues
  for (crop in cropsLPJmL) {
    # irrigated yields in irrigated growing period (in tDM/ha)
    irYlds[[crop]] <- calcOutput("LPJmLharmonize", subtype = "crops:pft_harvestc",
                                 # To Do (change once LPJmL runs ready): "cropsIr:pft_harvestc",
                                 subdata = c(crop, "irrigated"),
                                 version = lpjml, climatetype = climatetype,
                                 aggregate = FALSE)
    # rainfed yields in rainfed growing period (in tDM/ha)
    rfYlds[[crop]] <- calcOutput("LPJmLharmonize", subtype = "crops:pft_harvestc",
                                 # To Do (change once new LPJmL runs ready): "cropsRf:pft_harvestc",
                                 subdata = c(crop, "rainfed"),
                                 version = lpjml, climatetype = climatetype,
                                 aggregate = FALSE)
    # irrigated and rainfed yields in main growing period (in tDM/ha)
    yields[[crop]] <- mbind(rfYlds[[crop]], irYlds[[crop]])
    ### To do: check whether this is necessary or whether one mbind is sufficient
  }
  yields  <- mbind(yields)

  # Check for NA's
  if (any(is.na(yields))) {
    stop("calcYieldsLPJmL produced NA yields during run selection")
  }

  # For case of multiple cropping, off-season yield needs to be calculated
  if (multicropping) {
    # Multiple cropping yield increase factor
    increaseFactor <- calcOutput("MulticroppingYieldIncrease",
                                 lpjml = source[["lpjml"]], # nolint: undesirable_function_linter.
                                 climatetype = climatetype,
                                 selectyears = getItems(yields, dim = 2),
                                 aggregate = FALSE)

    # Main-season yield
    mainYield <- yields
    # Off-season yield
    offYield  <- yields * increaseFactor

    # LPJmL to MAgPIE crops
    mainYield <- toolAggregate(mainYield, lpj2mag, dim = 3.1, partrel = TRUE,
                               from = "LPJmL5", to = "MAgPIE")
    offYield  <- toolAggregate(offYield, lpj2mag, dim = 3.1, partrel = TRUE,
                               from = "LPJmL5", to = "MAgPIE")
    # Cap for off-season yield due to numerical reasons
    ### To Do (Feli): write as apply instead
    for (y in getItems(offYield, dim = 2)) {
      for (k in getItems(offYield, dim = 3)) {
        cap <- quantile(offYield[, y, k], 0.999, na.rm = TRUE)
        offYield[, y, k][offYield[, y, k] > cap] <- cap
      }
    }
    ### Question (Jens): Do we want to have a minimum second season yield cap here?
    ### (e.g., all crops/cells with second season yields smaller than 0.5 tDM/ha get 0 second season yield)

    # Multiple cropping suitability
    if (areaMask == "none") {
      suitMC <- calcOutput("MulticroppingCells",
                           sectoral = "lpj",
                           scenario = "potential:exogenous",
                           lpjml = c(crop = source[["lpjml"]]), # nolint: undesirable_function_linter.
                           climatetype = climatetype,
                           selectyears = getItems(yields, dim = 2),
                           ## To Do (Feli): double check which years are available
                           ## and either extrapolate or use constant suitability from iniyear
                           aggregate = FALSE)
      # multiple cropping is allowed everywhere
      suitMC[, , ] <- 1
      # Add grassland to suitMC object with suitability set to 0
      # Note: The grassland growing period is already the whole year, so no multiple
      #       cropping treatment necessary.
      suitMC <- add_columns(suitMC, dim = 3.1, addnm = "grassland", fill = 0)
    } else {
      suitMC <- calcOutput("MulticroppingCells", scenario = areaMask,
                           sectoral = "lpj",
                           lpjml =  c(crop = source[["lpjml"]]), # nolint: undesirable_function_linter.
                           climatetype = climatetype,
                           selectyears = getItems(yields, dim = 2),
                           ## To Do (Feli): double check which years are available
                           ## and either extrapolate or use constant suitability from iniyear
                           aggregate = FALSE)
      # Add grassland to suitMC object with suitability set to 0
      # Note: The grassland growing period is already the whole year, so no multiple
      #       cropping treatment necessary.
      suitMC <- add_columns(suitMC, dim = 3.1, addnm = "grassland", fill = 0)
    }

    # Whole year yields under multicropping (main-season yield + off-season yield)
    yields <- mainYield + offYield * suitMC

  } else {
    # Only main season yields are returned
    yields  <- yields
  }

  return(list(x            = yields,
              weight       = NULL,
              unit         = "tDM per ha",
              description  = "Yields for LPJmL crop types.",
              isocountries = FALSE))
}
