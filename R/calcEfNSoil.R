#' @title calcEmisNitrogenPast
#' @description Emission factors from cropland soils.
#'
#' @param method If IPCC, using the ipcc emission factors as share of applied N inputs.
#'               If Nloss, as share of cropland budget surplus.
#' @return List of magpie object with results on country level,
#'         weight on country level, unit and description.
#' @author Benjamin Leon Bodirsky
#' @examples
#'
#' \dontrun{
#' calcOutput("EmisNitrogenPast")
#' }

calcEfNSoil <- function(method = "IPCC_reg") {

  if (method == "Nloss") {
    surplus <- dimSums(calcOutput("NitrogenBudgetCropland")[, , "surplus"],
                       dim = 1)
    emis    <- dimSums(calcOutput("EmisNitrogenCroplandPast", method = "IPCC"),
                       dim = c(1, 3.2))
    efNSoil <- emis / surplus
    weight  <- surplus
    #extend beyond last year
    efNSoil <- toolHoldConstantBeyondEnd(efNSoil)
    weight <- toolHoldConstantBeyondEnd(weight)
  } else if (method == "IPCC") {
    efNSoil <- setYears(calcOutput("IPCCefNSoil", aggregate = "GLO")[, "y2010", ],
                        NULL)
    weight  <- NULL
  } else if (method == "IPCC_reg") {
    tmp     <- calcOutput("IPCCefNSoil", aggregate = FALSE,
                          supplementary = TRUE)
    efNSoil <- tmp$x
    weight  <- tmp$weight
    #extend beyond last year
    efNSoil <- toolHoldConstantBeyondEnd(efNSoil)
    weight <- toolHoldConstantBeyondEnd(weight)

  } else {
    stop("method unknown")
  }

  return(list(x = efNSoil,
              weight = weight,
              unit = "Share",
              description = paste0("Emission factors from cropland soils. ",
                                   "If IPCC, using the ipcc emission factors as
                                   share of applied N inputs. ",
                                   "If Nloss, as share of cropland budget surplus.")))
}
