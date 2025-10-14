#' @title calcNutritionAttributes
#'
#' @description Calculates nutrition attributes of food products, i.e. calorie and protein
#'              supply of a product dedicated to food use.
#'
#' @return magpie object
#' @author Benjamin Bodirsky
#'
#' @examples
#' \dontrun{
#' calcOutput("NutritionAttributes", aggregate = FALSE)
#' }
#'
#' @importFrom madrat calcOutput
#' @importFrom magclass dimSums

calcNutritionAttributes <- function() {

  massbalance <- calcOutput("FAOmassbalance", aggregate = FALSE)

  # global sum of household massbalance of reactive nitrogen (nr) and generalizable energy (ge)
  household   <- dimSums(massbalance[, , "households"][, , c("ge", "nr")], dim = c(1, 3.2), na.rm = TRUE)

  # milling still includes brans, which have to be considered for massbalance correction calcs
  # therefore we use flour for the correction
  fooduseFlour <- dimSums(massbalance[, , c("food", "flour1")][, , c("ge", "nr")], dim = c(1, 3.2), na.rm = TRUE)

  if (any(household > fooduseFlour)) {
    message(paste("The following items violate massbalance constraints: ",
                  paste(unique(dimnames(household)[[3]][which(household > fooduseFlour, arr.ind = TRUE)[, 3]]),
                        collapse = " "),
                  ". Violating items are corrected through household balance flow.",
                  "Epecially for livestock products it seems like FAOSTAT's nutrient conversion may be very",
                  "different from our reasonable values."))
    household[household > fooduseFlour] <- fooduseFlour[household > fooduseFlour]
  }

  # but now we want them as share of the milling quantity
  fooduseMilling <- dimSums(massbalance[, , c("food", "milling")][, , c("dm")],
                            dim = c(1, 3.2, 3.3),
                            na.rm = TRUE)
  out <- household / fooduseMilling
  out[out == Inf] <- 0
  out[is.nan(out)] <- 0
  out[, , "ge"] <- out[, , "ge"] / 4.184
  out[, , "nr"] <- out[, , "nr"] * 6.258
  dimnames(out)[[3]] <- gsub(".nr", dimnames(out)[[3]], replacement = ".protein")
  dimnames(out)[[3]] <- gsub(".ge", dimnames(out)[[3]], replacement = ".kcal")

  # Treatment of special cases:
  # replace zeros that have values before and/or after
  # with the mean value over the whole period
  for (i in getItems(out, dim = 3)) {
    if (any(out[, , i] == 0) && !all(out[, , i] == 0)) {
      tmp             <- out[, , i]
      tmp[tmp == 0]   <- NA
      tmp[is.na(tmp)] <- mean(tmp, na.rm = TRUE)
      out[, , i]      <- tmp
    }
  }

  # add missing products and delete unconsidered ones
  kall <- findset("kall")
  out  <- out[, , kall]
  out[, , "scp"][, , "kcal"]    <- 20.9 / 4.184
  out[, , "scp"][, , "protein"] <- 0.45

  # add missing years
  years        <- findset("time")
  lastyear     <- paste0("y", max(getYears(out, as.integer = TRUE)))
  missingyears <- setdiff(years, getYears(out))
  out          <- add_columns(x = out, addnm = missingyears, dim = 2.1)
  out[, missingyears, ] <- setYears(out[, lastyear, ], NULL)

  return(list(x = out, weight = NULL,
              unit = "mio kcal / t DM, t Protein / t DM",
              description = paste("Values from FAO Food Supply. ",
                                  "Describe final calory and protein supply ",
                                  "of a product dedicated for fooduse.")))
}
