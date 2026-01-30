#' @title calcTradeBilateralFAOHarmonized
#'
#' @description Harmonizes bilateral trade values to FAO mass balance imports
#' often termed "Import Dependency Ratio"
#' @param yearly whether to calculate yearly data or only magpie 5year timesteps
#' @return Self import to dupply ratio
#' @author David M Chen
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link[mrcommons]{calcFAOmassbalance}}
#' @examples
#' \dontrun{
#' a <- calcTradeBilateralFAOHarmonized()
#' }
#'
calcTradeBilateralFAOHarmonized <- function(yearly = TRUE) {

  mb <- calcOutput("FAOmassbalance", yearly = TRUE, aggregate = FALSE)
  tm <- calcOutput("FAOBilateralTrade", output = "qty",
                   products = "kcr", fiveYear = FALSE,
                   aggregate = FALSE)
  tm1 <- calcOutput("FAOBilateralTrade", output = "qty",
                    products = "kli", fiveYear = FALSE,
                    aggregate = FALSE)
  tm2 <- calcOutput("FAOBilateralTrade", output = "qty",
                    products = "kothers", fiveYear = FALSE,
                    aggregate = FALSE)
  tmfo <- calcOutput("TimberTradeBilateral", aggregate = FALSE)
  tm <- mbind(tm, tm1)
  tm <- mbind(tm, tm2)
  rm(tm1, tm2)

  ## hold the timber bilat constant into past and future, as they just missing a few years
  tmfo <- toolHoldConstant(tmfo,
                           years = c(tail(getYears(tmfo, as.integer = TRUE), n = 1):
                                       tail(getYears(tm, as.integer = TRUE), n = 1)))
  tmfo <- toolHoldConstant(tmfo, years = c(1995:1996))
  tmfo[, 1995, ] <- tmfo[, 1996, ] <- tmfo[, 1997, ]
  cyears <- intersect(getYears(tm), getYears(tmfo))

  tm <- mbind(tm, tmfo)
  # add fish and ethanol to trade matrix with 0's to be filled by mb
  tm <- add_columns(tm, addnm = c("fish", "ethanol"), dim = 3.1, fill = 0)

  tmi <- dimSums(tm, dim = 1.2)

  cyears <- intersect(getYears(tm), getYears(mb))
  citems <- intersect(getItems(tm, dim = 3), getItems(mb, dim = 3.1))

  mb <- mb[, cyears, citems]
  mbx <- mb[, , "export"][, , "dm", drop = TRUE]
  mbi <- mb[, , "import"][, , "dm", drop = TRUE]

  # sodistillers grain, scp, still  missing completely from the bilateral trade matrix
  # the first is a made up producsts so don't exist in FAO either, but
  # future trade should still be allowed?

  # here we look at difference between trade matrix and mass balance and scale the trade matrix to mb

  imDiff <- tmi[, cyears, citems] / mbi[, cyears, citems]

  imDiff[is.na(imDiff)] <- 1

  imDiff[is.infinite(imDiff)] <- NA

  # scale up imports

  # rename exporters so that we scale up  the Wimporters and not exporters
  getItems(tm, dim = 1.2) <- paste0(getItems(tm, dim = 1.2), "2")

  tmS  <- tm[, cyears, citems] / imDiff
  tmS[is.infinite(tmS)] <- 0
  tmS[is.na(tmS)] <- 0
  # for cases where there is no tmi but trade reported in mbi we still want to match mbi
  # so as assumption split the mbi imports across the bilateral patterns first by taking average of adjacent years
  # then by trade patterns of all products

  # fill based on function wwhich first tries with average between both adjacent years,
  # then tries one or the other, and then both again
  #  only for where mbi >0 and sum(tm, dim=1.2) i.e. total imports  =0

  .fillMiss <- function(tm, year = "both") {
    # amounts missing
    totMiss <- ifelse(dimSums(tm, dim = 1.2) == 0 & mbi[, cyears, citems] > 0,
                      mbi[, cyears, citems],
                      0)
    yearAdjTrPatterns <- new.magpie(cells_and_regions = getItems(tm, dim = 1),
                                    years = getYears(tm),
                                    names = getNames(tm),
                                    fill = 0)

    for (t in getYears(tm, as.integer = TRUE)[-c(1, length(getYears(tm)))]) {

      if (year == "both") {
        yearAdjTrPatterns[, t, ] <-      (setYears(tm[, t - 1, ] / dimSums(tm[, t - 1, ], dim = 1.2), NULL) +
                                            setYears(tm[, t + 1, ] / dimSums(tm[, t + 1, ], dim = 1.2), NULL)) / 2
      } else if (year == "after") {
        yearAdjTrPatterns[, t, ] <-      setYears(tmS2[, t + 1, ] / dimSums(tmS2[, t + 1, ], dim = 1.2), NULL)
      } else if (year == "before") {
        yearAdjTrPatterns[, t, ] <-      setYears(tmS2[, t - 1, ] / dimSums(tmS2[, t - 1, ], dim = 1.2), NULL)
      }
    }

    toFill <- totMiss * yearAdjTrPatterns
    toFill[is.na(toFill)] <- 0

    out <- tm + toFill
    return(out)
  }

  # make out object
  tmSout <- tmS

  for (i in seq(1:5)) {
    # first both
    tmS2 <- .fillMiss(tmSout, year = "both")
    # then after
    tmS3 <- .fillMiss(tmS2, year = "after")
    # then before
    tmSout <- .fillMiss(tmS3, year = "before")

  }

  # give exporters the 2 suffix
  getItems(mbx, dim = 1) <- paste0(getItems(mbx, dim = 1), "2")

  # for ethanol we will sum the total amount imported globally and distribute according to global exporters ratio
  # also for ethanol we know that exports are 0 globally for 1995-1997 but not 0 for imports
  # so for these years give the same ratio as 1998
  mbx[, c(1995:1997), "ethanol"] <- mbx[, 1998, "ethanol"]
  tmSout[, , "ethanol"] <- mbi[, , "ethanol"] * mbx[, , "ethanol"] /
    dimSums(mbx[, , "ethanol"], dim = 1)

  nonmatch <- dimSums(tmSout, dim = 1.2) - mbi[, cyears, citems]
  missing <- magclass::where(round(nonmatch, 3) < 0)$true$individual
  missing <- as.data.frame(missing)


  # Now we will fill all the others also with global exporters ratio
  for (i in seq_along(missing[, 1])) {
    tmSout[missing$im[i], missing$Year[i], missing$ItemCodeItem[i]] <-
      mbi[missing$im[i], missing$Year[i], missing$ItemCodeItem[i]] *
      mbx[, missing$Year[i], missing$ItemCodeItem[i]] /
      dimSums(mbx[, missing$Year[i], missing$ItemCodeItem[i]], dim = 1)
  }

  getItems(tmSout, dim = 1.2) <- gsub("[0-9]+", "", getItems(tmSout, dim = 1.2))

  return(list(x = tmSout,
              weight = NULL,
              unit = "Mt",
              description = "Trade matrix harmonized with massbalance imports"))
}
