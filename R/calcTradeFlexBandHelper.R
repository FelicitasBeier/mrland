#' calcTradeFlexBandHelper
#'
#' Calculate the rolling-range flexibility band of trade ratios observed over a window of years.
#'
#' Given a tidy data.frame/tibble with columns for exporter (ex), importer (im),
#' item (ItemCodeItem), year and a numeric value column (default `.value`),
#' this function computes the rolling RANGE (max - min) over `windowYears` for each
#' group and returns the mean, max and min of that rolling range per group. The mean
#' variant (mean over windows of the rolling range) is used by the model as the
#' flexibility window; it is empirically monotone in window length on the historical
#' data. The max variant is provably monotone (range is a monotone set function) and
#' wider; min is the narrowest. All three are written to the input file.
#'
#' @param dataIn data.frame or tibble; tidy table with columns specified in `groupVars`, `yearCol` and `valueCol`.
#' @param windowYears integer; window length (number of years) for the rolling range (must be >= 1).
#' @param groupVars character; grouping columns. Default: c("ex", "im", "ItemCodeItem").
#' @param yearCol character; name of the year column. Default: "Year".
#' @param valueCol character; name of the numeric value column. Default: ".value".
#' @return tibble with one row per group and columns: mean<years>, max<years>, min<years>.
#' @examples
#' \dontrun{
#' ratiodf <- as.data.frame(collapseNames(ratio), rev = 2)
#' calcTradeFlexBandHelper(ratiodf, windowYears = 5)
#' }
#' @author David M Chen
#' @export
#' @importFrom dplyr arrange group_by mutate summarise across all_of .data
#' @importFrom zoo rollapply
calcTradeFlexBandHelper <- function(dataIn,
                                    windowYears,
                                    groupVars = c("ex", "im", "ItemCodeItem"),
                                    yearCol = "Year",
                                    valueCol = ".value") {
  stopifnot(length(windowYears) == 1 && is.numeric(windowYears) && windowYears >= 1)

  meanCol <- paste0("mean", windowYears)
  maxCol  <- paste0("max", windowYears)
  minCol  <- paste0("min", windowYears)

  dataIn <- as.data.frame(collapseNames(dataIn), rev = 2)

  out <- dataIn %>%
    dplyr::arrange(.data[[yearCol]]) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(groupVars))) %>%
    dplyr::mutate(
      rollrange = zoo::rollapply(.data[[valueCol]],
                                 width = windowYears,
                                 FUN = function(x) diff(range(x)),
                                 fill = NA,
                                 align = "right")
    ) %>%
    dplyr::summarise(
      !!meanCol := mean(.data[["rollrange"]], na.rm = TRUE),
      !!maxCol  := max(.data[["rollrange"]],  na.rm = TRUE),
      !!minCol  := min(.data[["rollrange"]],  na.rm = TRUE),
      .groups = "drop"
    )

  out <- as.magpie(out,  spatial = c(1, 2), temporal = NULL)

  return(list(x = out,
              weight = NULL,
              unit = "ratio",
              isocountries = FALSE,
              description = "Flexibility band (rolling range) of historical import supply ratios"))
}
