#' Summarize a contrast list
#'
#' Takes a contrast list produced by runContrasts. Defaults are provided to specify columns to
#' summarize and thresholds for each column, though they can be adjusted. A fold change
#' threshold can optionally be specified. The function queries the topTable results and
#' returns a dataframe with the summary results, but only includes gene counts that meet
#' the specified conditions.
#'
#' Any specified column names that don't exist will be ignored. Normally the
#' defaults cover all the p-value and FDR related columns. However, a fcThreshold
#' can be added and the p-value/FDR thresholds can be modified using the fcThreshold
#' and sigThresholds arguments, respectively.
#'
#' @param contrastList A list of topTable dataframes.
#' @param columns Vector of column names to summarize from topTable dataframes.
#'   Default = c("P.Value", "adj.P.Val", "Qvalue", "qvalue.lfdr", "ihw.adj_pvalue")
#' @param sigThresholds Thresholds to use for each column specified in columns
#'   Must be same length at columns argument.
#'   Default = c(0.01, 0.05, 0.05, 0.05, 0.05)
#' @param fcThreshold Fold-change threshold (absolute value, not logged.)
#'
#' @return data.frame with one summary row per contrast.
#'
#' @importFrom assertthat assert_that
#'
#' @examples
#'    dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'    contrastList <- DGEobj::getType(dgeObj, type = "topTable")
#'
#'    #all defaults
#'    sigSummary <- summarizeSigCounts(contrastList)
#'
#'    #add the fold-chage threshold
#'    sigSummary <- summarizeSigCounts(contrastList, fcThreshold = 2)
#'
#'    #change the significance thresholds
#'    sigSummary <- summarizeSigCounts(contrastList,
#'                                     sigThresholds = c(0.01, 0.1, 0.1, 0.1, 0.1))
#'
#' @export
summarizeSigCounts <- function(contrastList,
                               columns       = c("P.Value", "adj.P.Val", "Qvalue", "qvalue.lfdr", "ihw.adj_pvalue"),
                               sigThresholds = c(0.01, 0.05, 0.05, 0.05, 0.05),
                               fcThreshold   = 0) {

    assertthat::assert_that(length(columns) == length(sigThresholds),
                            msg = "Supplied sigThresholds should be same length as supplied columns.")

    # Functions
    getSigCounts <- function(df, columns, thresholds, fcThreshold) {
        counts <- list()
        for (i in 1:length(columns)) {
            idx <- df[columns[i]] <= thresholds[i]
            if (fcThreshold > 0) {
                fcidx <- abs(df$logFC) >= log2(fcThreshold)
                idx <- idx & fcidx
            }

            counts[columns[i]] <- sum(idx)
        }
        return(unlist(counts))
    }

    columns <- columns[columns %in% colnames(contrastList[[1]])]
    sigThresholds <- sigThresholds[columns %in% colnames(contrastList[[1]])]

    myrows <- list()
    for (i in 1:length(contrastList)) {
        myrows[[i]] <- getSigCounts(contrastList[[i]], columns, sigThresholds, fcThreshold)
    }

    DF <- do.call("rbind", myrows)

    rownames(DF) <- names(contrastList)
    colnames(DF) <- columns

    return(DF)
}
