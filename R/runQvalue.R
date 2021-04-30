#' Calculate and add q-value and lFDR to dataframe
#'
#' Takes an list of contrasts (e.g. topTable output or other dataframes that contain
#' a p-value column).  Adds a q-value and local FDR (lFDR) column to each dataframe.
#'
#' The qvalue package from John Storey at Princeton takes a list of p-values and
#' calculates a q-value and a Local FDR (lFDR). The q-value is essentially a less
#' conservative FDR estimate compared to the default Benjamini-Hochberg FDR
#' produced by topTable analysis (i.e. will give more differential genes at
#' the same nominal cutoff). The q-value function also produces a Local FDR
#' (lFDR) column which answers a slightly different and possibly more relevant
#' question. The BH FDR (adj.P.Val in topTable data.frames) and q-value gives
#' the false discovery rate is for a list of genes at a given threshold.
#' The local FDR attempts to answer the question: what is the probability that
#' this particular gene is a false discovery?
#' See \doi{10.1007/978-3-642-04898-2_248} for a brief introduction
#' to FDRs and q-values.
#'
#' @param contrastList A list of dataframes with a p-value column (all tables
#'   must use the same colname for the p-value column.)
#' @param pvalField Define the colname of the p-value field in
#'   each dataframe. Not needed if using topTable output. (Optional. Default = "P.Value")
#' @param ... Optional arguments passed to the qvalue function (See ?qvalue)
#'
#' @return The input contrastList now containing q-value and lFDR columns
#'   in each dataframe.
#'
#' @examples
#'    dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'    contrastList <- DGEobj::getType(dgeObj, type = "topTable")
#'    contrastList <- lapply(contrastList, dplyr::select,
#'                           -Qvalue,
#'                           -qvalue.lfdr)
#'    colnames(contrastList[[1]])
#'
#'    contrastList <- runQvalue(contrastList)
#'
#'    # note new columns added
#'    colnames(contrastList[[1]])
#'
#' @importFrom qvalue qvalue
#' @importFrom assertthat assert_that
#'
#' @export
runQvalue <- function(contrastList, pvalField = "P.Value", ...){
    # Add Q-values to each topTable dataframe in contrastList
    assertthat::assert_that(!missing(contrastList),
                            !is.null(contrastList),
                            "list" %in% class(contrastList),
                            msg = "contrastList must be of class 'list'.")

    contrastNames = names(contrastList)

    for (i in 1:length(contrastList)) {
        assertthat::assert_that(exists(pvalField, contrastList[[i]]),
                                msg = "pvalField must exist as an item in contrastList.")
        p = contrastList[[i]][, pvalField]
        q = qvalue::qvalue(p, lambda = 0, ...)
        # Add the q-value and lFDR columns to the topTable df
        contrastList[[i]]$Qvalue = q$qvalues
        contrastList[[i]]$qvalue.lfdr = q$lfdr
        # Add documentation
        attr(contrastList[[i]], "qvalue") = TRUE
    }

    return(contrastList)
}
