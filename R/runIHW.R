#' Apply Independent Hypothesis Weighting (IHW) to a list of topTable dataframes
#'
#' This is a wrapper around the independent hypothesis weighting package that
#' takes a list of topTable data frames and applies Independent Hypothesis
#' Weighting (IHW) to each topTable data frame in the list.
#'
#' IHW is a method developed by N. Ignatiadis (http://dx.doi.org/10.1101/034330)
#' to weight FDR values based on a covariate (AveExpr in this case).
#'
#' The IHW FDR values are added as additional columns to the topTable data frames.
#'
#' Function runIHW is normally called by runContrasts with argument IHW=T.  It
#' can also be used independently on a list of topTable dataframes.  A list of
#' topTable dataframes is conveniently retrieved with the DGEobj::getType
#' function with the type argument set to "topTable".
#'
#' This function expects the following columns are present in each data frame:
#' P.value, adj.P.Val, AveExpr.
#'
#' Note that it is impractical to run IHW on a list of genes less than ~5000.
#' Operationally, IHW breaks the data into bins of 1500 genes for the analysis.
#' If bins = 1, IHW converges on the BH FDR value. Instead, run IHW on the
#' complete set of detected genes from topTable (not topTreat) results.
#'
#' @param contrastList A named list of topTable dataframes.
#' @param alpha Alpha should be the desired FDR level to interrogate (range 0-1; Default = 0.1)
#' @param FDRthreshold Threshold value for the p-values of a dataframe (Default = 0.1)
#' @param ... other arguments are passed directly to the ihw function (see ?ihw)
#'
#' @return A list of lists.  The first element is the original contrastList with
#'   additional IHW columns added to each dataframe. The topTable dataframes
#'   will contain additional columns added by the IHW analysis and prefixed with
#'   "ihw." The second list element is the IHW result dataframe.
#'
#' @examples
#' if (requireNamespace("IHW", quietly = TRUE)) {
#'    dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'    contrastList <- DGEobj::getType(dgeObj, type = "topTable")
#'    contrastList <- lapply(contrastList, dplyr::select,
#'                           -ihw.adj_pvalue,
#'                           -ihw.weight,
#'                           -ihw.weighted_pvalue)
#'    colnames(contrastList[[1]])
#'    contrastList <- runIHW(contrastList)
#'    # note new columns added
#'    colnames(contrastList[["contrasts"]][[1]])
#' }
#'
#' @export
runIHW <- function(contrastList,
                   alpha = 0.1,
                   FDRthreshold = 0.1,
                   ...){
    assertthat::assert_that(requireNamespace("IHW", quietly = TRUE),
                            msg = "IHW package is required to apply Independent Hypothesis Weighting (IHW) to the given list of topTable dataframes")

    assertthat::assert_that(!missing(contrastList),
                            is.list(contrastList),
                            msg = "contrastList must be specified and should be of class 'List'.")
    if (any(is.null(alpha),
            !is.numeric(alpha),
            length(alpha) != 1,
            alpha < 0,
            alpha > 1)) {
        warning("alpha must be a singular numeric value between 0 and 1. Assigning default value 0.1")
        alpha = 0.1
    }

    if (any(is.null(FDRthreshold),
            !is.numeric(FDRthreshold),
            length(FDRthreshold) != 1,
            FDRthreshold < 0,
            FDRthreshold > 1)) {
        warning("FDRthreshold must be a singular numeric value between 0 and 1. Assigning default value 0.1")
        FDRthreshold = 0.1
    }

    getProportion <- function(ttdf, threshold) {
        # Get the proportion for one df
        bhfdrproportion <- (sum(ttdf$adj.P.Val <= threshold)) / nrow(ttdf)
    }

    runIHWon1DF <- function(ttdf, alpha, proportion, ...){
        assertthat::assert_that(!is.null(ttdf$P.Value),
                                !is.null(ttdf$AveExpr),
                                msg = "The topTable dataframes in contrastList must have both P.Value and AveExpr columns.")
        # Run IHW on one df
        # Return an ihwResult object
        do.call("require", list("IHW"))
        IHWresult <- tryCatch({
            do.call("ihw", list(ttdf$P.Value,
                                covariates = ttdf$AveExpr,
                                alpha = alpha,
                                ...))
        },
        error = function(e) {
            message("Unexpected error: ", e$message, " happened during runIHWon1DF execution")
            return(NULL)
        })

    }

    # Run IHW on each dataframe, collect the result objects in a list which
    # is added to the result object.
    proportion    <- sapply(contrastList, getProportion, threshold = FDRthreshold)
    ihwList       <- list()
    contrastNames <- names(contrastList)

    for (i in 1:length(contrastList)) {
        ihwResult <- runIHWon1DF(contrastList[[i]],
                                 alpha = alpha,
                                 proportion = proportion[i], ...)
        if (!is.null(ihwResult)) {
            # Capture the ihwResult object
            ihwList[[i]] <- ihwResult
            contrastList[[i]] <- cbind(contrastList[[i]],
                                       ihwResult@df[,2:4])
            # Prefix the colnames of those three columns with "ihw."
            cnames <- colnames(contrastList[[i]])
            numcol <- length(cnames)
            cnames[(numcol - 2):numcol] <- paste("ihw.", cnames[(numcol - 2):numcol], sep = "")
            colnames(contrastList[[i]]) <- cnames

            # Add documentation
            attr(contrastList[[i]], "ihw") <-  TRUE
        }
    }

    list(contrasts = contrastList, ihwObj = ihwList)
}
