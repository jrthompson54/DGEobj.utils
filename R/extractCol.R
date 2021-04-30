#' Extract a named column from a series of df or matrices
#'
#' Take a named list of dataframes where each dataframe has the same
#' column names (e.g. a list of topTable dataframes), then extract
#' the named column from each dataframe and return a matrix.  The name of
#' each dataframe is used as the column name in the resulting table.
#'
#' The common use case for this is to provide a list of topTable
#' data frames and extract one column from each file to create
#' a matrix of LogRatios or P-values (genes x contrasts)..
#'
#' This should work as long as the requested column name is present in every
#' dataframe.  The default robust = TRUE should be used unless it has been
#' verified that each dataframe in the input list has the same row count and row
#' order.
#'
#' @param contrastList A list of data.frames which all have the same colnames and same row counts.
#'     The dataframes in the list should have geneIDs as rownames.
#' @param colName The name of the data column to extract into a matrix.
#' @param robust Default = TRUE; TRUE forces use of a joins to merge columns
#'     which is more reliable and allows for combination of contrasts from different
#'     projects, but may not return items in the same row order as the source
#'     table. Setting to FALSE invokes a cbind() approach that requires all
#'     data.frames to have the same row count and row order but preserves the
#'     original row order.
#'
#' @return A dataframe containing the extracted columns
#'
#' @examples
#'    dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'    TopTableList <- DGEobj::getType(dgeObj, type = "topTable")
#'    MyPvalues    <- extractCol(TopTableList, colName = "P.Value")
#'    head(MyPvalues)
#'
#' @export
extractCol <- function(contrastList, colName, robust = TRUE){
    assertthat::assert_that(!missing(contrastList),
                            !is.null(contrastList),
                            length(unique(lapply(contrastList, colnames))) == 1,
                            length(unique(lapply(contrastList, dim))) == 1,
                            msg = "contrastList must be a list of data.frames which all have the same colnames and same row counts.")
if (any(is.null(robust),
            !is.logical(robust),
            length(robust) != 1)) {
        warning("robust must be a singular logical value. Assigning default value TRUE.")
        robust = TRUE
    }
    ifelse(robust,
           return(.extractCol2(contrastList, colName)),
           return(.extractCol1(contrastList, colName))
    )
}


# Helper functions
.extractCol1 <- function(contrastList, colName){
    assertthat::assert_that("list" %in% class(contrastList),
                            !is.null(names(contrastList)),
                            msg = "contrastList must be a named list.")
    assertthat::assert_that("character" %in% class(colName),
                            msg = "colName must be a column in the data of class 'character'.")

    MyMatrix = do.call("cbind", lapply(contrastList, `[[`, colName))
    # Get gene ids from first df
    rownames(MyMatrix) <- rownames(contrastList[[1]])
    # Transfer the contrast names
    colnames(MyMatrix) <- names(contrastList)
    return(as.data.frame(MyMatrix))
}


.extractCol2 <- function(contrastList, colName){
    # Support combining topTable data from different DGEobjs
    assertthat::assert_that("list" %in% class(contrastList),
                            !is.null(names(contrastList)),
                            msg = "contrastList must be a named list.")
    assertthat::assert_that("character" %in% class(colName),
                            msg = "colName must be a column in the data of class 'character'.")

    for (i in 1:length(contrastList)) {
        newdat <- cbind(rowid = rownames(contrastList[[i]]), data.frame(contrastList[[i]], row.names = NULL))
        newdat <- newdat[, c("rowid", colName)]

        if (i == 1) {
            dat <- newdat
        } else {
            dat <- merge(x = dat, y = newdat, by = "rowid", all = TRUE, sort = FALSE)
        }
        colnames(dat)[i + 1] <- names(contrastList[i])
    }
    dat <- data.frame(dat[,-c(1)], row.names = dat[,c(1)])
    return(dat)
}
