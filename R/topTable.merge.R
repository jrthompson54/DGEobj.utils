#' Merge specified topTable df cols
#'
#' Take a named list of topTable dataframes and cbinds the requested columns
#' from each file.  To avoid column name conflicts the names are used as suffixes
#' to the colnames. Although written for topTable data, this should work on any
#' named list of dataframes where each member of the list has the same columns.
#'
#' @param contrastList A named list of topTable data.frames which all have the same colnames and same row counts.
#' The dataframes in the list should have rownames (geneIDs).
#' @param colNames The list of column names of the data column to extract to a
#'   matrix (Default = c("logFC", "AveExpr", "P.Value", "adj.P.Val"))
#' @param digits Number of decimal places for specified columns. Should be same
#' length as colNames. (Default = c(2, 2, 4, 3)). If one value supplied, it is used
#' for all columns.
#'
#' @return A matrix containing the extracted columns.
#'
#' @examples
#' dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#' contrastList <- DGEobj::getType(dgeObj, type = "topTable")
#'
#' mergedData <- topTable.merge(contrastList)
#' colnames(mergedData)
#'
#' @importFrom stringr str_c
#'
#' @export
topTable.merge <- function(contrastList,
                           colNames = c("logFC",
                                        "AveExpr",
                                        "P.Value",
                                        "adj.P.Val"),
                           digits = c(2, 2, 4, 3)) {

    assertthat::assert_that(!missing(contrastList),
                            "list" %in% class(contrastList),
                            "data.frame" %in% class(contrastList[[1]]),
                            !is.null(names(contrastList)),
                            length(unique(lapply(contrastList, colnames))) == 1,
                            length(unique(lapply(contrastList, dim))) == 1,
                            msg = "contrastList must be specified, be of class 'list' and be a named list specifically, and include items of class 'data.frame'.")
    assertthat::assert_that(length(digits) %in% c(1, length(colNames)),
                            msg = "digits must be either of length 1 or the same length as colNames.")

    if (length(digits) == 1) { # Expand the digits vector
        digits <- rep(digits, length(colNames))
    }

    # Add contrast suffix to each colname
    contrastNames <- names(contrastList)

    # Get the first set of columns
    dat <- as.data.frame(extractCol(contrastList, colName = colNames[1], robust = TRUE))
    colnames(dat) <- stringr::str_c(colNames[1], "_", colnames(dat))
    dat <- round(dat, digits[1])
    dat <- cbind(rowID = rownames(dat), data.frame(dat, row.names = NULL))

    if (length(colNames)  > 1) {
        for (i in 1:length(colNames)) {
            dat2 <- as.data.frame(extractCol(contrastList, colName = colNames[i], robust = TRUE))
            # Add datatype as prefix on colname e.g. logFC_contrastname
            colnames(dat2) <- stringr::str_c(colNames[i], "_", colnames(dat2))
            dat2 <- round(dat2, digits[i])
            dat2 <- cbind(rowID = rownames(dat2), data.frame(dat2, row.names = NULL))
            if (i == 1) {
                dat <- dat2
            } else {
                dat <- merge(x = dat, y = dat2, by = "rowID", all.x = TRUE, sort = FALSE)
            }
        }
    }

    dat <- data.frame(dat[,-c(1)], row.names = dat[,c(1)])
    return(dat)
}
