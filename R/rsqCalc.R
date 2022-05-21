#' Calculate R-squared for each gene fit
#'
#' Takes a Log2CPM numeric matrix and MArrayLM fit object from limma's lmFit()
#' and calculates R-squared for each gene fit.
#'
#' @param normMatrix A normalized log2cpm matrix
#' @param fit A MArrayLM object from limma's lmFit()
#'
#' @return A vector of R-squared values for each gene fit.
#'
#' @examples
#' \dontrun{
#'    # NOTE: Requires the edgeR package
#'
#'    dgeObj    <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'    log2cpm   <- convertCounts(dgeObj$counts, unit = "cpm", log=TRUE, normalize = "tmm")
#'    fitObject <- dgeObj$ReplicateGroupDesign_fit
#'    rsq       <- rsqCalc (log2cpm, fitObject)
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom stringr str_c
#'
#' @export
rsqCalc <- function(normMatrix, fit) {
    assertthat::assert_that(!missing(normMatrix),
                            !is.null(normMatrix),
                            any(c("data.frame", "matrix") %in% class(normMatrix)),
                            msg = "normMatrix must be of class 'data.frame' or 'matrix'.")
    assertthat::assert_that(!missing(fit),
                            !is.null(fit),
                            "MArrayLM" %in% class(fit),
                            msg = "fit must be of class 'MArrayLM'.")
    assertthat::assert_that(is.numeric(as.matrix(normMatrix)),
                            msg = "All of the entries in normMatrix must be numeric.")

    sst <- rowSums(normMatrix^2)
    ssr <- sst - fit$df.residual*fit$sigma^2
    return(ssr)
}
