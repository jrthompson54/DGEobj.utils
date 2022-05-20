#' Convert count matrix to CPM, FPKM, FPK, or TPM
#'
#' Takes a count matrix as input and converts to other desired units.  Supported
#' units include CPM, FPKM, FPK, and TPM.  Output units can be logged and/or
#' normalized.  Calculations are performed using edgeR functions except for the
#' conversion to TPM which is converted from FPKM using the formula provided by
#' \href{https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/}{Harold Pimental}.
#'
#' geneLength is a vector where length(geneLength) == nrow(countsMatrix). If a
#' RSEM effectiveLength matrix is passed as input, rowMeans(effectiveLength) is
#' used (because edgeR functions only accept a vector for effectiveLength).
#'
#' Note that log2 values for CPM, TPM, and FPKM employ edgeR's prior.count handling to avoid divide by zero.
#'
#' @param countsMatrix A numeric matrix or dataframe of N genes x M Samples.  All columns must be numeric.
#' @param unit  Required. One of CPM, FPKM, FPK or TPM.
#' @param geneLength A vector or matrix of gene lengths. Required for length-normalized units (TPM, FPKM or FPK).
#'    If geneLength is a matrix, the rowMeans are calculated and used.
#' @param log Default = FALSE.  Set TRUE to return Log2 values.
#'    Employs edgeR functions which use an prior.count of 0.25 scaled by the library size.
#' @param normalize Default = "none". Invokes edgeR's calcNormFactors() for normalization.
#'    Other options are: "TMM", "RLE", "upperquartile" (uses 75th percentile), "TMMwzp" and are case-insensitive.
#' @param prior.count Average count to be added to each observation to avoid taking log of zero.
#'    Used only if log = TRUE. (Default dependent on method; 0 for TPM, 0.25 for CPM and FPKM)
#'    The prior.count is passed to edgeR cpm and rpkm functions and applies to logTPM, logCPM, and logFPKM calculations.
#'
#' @return A matrix in the new unit space
#'
#' @examples
#' \dontrun{
#'   # NOTE: Requires the edgeR package
#'
#'   # Simulate some data
#'   counts <- trunc(matrix(runif(6000, min=0, max=2000), ncol=6))
#'   geneLength <- rowMeans(counts)
#'
#'   # TMM normalized Log2FPKM
#'   Log2FPKM <- convertCounts(counts,
#'                             unit       = "fpkm",
#'                             geneLength = geneLength,
#'                             log        = TRUE,
#'                             normalize  = "tmm")
#'
#'   # Non-normalized CPM (not logged)
#'   RawCPM <- convertCounts(counts,
#'                           unit      = "CPM",
#'                           log       = FALSE,
#'                           normalize = "none")
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr %>%
#'
#' @export
convertCounts <- function(countsMatrix,
                          unit,
                          geneLength,
                          log = FALSE,
                          normalize = "none",
                          prior.count = NULL) {
    assertthat::assert_that(requireNamespace("edgeR", quietly = TRUE),
                            msg = "edgeR package is required to apply edgeR normalization to the given DGEobj")

    assertthat::assert_that(!missing(countsMatrix),
                            !is.null(countsMatrix),
                            any(c("data.frame", "matrix") %in% class(countsMatrix)),
                            length(countsMatrix) != 0,
                            nrow(countsMatrix) != 0,
                            is.numeric(countsMatrix),
                            msg = "countsMatrix must be a numeric matrix or dataframe of N genes x M Samples. All columns must be numeric.")
    assertthat::assert_that(!missing(unit),
                            !is.null(unit),
                            is.character(unit),
                            length(unit) == 1,
                            toupper(unit) %in% c("CPM", "FPKM", "FPK", "TPM"),
                            msg = "unit must be specified and must be one of 'CPM', 'FPKM', 'FPK' or 'TPM'.")

    do.call("require", list("edgeR"))

    unit <- toupper(unit)
    if (unit %in% c('FPKM', 'TPM', 'FPK')) {
        # In these cases geneLength is required
        assertthat::assert_that(!missing(geneLength),
                                msg = "geneLength must be specified when unit is 'FPK', 'FPKM', or 'TPM.'")

        if ("matrix" %in% class(geneLength)) {# Flatten to a vector
            geneLength <- rowMeans(geneLength, na.rm = TRUE)
        }
    }
    # Make normalize method case insensitive (calcNormFactors is case sensitive)
    if (!is.null(normalize) &&
        any(!(is.character(normalize) || is.logical(normalize)),
            length(normalize) != 1,
            is.character(normalize) && length(normalize) == 1 && !tolower(normalize) %in% c('tmm', 'rle', 'upperquartile', 'tmmwzp', 'none'))) {
        warning("normalize must be only one of the following values 'TMM', 'RLE', 'upperquartile', 'TMMwzp', 'none', TRUE, FALSE or NULL. Assigning default values 'none'")
        normalize = "none"
    }

    if (is.null(normalize)) {
        normalize = "none"
    }

    if (toupper(normalize) %in% c("TMM", "RLE")) {
        normalize <- toupper(normalize)
    }
    if (toupper(normalize) %in% c("UPPERQUARTILE", "NONE")) {
        normalize <- tolower(normalize)
    }
    if (toupper(normalize) %in% c("TMMWZP")) {
        normalize <- "TMMwzp"
    }

    # Coerce countsMatrix to a matrix
    result <- countsMatrix <- as.matrix(countsMatrix)
    assertthat::assert_that("matrix" %in% class(countsMatrix),
                            msg = "countsMatrix must be able to be coerced to a matrix.")

    # Make sure geneLength is correct length
    if (!missing(geneLength)) {
        assertthat::assert_that(length(geneLength) == nrow(countsMatrix),
                                msg = "geneLength must be the same length of the number of rows in countsMatrix.")
    }

    # Set defaults
    if (any(is.null(log),
            !is.logical(log),
            length(log) != 1)) {
        warning("log must be a singular logical value. Assigning default value FALSE")
        log = FALSE
    }

    if (is.logical(normalize)) { # Don't encourage logicals; here for backward compatibility
        if (normalize == TRUE) {
            normalize <- 'TMM'
        }
        if (normalize == FALSE) {
            normalize <- 'none'
        }
    }

    if (is.null(prior.count)) {
        if (log == FALSE) {
            prior.count <- 0 # Not used when log = F
        }
        else if (unit == "TPM") {
            prior.count <- 0
        }
        else {
            prior.count <- 0.25
        }
    }

    result <- switch(toupper(unit),
                     "CPM"  = calcCPM(countsMatrix,  log, normalize, prior.count),
                     "FPKM" = calcFPKM(countsMatrix, log, normalize, geneLength, prior.count),
                     "FPK"  = calcFPK(countsMatrix,  log, normalize, geneLength, prior.count),
                     "TPM"  = calcTPM(countsMatrix,  log, normalize, geneLength, prior.count)
    )
    return(result)
}


#' Calculate TPM for a subsetted DGEobj
#'
#' Calculates TPM for a heavily subsetted DGEobj. The function will calculate TPM
#' using the original data but returns a DGEobj with the subset.
#'
#' TPM should be calculated on a full dataset with only low signal genes removed.
#' tpm.on.subset therefore allows calculation of TPM after heavy filtering of a DGEobj.
#'
#' Internally, convertCounts uses edgeR's fpkm() to calculate FPKM and converts to TPM
#' using the formula provided by [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
#'
#' @param dgeObj A DGEobj data structure
#' @param applyFilter Default = TRUE. If TRUE, reduces to the filtered gene list. FALSE returns
#'   all genes in the raw data.
#'
#' @return A matrix of TPM values
#'
#' @examples
#' \dontrun{
#'    # NOTE: Requires the edgeR package
#'
#'    dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'    tpm    <- tpm.on.subset(dgeObj)
#' }
#'
#' @import DGEobj
#' @importFrom assertthat assert_that
#' @export
tpm.on.subset <- function(dgeObj, applyFilter = TRUE){
    assertthat::assert_that(requireNamespace("edgeR", quietly = TRUE),
                            msg = "edgeR package is required to calculate FPKM and convert to TPM")

    assertthat::assert_that(!missing(dgeObj),
                            !is.null(dgeObj),
                            "DGEobj" %in% class(dgeObj),
                            msg = "dgeObj must be specified and should be of class 'DGEobj'.")
    assertthat::assert_that(attr(dgeObj, "level") %in% c("isoform", "gene"),
                            msg = "The level of dgeObj should be of type 'isoform' or type 'gene'.")

    do.call("require", list("edgeR"))


    if (any(is.null(applyFilter),
            !is.logical(applyFilter),
            length(applyFilter) != 1)) {
        warning("applyFilter must be a singular logical value. Assigning default value TRUE.")
        applyFilter = TRUE
    }

    # Default to gene level
    level <- "gene"
    if (attr(dgeObj, "level") == "isoform") {
        level <- "isoform"
    }

    if (level == "gene") {
        rowdata <- getItem(dgeObj, "geneData_orig")
    } else {
        rowdata <- getItem(dgeObj, "isoformData_orig")
    }

    # Get geneLength depending on source data
    if (attr(dgeObj, "source") == "Omicsoft") { # Omicsoft data
        geneLength <- rowdata$ExonLength
    } else if ("effectiveLength_orig" %in% names(dgeObj)) { # Use rowMeans(effectiveLength)
        geneLength <- rowMeans(getItem(dgeObj, "effectiveLength_orig"), na.rm = TRUE)
    }

    TPM <- convertCounts(getItem(dgeObj, "counts_orig"),
                         geneLength = geneLength,
                         unit = "tpm",
                         log = FALSE,
                         normalize = FALSE)

    # Remove filtered out genes
    if (applyFilter == TRUE) {
        idx <- rownames(TPM) %in% rownames(getItem(dgeObj, "counts"))
        TPM <- TPM[idx,]
    }
    return(TPM)
}


#' Convert countsMatrix and geneLength to TPM units
#'
#' Takes a countsMatrix and geneLength as input and converts to TPM units using the equation from
#' \href{https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/}{Harold Pimental}.
#'
#' The result should be the same as using convertCounts with normalize = 'tpm' and log = FALSE.
#'
#' geneLength can be a vector (length == nrow(countsMatrix)) or a matrix (same dim as countsMatrix).
#' The geneLength is used as is, or optionally collapsed to a vector by rowMeans.
#'
#' @param countsMatrix A numeric matrix of N genes x M samples. All columns must be numeric.
#' @param geneLength Numeric matrix of gene lengths. Often the ExonLength item of a DGEobj.
#' @param collapse Default = FALSE. TRUE or FALSE determines whether to use rowMeans on the geneLength matrix.
#'
#' @return A matrix of TPM values
#'
#' @examples
#'    dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'
#'    counts <- DGEobj::getItem(dgeObj, "counts")
#'    exonLength <- dgeObj$geneData$ExonLength
#'    tpm <- tpm.direct(counts, geneLength = exonLength)
#'
#' @importFrom assertthat assert_that
#'
#' @export
tpm.direct <- function(countsMatrix,
                       geneLength,
                       collapse = FALSE) {
    assertthat::assert_that(!missing(countsMatrix),
                            !is.null(countsMatrix),
                            "matrix" %in% class(countsMatrix) || "matrix" %in% class(as.matrix(countsMatrix)),
                            length(countsMatrix) != 0,
                            is.numeric(countsMatrix),
                            msg = "countsMatrix must be a numeric matrix of N genes x M Samples. All columns must be numeric.")
    if (!"matrix" %in% class(countsMatrix)) {
        countsMatrix <- as.matrix(countsMatrix)
    }
    assertthat::assert_that(!missing(geneLength),
                            !is.null(geneLength),
                            msg = "geneLength must be specified")
    if (is.vector(geneLength)) {
        assertthat::assert_that(length(geneLength) == nrow(countsMatrix),
                                msg = "geneLength should be of the same length as the number of rows in countsMatrix.")
    } else {
        if (!is.matrix(geneLength)) {
            result <- geneLength <- as.matrix(geneLength)
            assertthat::assert_that("matrix" %in% class(result),
                                    msg = "geneLength must be able to be coerced to a matrix.")
            assertthat::assert_that(all(dim(countsMatrix) == dim(geneLength)),
                                    msg = "The dimensions of countsMatrix and geneLength should match.")
        }
    }

    if (any(is.null(collapse),
            !is.logical(collapse),
            length(collapse) != 1)) {
        warning("collapse must be a singular logical value. Assigning default value FALSE.")
        collapse = FALSE
    }

    if (collapse & is.matrix(geneLength)) {
        geneLength <- rowMeans(geneLength, na.rm = TRUE)
    }

    # Calculation  (fpk / colsum(fpk) ) * 10e6
    fpb <- countsMatrix / geneLength
    sumfpb <- colSums(fpb)
    tpm <- fpb / expandAsMatrix(sumfpb, byrow = TRUE, dim = dim(fpb)) * 1e6
}

# Helper Functions
calcCPM <- function(countsMatrix, log, normalize, prior.count){
    tryCatch({
        do.call("cpm",
                list(y      = do.call("calcNormFactors",
                                      list(object = do.call("DGEList",
                                                            list(counts = countsMatrix)),
                                           method = normalize)),
                     log         = log,
                     prior.count = prior.count))
    },
    error = function(e) {
        message("Unexpected error: ", e$message, " happened during edgeR computation of CPM")
        return(NULL)
    })
}

calcFPKM <- function(countsMatrix, log, normalize, geneLength, prior.count){
    tryCatch({
        do.call("rpkm",
                list(y      = do.call("calcNormFactors",
                                      list(object = do.call("DGEList",
                                                            list(counts = countsMatrix)),
                                           method = normalize)),
                     log         = log,
                     gene.length = geneLength,
                     prior.count = prior.count))
    },
    error = function(e) {
        message("Unexpected error: ", e$message, " happened during edgeR computatin of RPKM")
        return(NULL)
    })
}

calcTPM <- function(countsMatrix, log, normalize, geneLength, prior.count){
    if (normalize != "none") {
        warning(paste('TPM normalization overides', normalize, 'normalization!'))
    }
    if (prior.count != 0 && log == TRUE) {
        warning("Using a prior.count for logTPM calculations is not recommended and may produce unpredictable results!")
    }

    fpkm <- calcFPKM(countsMatrix, log = log, normalize = normalize,
                     geneLength = geneLength, prior.count = prior.count)

    # Helper function
    fpkmToTpm <- function(fpkm) {
        colSumMat <- expandAsMatrix(colSums(fpkm, na.rm = TRUE), byrow = TRUE, dim = dim(fpkm))
        fpkm / colSumMat * 1e6
    }

    if (log == FALSE) {
        TPM <- fpkmToTpm(fpkm)
    } else {
        TPM <- log2(fpkmToTpm(2^fpkm))
    }
    return(TPM)
}

calcFPK <- function(countsMatrix, log, normalize, geneLength, prior.count){
    if (tolower(normalize) == 'none') {
        # Check for zero geneLength just in case
        if (min(geneLength) == 0) {
            geneLength <- geneLength + 1
        }
    }
    FPK <- countsMatrix / (geneLength / 1000)
    if (log == TRUE) {
        FPK <- log2(FPK + prior.count)
    }
    return(FPK)
}

# brought in from edgeR
#' @importFrom methods is
expandAsMatrix <- function(x, dim = NULL, byrow = TRUE) {
    if (is.null(dim))
        return(as.matrix(x))
    if (length(dim) < 2)
        stop("dim must be numeric vector of length 2")
    dim <- round(dim[1:2])
    if (any(dim < 0))
        stop("negative dimensions not allowed")
    dx <- dim(x)
    if (is.null(dx)) {
        lx <- length(x)
        if (lx == 1)
            return(matrix(x, dim[1], dim[2]))
        if (lx == dim[1] & lx == dim[2])
            return(matrix(x, dim[1], dim[2], byrow = byrow))
        if (lx == dim[2])
            return(matrix(x, dim[1], dim[2], byrow = TRUE))
        if (lx == dim[1])
            return(matrix(x, dim[1], dim[2], byrow = FALSE))
        stop("x of unexpected length")
    }
    if (length(dx) < 2)
        stop("x has less than 2 dimensions")
    if (length(dx) > 2)
        stop("x has more than 2 dimensions")
    if (all(dx == dim))
        return(as.matrix(x))
    if (methods::is(x, "CompressedMatrix")) {
        return(Recall(as.matrix(x), dim = dim, byrow = byrow))
    }
    stop("x is matrix of wrong size")
}
