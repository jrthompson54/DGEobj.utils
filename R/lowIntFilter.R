#' Apply low intensity filters to a DGEobj
#'
#' Takes a DGEobj as input and applies a combination of low intensity filters as
#' specified by the user. Raw count, zFPKM, TPM, and/or FPK filters are
#' supported.  A gene must pass all active filters.  Not setting a threshold
#' argument inactivates that threshold.
#'
#' @param dgeObj A DGEobj with RNA-Seq (counts) data (Required)
#' @param countThreshold Genes below this threshold are removed (10 is recommended).
#' @param zfpkmThreshold Genes below this threshold are removed. (-3.0 is recommended)
#' @param fpkThreshold Genes below this threshold are removed. (5 is recommended)
#' @param tpmThreshold Genes below this threshold are removed.
#' @param sampleFraction The proportion of samples that must meet the thresholds
#'   (Default = 0.5). Range > 0 and <= 1.
#' @param geneLength Vector of geneLengths for rows of dgeObj. Required for FPK and
#'   zFPKM filters, unless dgeObj is a DGEobj.  If a DGEobj is supplied, geneLength is
#'   retrieved from the DGEobj, unless supplied by the geneLength argument.
#' @param verbose Prints messages about the filtering process.
#'
#' @return Same class as input object with low intensity rows removed
#'
#' @examples
#'   myDGEobj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'   dim(myDGEobj)
#'
#'   # Simple count threshold in at least 3/4ths the samples
#'   myDGEobj <- lowIntFilter(myDGEobj,
#'                            countThreshold = 10,
#'                            sampleFraction = 0.5)
#'   dim(myDGEobj)
#'
#'   # Count and FPK thresholds
#'   myDGEobj <- lowIntFilter(myDGEobj,
#'                            countThreshold = 10,
#'                            fpkThreshold = 5,
#'                            sampleFraction = 0.5)
#'   dim(myDGEobj)
#'
#' @importFrom assertthat assert_that
#' @importFrom DGEobj getItem
#' @importFrom zFPKM zFPKM
#' @importFrom stringr str_c
#' @importFrom stats complete.cases
#'
#' @export
lowIntFilter <- function(dgeObj,
                         countThreshold,
                         zfpkmThreshold,
                         fpkThreshold,
                         tpmThreshold,
                         sampleFraction = 0.5,
                         geneLength,
                         verbose = FALSE)
{
    assertthat::assert_that(!missing(dgeObj),
                            !is.null(dgeObj),
                            "DGEobj" %in% class(dgeObj),
                            msg = "dgeObj must be of class 'DGEobj'.")

    if (!missing(zfpkmThreshold) && !missing(tpmThreshold)) {
        assertthat::assert_that(is.null(zfpkmThreshold) & !is.null(tpmThreshold),
                                is.null(tpmThreshold) & !is.null(zfpkmThreshold),
                                msg = "Must use zfpkmThreshold or tpmThreshold, but not both.")
    }

    if (any(is.null(sampleFraction),
            !is.numeric(sampleFraction),
            length(sampleFraction) != 1)) {
        warning("sampleFraction must be a singular numeic value. Assigning default value 0.5")
        sampleFraction = 0.5
    }

    if (any(is.null(verbose),
            !is.logical(verbose),
            length(verbose) != 1)) {
        warning("verbose must be a singular logical value. Assigning default value FALSE")
        verbose = FALSE
    }
    counts <- getItem(dgeObj, "counts")
    starting_rowcount <- nrow(counts)

    # Get geneLength from DGEobj
    if (missing(geneLength) || is.null(geneLength)) {  # User supplied geneLength supercedes geneLength from DGEobj
        geneLength <- dgeObj$geneData$ExonLength
    }

    # Apply zFPKM threshold
    if (!missing(zfpkmThreshold)) {
        assertthat::assert_that(!is.null(geneLength),
                                length(geneLength) == nrow(dgeObj$counts),
                                msg = "geneLength must be specified and should be the same length as the number of rows in dgeObj's counts.")
        fpkm <- convertCounts(dgeObj$counts, unit = "fpkm", geneLength = geneLength)

        # Need to filter out rows filled with NaNs first before calculating zFPKM
        idx <- complete.cases(fpkm)
        dgeObj <- dgeObj[row = idx,]
        fpkm <- fpkm[idx,]
        geneLength <- geneLength[idx]
        nan_genes <- sum(!idx)
        if (nan_genes > 0 & verbose == TRUE) {
            message(stringr::str_c(nan_genes, " genes with NaN FPKM values removed."))
        }

        # Calculate zFPKM
        zfpkm <- as.matrix(zFPKM::zFPKM(as.data.frame(fpkm)))

        # Create index for zFPKM >= zFPKMThreshold in fracThreshold of samples
        idx_zfpkm <- zfpkm >= zfpkmThreshold
        frac <- rowSums(idx_zfpkm) / ncol(idx_zfpkm)
        fpkmidx <- frac >= sampleFraction
        dgeObj <- subset(dgeObj, row = fpkmidx)
        geneLength <- geneLength[fpkmidx]

        if (verbose == TRUE) {
            message(stringr::str_c(sum(fpkmidx), " of ", starting_rowcount, " genes retained by the zFPKM filter."))
        }
    }

    # Apply TPM threshold
    if (!missing(tpmThreshold)) {
        assertthat::assert_that(!is.null(geneLength),
                                length(geneLength) == nrow(counts),
                                msg = "geneLength must be specified and should be the same length as the number of rows in counts.")
        tpm <- convertCounts(dgeObj$counts, unit = "tpm", geneLength = geneLength)

        # Need to filter out rows filled with NaNs first before calculating zFPKM
        idx <- complete.cases(tpm)
        dgeObj <- subset(dgeObj, row = idx)
        tpm <- tpm[idx,]
        geneLength <- geneLength[idx]
        nan_genes <- sum(!idx)
        if (nan_genes > 0 & verbose == TRUE) {
            message(stringr::str_c(nan_genes, " genes with NaN TPM values removed."))
        }

        #create index for tpm >=tpmThreshold in fracThreshold of samples
        idx_tpm <- tpm >= tpmThreshold
        frac <- rowSums(idx_tpm) / ncol(idx_tpm)
        tpmidx <- frac >= sampleFraction
        dgeObj <- subset(dgeObj, row = tpmidx)
        geneLength <- geneLength[tpmidx]

        if (verbose == TRUE) {
            message(stringr::str_c(sum(tpmidx), " of ", starting_rowcount, " genes retained by the TPM filter."))
        }
    }

    # Apply FPK threshold
    if (!missing(fpkThreshold)) {
        assertthat::assert_that(!is.null(geneLength),
                                length(geneLength) == nrow(dgeObj),
                                msg = "geneLength must be specified and should be the same length as the number of rows in dgeObj.")
        fpk <- convertCounts(dgeObj$counts, unit = "fpk", geneLength = geneLength)
        # Keep FPK >= fpkThreshold in fracThreshold of samples
        idxfpk <- fpk >= fpkThreshold
        frac <- rowSums(idxfpk)/ncol(idxfpk)
        idx <- frac >= sampleFraction
        dgeObj <- subset(dgeObj, row = idx)
        geneLength <- geneLength[idx]

        if (verbose == TRUE) {
            message(stringr::str_c(sum(idx), " of ", length(idx), " genes retained by the FPK filter."))
        }
    }

    # Apply count threshold
    if (!missing(countThreshold)) {
        # Overlay a min count filter
        idxmin <- dgeObj$counts >= countThreshold
        frac <- rowSums(idxmin) / ncol(idxmin)
        idx <- frac >= sampleFraction
        dgeObj <- subset(dgeObj, row = idx)
        geneLength <- geneLength[idx]

        if (verbose == TRUE) {
            message(stringr::str_c(sum(idx), " of ", length(idx), " genes retained by the low count filter."))
        }
    }

    return(dgeObj)
}
