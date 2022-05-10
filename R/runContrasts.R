#' Build contrast matrix and calculate contrast fits
#'
#' Takes a DGEobj and a named list of contrasts to build. The DGEobj must
#' contain a limma Fit object and associated designMatrix. Returns the DGEobj with
#' contrast fit(s), contrast matrix, and topTable/topTreat dataframes added.
#'
#' The contrastList is a named list.  The values are composed of column names
#' from the designMatrix of the DGEobj.  Each contrast is named to give it a
#' short, recognizable name to be used for display purposes.
#'
#' Example contrastList \cr
#'
#' contrastList = list( \cr
#'    T1 = "treatment1 - control", \cr
#'    T2 = "treatment2 - control" \cr
#' ) \cr
#'
#' where treatment1, treatment2, and control are column names in the designMatrix.
#'
#' The returned DGEobj list contains the new items:
#' \itemize{
#'     \item{"contrastMatrix"} {a matrix}
#'     \item{"Fit.Contrasts"} {a Fit object}
#'     \item{"topTableList"} {a List of dataframes}
#'     \item{"topTreatList"} {a List of dataframes: if enabled}
#' }
#'
#' @param dgeObj A DGEobj object containing a Fit object and design matrix. (Required)
#' @param designMatrixName The name of the design matrix within dgeObj to use for
#'    contrast analysis. (Required)
#' @param contrastList A named list of contrasts. (Required)
#' @param contrastSetName Name for the set of contrasts specified in
#'   contrastList.  Defaults to "<fitName>_cf". Only needed to create 2 or more
#'   contrast sets from the same initial fit.
#' @param runTopTable Runs topTable on the specified contrasts. (Default = TRUE)
#' @param runTopTreat Runs topTreat on the specified contrasts. (Default = FALSE)
#' @param foldChangeThreshold Only applies to topTreat (Default = 1.5)
#' @param runEBayes Runs eBayes after contrast.fit (Default = TRUE)
#' @param robust eBayes robust option. 'statmod' package must be installed to perform required
#' calculations if robust = TRUE (Default = TRUE)
#' @param proportion Proportion of genes expected to be differentially expressed.
#'   (used by eBayes) (Default = 0.01)
#' @param qValue Set TRUE to include Q-values in topTable output. (Default = FALSE)
#' @param IHW Set TRUE to add FDR values from the IHW package. (Default = FALSE)
#' @param verbose Set TRUE to print some information during processing. (Default = FALSE)
#'
#' @return The DGEobj with contrast matrix, fit and topTable/topTreat dataframes added.
#'
#' @examples
#' \dontrun{
#'    myDGEobj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'
#'    # Name the design matrix to be used (see inventory(myDGEobj))
#'    designMatrixName <- "ReplicateGroupDesign"
#'
#'    # Define the named contrasts from design matrix column names
#'    contrastList  <- list(BDL_v_Sham = "ReplicateGroupBDL - ReplicateGroupSham",
#'                          EXT1024_v_BDL = "ReplicateGroupBDL_EXT.1024  - ReplicateGroupBDL",
#'                          Nint_v_BDL = "ReplicateGroupBDL_Nint - ReplicateGroupBDL",
#'                          Sora_v_BDL = "ReplicateGroupBDL_Sora - ReplicateGroupBDL")
#'
#'
#'    myDGEobj <- runContrasts(myDGEobj,
#'                             designMatrixName=designMatrixName,
#'                             contrastList=contrastList,
#'                             contrastSetName = "SecondContrastSet",
#'                             qValue = TRUE,
#'                             IHW = FALSE,
#'                             runTopTable = TRUE,
#'                             runTopTreat = TRUE,
#'                             foldChangeThreshold = 1.25)
#'    DGEobj::inventory(myDGEobj)
#'}
#'
#' @importFrom limma contrasts.fit eBayes makeContrasts topTable topTreat treat
#' @importFrom DGEobj addItem getItem
#' @importFrom assertthat assert_that
#' @importFrom stringr str_c
#'
#' @export
runContrasts <- function(dgeObj,
                         designMatrixName,
                         contrastList,
                         contrastSetName = NULL,
                         runTopTable = TRUE,
                         runTopTreat = FALSE,
                         foldChangeThreshold = 1.5,
                         runEBayes = TRUE,
                         robust = TRUE,
                         proportion = 0.01,
                         qValue = FALSE,
                         IHW = FALSE,
                         verbose = FALSE) {
    assertthat::assert_that(!missing(dgeObj),
                            !is.null(dgeObj),
                            "DGEobj" %in% class(dgeObj),
                            msg = "dgeObj must be specified and should be of class 'DGEobj'.")
    assertthat::assert_that(!missing(designMatrixName),
                            !is.null(designMatrixName),
                            is.character(designMatrixName),
                            length(designMatrixName) == 1,
                            designMatrixName %in% names(dgeObj),
                            msg = "designMatrixName must be a signular character value and one of dgeobj names.")
    assertthat::assert_that("list" %in% class(contrastList),
                            !missing(contrastList),
                            !is.null(names(contrastList)),
                            msg = "contrastList must specified and must be a named list.")
    assertthat::assert_that(foldChangeThreshold >= 0,
                            msg = "foldChangeThreshold must be greater than or equal to 0.")
    assertthat::assert_that(!(runTopTable == FALSE & runTopTreat == FALSE),
                            msg = "One of runTopTable or runTopTreat must be TRUE.")
    assertthat::assert_that(designMatrixName %in% names(dgeObj),
                            msg = "The specified designMatrixName not found in dgeObj.")

    fitName <- paste(designMatrixName, "_fit", sep = "")
    assertthat::assert_that(fitName %in% names(dgeObj),
                            msg = "The specified fitName object not found in dgeObj.")

    if (any(is.null(runEBayes),
            !is.logical(runEBayes),
            length(runEBayes) != 1)) {
        warning("runEBayes must be a singular logical value. Assigning default value TRUE")
        runEBayes = TRUE
    }

    if (any(is.null(robust),
            !is.logical(robust),
            length(robust) != 1)) {
        warning("robust must be a singular logical value. Assigning default value TRUE")
        robust = TRUE
    }

    if (robust) {
        assertthat::assert_that(requireNamespace("statmod", quietly = TRUE),
                                msg = "'statmod' package is required to run eBayes calculations")
    }

    if (any(is.null(proportion),
            !is.numeric(proportion),
            length(proportion) != 1)) {
        warning("proportion must be a singular numeric value. Assigning default value 0.01")
        proportion = 0.01
    }

    if (any(is.null(qValue),
            !is.logical(qValue),
            length(qValue) != 1)) {
        warning("qValue must be a singular logical value. Assigning default value FALSE")
        qValue = FALSE
    }

    if (any(is.null(IHW),
            !is.logical(IHW),
            length(IHW) != 1)) {
        warning("IHW must be a singular logical value. Assigning default value FALSE")
        IHW = FALSE
    }

    if (any(is.null(verbose),
            !is.logical(verbose),
            length(verbose) != 1)) {
        warning("verbose must be a singular logical value. Assigning default value FALSE")
        verbose = FALSE
    }

    if (any(is.null(contrastSetName),
            !is.character(contrastSetName),
            length(contrastSetName) != 1)) {
        warning("contrastSetName must be a character value. Assigning default value: ", fitName)
            contrastSetName <- fitName
    }

    assertthat::assert_that(!(paste0(contrastSetName, '_cm') %in% names(dgeObj)),
                            msg = "The contrastSetName already exists in dgeObj.")

    funArgs <- match.call()

    # Retrieve designMatrix & fit
    designMatrix <- DGEobj::getItem(dgeObj, designMatrixName)
    fit <- DGEobj::getItem(dgeObj, fitName)

    # Run the contrast fit
    ContrastMatrix <- limma::makeContrasts(contrasts = contrastList, levels = designMatrix)
    MyFit.Contrasts <- limma::contrasts.fit(fit, ContrastMatrix)

    # Run eBayes
    if (runEBayes) {
        if (verbose) {
            .tsmsg(stringr::str_c("running EBayes: proportion = ", proportion))
        }
        MyFit.Contrasts = limma::eBayes(MyFit.Contrasts, robust = robust, proportion = proportion)
        MyFit.Contrasts.treat = limma::treat(MyFit.Contrasts, lfc = log2(foldChangeThreshold),
                                             robust = robust)
    }

    # Run topTable on each contrast and add each DF to a list
    if (runTopTable) {
        if (verbose) {
            .tsmsg("Running topTable...")
        }
        # Run topTable via lapply to generate a bunch of contrasts.
        MyCoef <- as.list(1:length(contrastList))
        topTableList <- lapply(MyCoef, function(x) (limma::topTable(MyFit.Contrasts, coef = x,
                                                                    confint = T, number = Inf, p.value = 1, sort.by = "none")))

        # Transfer the contrast names
        names(topTableList) = names(contrastList)

        if (qValue) {
            topTableList <- runQvalue(topTableList)
        }
        if (IHW) {
            IHW_result <- runIHW(topTableList)
            topTableList <- IHW_result[[1]]
        }
    }

    if (runTopTreat) {
        if (verbose) {
            .tsmsg("Running topTreat...")
        }
        # Run topTreat via lapply to generate a bunch of contrasts.
        LFC    <- log2(foldChangeThreshold)
        MyCoef <- as.list(1:length(contrastList))
        topTreatList = lapply(MyCoef, function(x) (limma::topTreat(MyFit.Contrasts.treat, coef = x,
                                                                   confint = T, lfc = LFC, number = Inf, p.value = 1, sort.by = "none")))
        # Transfer the contrast names
        names(topTreatList) = names(contrastList)
    }

    # Capture the contrast matrix
    dgeObj <- DGEobj::addItem(dgeObj,
                              item = ContrastMatrix,
                              itemName = paste(contrastSetName, "_cm", sep = ""),
                              itemType = "contrastMatrix",
                              funArgs = funArgs,
                              parent = fitName)

    if (runTopTable) {
        # Add the contrast fit
        dgeObj <- DGEobj::addItem(dgeObj,
                                  item = MyFit.Contrasts,
                                  itemName = paste(contrastSetName, "_cf", sep = ""),
                                  itemType = "contrast_fit",
                                  funArgs = funArgs,
                                  parent = fitName)

        # Add the topTable DFs
        listNames <- names(topTableList)
        for (i in 1:length(topTableList)) {
            dgeObj <- DGEobj::addItem(dgeObj,
                                      item = topTableList[[i]],
                                      itemName = listNames[[i]],
                                      itemType = "topTable",
                                      funArgs = funArgs,
                                      parent = paste(contrastSetName, "_cf", sep = ""))
        }
    }

    if (runTopTreat) {
        dgeObj <- DGEobj::addItem(dgeObj,
                                  item = MyFit.Contrasts.treat,
                                  itemName = paste(contrastSetName, "_cft", sep = ""),
                                  itemType = "contrast_fit_treat",
                                  funArgs = funArgs,
                                  parent = fitName)

        listNames <- names(topTreatList)
        for (i in 1:length(topTreatList)) {
            dgeObj <- DGEobj::addItem(dgeObj,
                                      item = topTreatList[[i]],
                                      itemName = paste(listNames[i], "_treat", sep = ""),
                                      itemType = "topTreat",
                                      funArgs = funArgs,
                                      parent = paste(contrastSetName, "_cft", sep = ""))
        }
    }
    dgeObj
}
