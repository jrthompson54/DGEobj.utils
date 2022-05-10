#' Run functions in a typical voom/lmFit workflow
#'
#' In the default workflow, this function runs voomWithQualityWeights followed by
#' lmFit and optionally eBayes. If the contrasts of interest are already represented
#' in the model, enable eBayes. To use contrasts.fit downstream, run eBayes
#' after that step instead. eBayes should always be run last.
#'
#' Input is minimally a DGEobj containing a DGEList (typically TMM-normalized),
#' and a formula (character representation).  Other arguments can invoke
#' the duplicateCorrelation method and modify use of quality weights.
#'
#' Returns a DGEobj class object containing the VoomElist (voom
#' output), and Fit object (lmFit output).
#'
#' Quality weights should be enabled unless there is a good reason to turn them
#' off. If all samples are equal quality, the weights will all approach 1.0 with no
#' consequence on the results. More typically, some samples are better than others
#' and using quality weights improves the overall result.
#'
#' Use var.design if the quality weights are correlated with some factor in the experiment.
#' This will cause the quality weights to be calculated as a group instead of individually.
#'
#' Use duplicate correlation (dupCorBlock) when there are subjects that have been sampled more
#' than once (e.g. before and after some treatment).  This calculates a within-
#' subject correlation and includes this in the model.
#'
#' @param dgeObj A DGEobj containing a DGEList (e.g. from runEdgeRNorm.) (Required)
#' @param designMatrixName Name of a design matrix within dgeObj. (Required)
#' @param dupCorBlock Supply a block argument to trigger duplicateCorrelation. (Optional)
#'    Should be a vector the same length as ncol with values to indicate common
#'    group membership for duplicateCorrelation.
#'    Also, 'statmod' package must be installed to run duplicate correlation calculations.
#' @param runDupCorTwice Default = TRUE. Gordon Smyth recommends running duplicateCorrelation
#'   twice. Set this to false to run duplicateCorrelation just once.
#' @param qualityWeights Runs limma::voomWithQualityWeights if set to TRUE (Default = TRUE).
#'    This should normally be set to TRUE.
#' @param var.design Provide a design matrix (from model.matrix) to identify
#'    replicate groups (e.g. "~ ReplicateGroup") for quality weight determination.
#'    Causes quality weights to be determined on a group basis.  If omitted
#'    limma::voomWithQualityWeights treats each sample individually.
#' @param runEBayes Runs eBayes after lmFit. (Default = TRUE)
#'    Note, 'statmod' package must be installed to run eBayes calculations.
#' @param robust Used by eBayes. (Default = TRUE)
#'    Note, 'statmod' package must be installed to run eBayes calculations.
#' @param proportion Proportion of genes expected to be differentially expressed
#'   (used by eBayes) (Default = 0.01) Modify the prior accordingly if the 1st pass analysis shows
#'   a significantly higher or lower proportion of genes regulated than the default.
#' @param mvPlot Enables the voom mean-variance plot. (Default = TRUE)
#'
#' @return A DGEobj now containing designMatrix, Elist, and fit object.
#'
#' @examples
#' \dontrun{
#'    dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'    for (name in names(dgeObj)[11:length(dgeObj)]) {
#'        dgeObj <- DGEobj::rmItem(dgeObj, name)
#'    }
#'
#'    dgeObj <- runVoom(dgeObj,
#'                      designMatrixName = "ReplicateGroupDesign",
#'                      mvPlot = TRUE)
#'
#'    # Note the Elist and fit objects have been added
#'    DGEobj::inventory(dgeObj)
#'}
#'
#' @importFrom limma voom lmFit eBayes voomWithQualityWeights duplicateCorrelation
#' @importFrom stringr str_c
#' @importFrom DGEobj getItem addItem
#' @importFrom dplyr %>%
#' @importFrom assertthat assert_that
#'
#' @export
runVoom <- function(dgeObj,
                    designMatrixName,
                    dupCorBlock,
                    runDupCorTwice = TRUE,
                    qualityWeights = TRUE,
                    var.design,
                    mvPlot = TRUE,
                    runEBayes = TRUE,
                    robust = TRUE,
                    proportion = 0.01) {
    assertthat::assert_that(!missing(dgeObj),
                            !is.null(dgeObj),
                            "DGEobj" %in% class(dgeObj),
                            msg = "dgeObj must be specified and must be of class 'DGEobj'.")
    assertthat::assert_that(!missing(designMatrixName),
                            !is.null(designMatrixName),
                            is.character(designMatrixName),
                            length(designMatrixName) == 1,
                            designMatrixName %in% names(dgeObj),
                            msg = "designMatrixName must be specified and must be one of the items in dgeObj. Use names(dgeObj) to check for available options.")
    assertthat::assert_that("DGEList" %in% DGEobj::showTypes(dgeObj)$Type,
                            msg = "No DGEList found in dgeObj. Specify a DGEobj that contains a DGEList.")
    designMatrix <- DGEobj::getItem(dgeObj, designMatrixName)

    if ("DGEList" %in% attr(dgeObj, "type")) {
        dgelist <- DGEobj::getItem(dgeObj, "DGEList")
    }

    if (any(is.null(runDupCorTwice),
            !is.logical(runDupCorTwice),
            length(runDupCorTwice) != 1)) {
        warning("runDupCorTwice must be a singular logical value. Assigning default value TRUE")
        runDupCorTwice = TRUE
    }

    if (any(is.null(qualityWeights),
            !is.logical(qualityWeights),
            length(qualityWeights) != 1)) {
        warning("qualityWeights must be a singular logical value. Assigning default value TRUE")
        qualityWeights = TRUE
    }

    if (any(is.null(mvPlot),
            !is.logical(mvPlot),
            length(mvPlot) != 1)) {
        warning("mvPlot must be a singular logical value. Assigning default value TRUE")
        mvPlot = TRUE
    }

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

    if (any(is.null(proportion),
            !is.numeric(proportion),
            length(proportion) != 1)) {
        warning("proportion must be a singular numeric value. Assigning default value 0.01")
        proportion = 0.01
    }

    # Collect calling args for documentation
    funArgs <- match.call()

    # Set run parameters
    dupcor <- FALSE
    if (!missing(dupCorBlock)) {
        assertthat::assert_that(requireNamespace("statmod", quietly = TRUE),
                                msg = "'statmod' package is required to run duplicate correlation calculations")
        dupcor <- TRUE
    }

    if (robust) {
        assertthat::assert_that(requireNamespace("statmod", quietly = TRUE),
                                msg = "'statmod' package is required to run eBayes calculations")
    }

    blockQW <- FALSE
    if (qualityWeights && !missing(var.design)) {
        blockQW <- TRUE
    }

    # Main Calculations (one of six blocks will be run)
    tryCatch({
        # Set type of analysis
        if (!dupcor && !qualityWeights && !blockQW) {
            # Voom squeezes the variance (borrowing from other genes) to deal
            # with the heteroskedasticity problem
            VoomElist <- limma::voom(dgelist, designMatrix, plot = mvPlot)
            fit <- limma::lmFit(VoomElist, designMatrix)
        } else if (!dupcor && qualityWeights && !blockQW) { # indQW analysis
            VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix, plot = mvPlot, col = "blue")
            fit <- limma::lmFit(VoomElist, designMatrix)
        } else if (!dupcor && qualityWeights && blockQW) { # blockedQW analysis
            VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix, plot = mvPlot, col = "blue", var.design = var.design)
            fit <- limma::lmFit(VoomElist, designMatrix)
        } else if (dupcor && !qualityWeights && !blockQW) { # dupcor_base analysis
            VoomElist <- limma::voom(dgelist, designMatrix)
            corfit <- limma::duplicateCorrelation(VoomElist,
                                                  designMatrix,
                                                  block = dupCorBlock)
            if (runDupCorTwice) {
                VoomElist <- limma::voom(dgelist, designMatrix,
                                         correlation = corfit$consensus.correlation,
                                         plot = mvPlot)
                corfit <- limma::duplicateCorrelation(VoomElist,
                                                      designMatrix,
                                                      block = dupCorBlock)
            }
            fit <- limma::lmFit(VoomElist, designMatrix, block = dupCorBlock,
                                correlation = corfit$consensus.correlation)
        } else if (dupcor && qualityWeights && !blockQW) { # dupcor_indQW analysis
            VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix)
            corfit <- limma::duplicateCorrelation(VoomElist,
                                                  designMatrix,
                                                  block = dupCorBlock)
            if (runDupCorTwice) {
                VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix,
                                                           plot = mvPlot, col = "blue",
                                                           correlation = corfit$consensus.correlation)
                corfit <- limma::duplicateCorrelation(VoomElist,
                                                      designMatrix,
                                                      block = dupCorBlock)
            }
            fit <- limma::lmFit(VoomElist, designMatrix, block = dupCorBlock,
                                correlation = corfit$consensus.correlation)
        } else if (dupcor && qualityWeights && blockQW) { # dupcor_vdQW analysis
            VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix,
                                                       var.design = var.design)
            corfit <- limma::duplicateCorrelation(VoomElist,
                                                  designMatrix,
                                                  block = dupCorBlock)
            if (runDupCorTwice) {
                VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix,
                                                           plot = mvPlot, col = "blue",
                                                           correlation = corfit$consensus.correlation,
                                                           var.design = var.design)
                corfit <- limma::duplicateCorrelation(VoomElist,
                                                      designMatrix,
                                                      block = dupCorBlock)
            }
            fit <- limma::lmFit(VoomElist, designMatrix, block = dupCorBlock,
                                correlation = corfit$consensus.correlation)
        }

        # Run eBayes
        if (runEBayes) {
            fit = limma::eBayes(fit, robust = robust, proportion = proportion)
            itemAttr <- list(eBayes = TRUE)
        } else {
            itemAttr <- list(eBayes = FALSE)
        }

        if (exists("corfit")) { # Duplicate correlation was used; capture the correlation value
            cat(stringr::str_c("Duplicate Correlation = ", round(corfit$consensus.correlation, 4), "   \n"))
            attr(VoomElist, "DupCor") <- corfit$consensus.correlation
            attr(fit, "DupCor") <- corfit$consensus.correlation
        }

        VoomElistName = paste(designMatrixName, "_Elist", sep = "")
        dgeObj <- dgeObj %>%
            DGEobj::addItem(VoomElist, VoomElistName,
                            "Elist",
                            funArgs = funArgs,
                            parent = list("DGEList", designMatrixName))
        # Add corfit if present
        if (exists("corfit")) {
            dgeObj <- dgeObj %>%
                DGEobj::addItem(corfit, paste(designMatrixName, "_corFit", sep = ""),
                                "corFit",
                                funArgs = funArgs,
                                parent = paste(designMatrixName, "_Elist", sep = ""))
        }
        dgeObj <- dgeObj %>%
            DGEobj::addItem(fit, paste(designMatrixName, "_fit", sep = ""),
                            "fit",
                            funArgs = funArgs,
                            itemAttr = itemAttr,
                            parent = list(VoomElistName, designMatrixName))
        },
        error = function(e) {
            message(paste("runVoom did not complete successfully due to: ", e))
        })
    dgeObj
}
