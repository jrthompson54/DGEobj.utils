#' Run functions in a typical voom/lmFit workflow
#'
#' In the default workflow, this function runs voomWithQualityWeights followed by
#' lmFit and optionally eBayes. If the contrasts of interest are already represented
#' in the model, enable eBayes. To use contrasts.fit downstream, run eBayes
#' after that step instead. eBayes should always be run last.
#'
#' Input is minimally a DGEobj containing a DGEList (which is typically TMM-normalized)
#' and a formula (character representation).  If a DGEList is missing on the object the
#' counts are used as-is.  Other arguments can invoke the duplicateCorrelation method and
#' modify use of quality weights.
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
#' @param dgeObj A DGEobj containing a DGEList (e.g. from runEdgeRNorm) or counts (Required)
#' @param designMatrixName Name of a design matrix within dgeObj. (Required)
#' @param dupCorBlock Supply a block argument to trigger duplicateCorrelation. (Optional)
#'    Should be a vector the same length as ncol with values to indicate common
#'    group membership for duplicateCorrelation.
#'    Also, 'statmod' package must be installed to run duplicate correlation calculations.
#' @param runDupCorTwice Default = TRUE. Gordon Smyth recommends running duplicateCorrelation
#'   twice. Set this to false to run duplicateCorrelation just once.
#' @param qualityWeights Runs limma's voomWithQualityWeights() if set to TRUE (Default = TRUE).
#'    This should normally be set to TRUE.
#' @param var.design Provide a design matrix (from model.matrix) to identify
#'    replicate groups (e.g. "~ ReplicateGroup") for quality weight determination.
#'    Causes quality weights to be determined on a group basis.  If omitted
#'    limma's voomWithQualityWeights() treats each sample individually.
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
#'   # NOTE: Requires the limma package
#'
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
    assertthat::assert_that(requireNamespace("limma", quietly = TRUE),
                            msg = "limma package is required to prepare data for linear modeling")

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
    assertthat::assert_that(any(c("DGEList", "counts") %in% attr(dgeObj, "type")),
                            msg = "No counts or DGEList found in dgeObj. Specify a DGEobj that contains a DGEList or counts.")

    do.call("require", list("limma"))

    designMatrix <- DGEobj::getItem(dgeObj, designMatrixName)

    if ("DGEList" %in% attr(dgeObj, "type")) {
        dgelist <- DGEobj::getItem(dgeObj, "DGEList")
    } else {
        dgelist <- DGEobj::getItem(dgeObj, "counts")
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
    do.call("require", list("limma"))

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
            VoomElist <- tryCatch({
                do.call("voom",
                        list(counts = dgelist,
                             design = designMatrix,
                             plot   = mvPlot))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during limma voom()")
                return(NULL)
            })

            fit <- tryCatch({
                do.call("lmFit",
                        list(object = VoomElist,
                             design = designMatrix))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during linear modeling")
                return(NULL)
            })
        } else if (!dupcor && qualityWeights && !blockQW) {
            # indQW analysis
            VoomElist <- tryCatch({
                do.call("voomWithQualityWeights",
                        list(counts = dgelist,
                             design = designMatrix,
                             plot   = mvPlot,
                             col    = 'blue'))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during limma voomWithQualityWeights()")
                return(NULL)
            })

            fit <- tryCatch({
                do.call("lmFit",
                        list(object = VoomElist,
                             design = designMatrix))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during linear modeling")
                return(NULL)
            })
        } else if (!dupcor && qualityWeights && blockQW) {
            # blockedQW analysis
            VoomElist <- tryCatch({
                do.call("voomWithQualityWeights",
                        list(counts     = dgelist,
                             design     = designMatrix,
                             plot       = mvPlot,
                             col        = 'blue',
                             var.design = var.design))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during limma voomWithQualityWeights()")
                return(NULL)
            })

            fit <- tryCatch({
                do.call("lmFit",
                        list(object = VoomElist,
                             design = designMatrix))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during linear modeling")
                return(NULL)
            })
        } else if (dupcor && !qualityWeights && !blockQW) {
            # dupcor_base analysis
            VoomElist <- tryCatch({
                do.call("voom",
                        list(counts = dgelist,
                             design = designMatrix))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during limma voom()")
                return(NULL)
            })

            corfit <- tryCatch({
                do.call("duplicateCorrelation",
                        list(object = VoomElist,
                             design = designMatrix,
                             block  = dupCorBlock))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during intra-block correlation estimation")
                return(NULL)
            })

            if (runDupCorTwice) {
                VoomElist <- tryCatch({
                    do.call("voom",
                            list(counts = dgelist,
                                 design = designMatrix,
                                 plot   = mvPlot))
                },
                error = function(e) {
                    message("Unexpected error: ", e$message, " happened during limma voomWithQualityWeights()")
                    return(NULL)
                })

                corfit <- tryCatch({
                    do.call("duplicateCorrelation",
                            list(object = VoomElist,
                                 design = designMatrix,
                                 block  = dupCorBlock))
                },
                error = function(e) {
                    message("Unexpected error: ", e$message, " happened during intra-block correlation estimation")
                    return(NULL)
                })
            }

            fit <- tryCatch({
                do.call("lmFit",
                        list(object = VoomElist,
                             design = designMatrix,
                             block  = dupCorBlock,
                             correlation = corfit$consensus.correlation))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during linear modeling")
                return(NULL)
            })
        } else if (dupcor && qualityWeights && !blockQW) {
            # dupcor_indQW analysis
            VoomElist <- tryCatch({
                do.call("voomWithQualityWeights",
                        list(counts     = dgelist,
                             design     = designMatrix))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during limma voomWithQualityWeights()")
                return(NULL)
            })

            corfit <- tryCatch({
                do.call("duplicateCorrelation",
                        list(object = VoomElist,
                             design = designMatrix,
                             block  = dupCorBlock))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during intra-block correlation estimation")
                return(NULL)
            })

            if (runDupCorTwice) {
                VoomElist <- tryCatch({
                    do.call("voomWithQualityWeights",
                            list(counts     = dgelist,
                                 design     = designMatrix,
                                 plot       = mvPlot,
                                 col        = "blue",
                                 correlation = corfit$consensus.correlation))
                },
                error = function(e) {
                    message("Unexpected error: ", e$message, " happened during limma voomWithQualityWeights()")
                    return(NULL)
                })

                corfit <- tryCatch({
                    do.call("duplicateCorrelation",
                            list(object = VoomElist,
                                 design = designMatrix,
                                 block  = dupCorBlock))
                },
                error = function(e) {
                    message("Unexpected error: ", e$message, " happened during intra-block correlation estimation")
                    return(NULL)
                })
            }

            fit <- tryCatch({
                do.call("lmFit",
                        list(object = VoomElist,
                             design = designMatrix,
                             block  = dupCorBlock,
                             correlation = corfit$consensus.correlation))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during linear modeling")
                return(NULL)
            })
        } else if (dupcor && qualityWeights && blockQW) {
            # dupcor_vdQW analysis

            VoomElist <- tryCatch({
                do.call("voomWithQualityWeights",
                        list(counts     = dgelist,
                             design     = designMatrix,
                             var.design = var.design))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during limma voomWithQualityWeights()")
                return(NULL)
            })

            corfit <- tryCatch({
                do.call("duplicateCorrelation",
                        list(object = VoomElist,
                             design = designMatrix,
                             block  = dupCorBlock))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during intra-block correlation estimation")
                return(NULL)
            })

            if (runDupCorTwice) {
                VoomElist <- tryCatch({
                    do.call("voomWithQualityWeights",
                            list(counts      = dgelist,
                                 design      = designMatrix,
                                 plot        = mvPlot,
                                 col         = "blue",
                                 correlation = corfit$consensus.correlation,
                                 var.design  = var.design))
                },
                error = function(e) {
                    message("Unexpected error: ", e$message, " happened during limma voomWithQualityWeights()")
                    return(NULL)
                })

                corfit <- tryCatch({
                    do.call("duplicateCorrelation",
                            list(object = VoomElist,
                                 design = designMatrix,
                                 block  = dupCorBlock))
                },
                error = function(e) {
                    message("Unexpected error: ", e$message, " happened during intra-block correlation estimation")
                    return(NULL)
                })
            }

            fit <- tryCatch({
                do.call("lmFit",
                        list(object = VoomElist,
                             design = designMatrix,
                             block  = dupCorBlock,
                             correlation = corfit$consensus.correlation))
            },
            error = function(e) {
                message("Unexpected error: ", e$message, " happened during linear modeling")
                return(NULL)
            })
        }

        # Run eBayes
        if (runEBayes) {
            fit <- do.call("eBayes",
                   list(fit        = fit,
                        robust     = robust,
                        proportion = proportion))
            itemAttr <- list(eBayes = TRUE)
        } else {
            itemAttr <- list(eBayes = FALSE)
        }

        if (exists("corfit")) { # Duplicate correlation was used; capture the correlation value
            cat(stringr::str_c("Duplicate Correlation = ", round(corfit$consensus.correlation, 4), "   \n"))
            attr(VoomElist, "DupCor") <- corfit$consensus.correlation
            if (!is.null(fit)) {
                attr(fit, "DupCor") <- corfit$consensus.correlation
            }
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
