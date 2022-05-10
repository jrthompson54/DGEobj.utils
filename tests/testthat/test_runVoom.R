context("DGEobj.utils - tests for runVoom.R functions")

skip_on_cran()
skip_on_ci()
skip_if_not_installed("statmod")

test_that('runVoom.R: runVoom()', {
    dgeObj <- t_obj1
    design <- getItem(dgeObj, "design")
    designMatrix <- model.matrix(~ 0 + ReplicateGroup, design)
    dgeObj <- addItem(dgeObj   = dgeObj,
                      item     = designMatrix,
                      itemName = "designMat",
                      itemType = "designMatrix")

    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           qualityWeights   = FALSE,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing indQW analysis
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           qualityWeights   = TRUE,
                           mvPlot           = FALSE,
                           runEBayes        = FALSE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing blockedQW analysis
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           qualityWeights   = TRUE,
                           var.design       = designMatrix,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing dupcor_indQW analysis
    dupcorBlock <- rep(1:6, 8)
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           dupCorBlock      = dupcorBlock,
                           qualityWeights   = TRUE,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing dupcor_base analysis
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           runDupCorTwice   = TRUE,
                           qualityWeights   = FALSE,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing dupcor_vdQW analysis
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           dupCorBlock      = dupcorBlock,
                           qualityWeights   = TRUE,
                           var.design       = designMatrix,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing assert statements
    ## dgeObj
    msg <- "dgeObj must be specified and must be of class 'DGEobj'."
    expect_error(runVoom(),
                 regexp = msg)
    expect_error(runVoom(dgeObj = NULL),
                 regexp = msg)
    ## designMatrixName
    msg <- "designMatrixName must be specified and must be one of the items in dgeObj. Use names(dgeObj) to check for available options."
    expect_error(runVoom(dgeObj = dgeObj, designMatrixName = "xyz"),
                 regexp = msg,
                 fixed = TRUE)
    expect_error(runVoom(dgeObj = dgeObj),
                 regexp = msg,
                 fixed = TRUE)
    expect_error(runVoom(dgeObj = dgeObj, designMatrixName = NULL),
                 regexp = msg,
                 fixed = TRUE)
    expect_error(runVoom(dgeObj = dgeObj, designMatrixName = 123),
                 regexp = msg,
                 fixed = TRUE)
    expect_error(runVoom(dgeObj = dgeObj, designMatrixName = c("designMat", "designMat")),
                 regexp = msg,
                 fixed = TRUE)
    ## runEBayes
    msg <- "runEBayes must be a singular logical value. Assigning default value TRUE"
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           runEBayes        = NULL),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           runEBayes        = "FALSE"),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           runEBayes        = c(FALSE, FALSE)),
                   regexp = msg)
    ## runDupCorTwice
    msg <- "runDupCorTwice must be a singular logical value. Assigning default value TRUE"
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           runDupCorTwice   = NULL),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           runDupCorTwice        = "FALSE"),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           runDupCorTwice   = c(FALSE, FALSE)),
                   regexp = msg)
    ## qualityWeights
    msg <- "qualityWeights must be a singular logical value. Assigning default value TRUE"
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           qualityWeights   = NULL),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           qualityWeights        = "FALSE"),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           qualityWeights   = c(FALSE, FALSE)),
                   regexp = msg)
    ## mvPlot
    msg <- "mvPlot must be a singular logical value. Assigning default value TRUE"
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           mvPlot           = NULL),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           mvPlot           = "FALSE"),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           mvPlot           = c(FALSE, FALSE)),
                   regexp = msg)
    ## robust
    msg <- "robust must be a singular logical value. Assigning default value TRUE"
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           robust           = NULL),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           robust           = "FALSE"),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           robust           = c(FALSE, FALSE)),
                   regexp = msg)
    ## runEBayes
    msg <- "proportion must be a singular numeric value. Assigning default value 0.01"
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           proportion        = NULL),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           proportion       = "0.01"),
                   regexp = msg)
    expect_warning(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           proportion       = c(0.01, 0.01)),
                   regexp = msg)
    ## dubCorBlock & var.design
    msg <- "runVoom did not complete successfully due to: "
    expect_message(runVoom(dgeObj             = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = "abc"),
                   regexp = msg)
    expect_message(runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           var.design       = "abc"),
                   regexp = msg)
})
