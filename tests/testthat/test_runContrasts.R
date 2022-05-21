context("DGEobj.utils - tests for runContrasts.R functions")
skip_if(setup_failed)


test_that('runContrasts.R: runContrasts()', {
    # Define the named contrasts from design matrix column names
    contrastList  <- list(Sham_vs_BDL     = "ReplicateGroupSham - ReplicateGroupBDL",
                          Sham_vs_EXT1024 = "ReplicateGroupSham  - ReplicateGroupBDL_EXT.1024",
                          Sham_vs_Nint    = "ReplicateGroupSham - ReplicateGroupBDL_Nint",
                          Sham_vs_Sora    = "ReplicateGroupSham - ReplicateGroupBDL_Sora")

    if (requireNamespace("statmod", quietly = TRUE)) {
        dgeObj_output <- runContrasts(dgeObj              = t_obj1,
                                      designMatrixName    = "ReplicateGroupDesign",
                                      contrastList        = contrastList,
                                      contrastSetName     = "ReplicateGroup_Contrasts")
        expect_s3_class(dgeObj_output, "DGEobj")

        contrastList  <- list(EXT1024_vs_Sham    = "ReplicateGroupBDL_EXT.1024 - ReplicateGroupSham",
                              BDL_vs_EXT1024     = "ReplicateGroupBDL  - ReplicateGroupBDL_EXT.1024",
                              EXT1024_vs_Nint    = "ReplicateGroupBDL_EXT.1024 - ReplicateGroupBDL_Nint",
                              EXT1024_vs_Sora    = "ReplicateGroupBDL_EXT.1024 - ReplicateGroupBDL_Sora")
        if (requireNamespace("IHW", quietly = TRUE)) {
            dgeObj_output <- runContrasts(dgeObj              = t_obj1,
                                          designMatrixName    = "ReplicateGroupDesign",
                                          contrastList        = contrastList,
                                          contrastSetName     = "ReplicateGroup_Contrasts",
                                          runTopTreat         = TRUE,
                                          qValue              = TRUE,
                                          IHW                 = TRUE,
                                          verbose             = TRUE)
            expect_s3_class(dgeObj_output, "DGEobj")
        }
    } else {
        expect_error(runContrasts(dgeObj              = t_obj1,
                                  designMatrixName    = "ReplicateGroupDesign",
                                  contrastList        = contrastList,
                                  contrastSetName     = "ReplicateGroup_Contrasts"),
                     "'statmod' package is required to run eBayes calculations")
    }

    # robust is disabled
    dgeObj_output <- runContrasts(dgeObj              = t_obj1,
                                  designMatrixName    = "ReplicateGroupDesign",
                                  contrastList        = contrastList,
                                  contrastSetName     = "ReplicateGroup_Contrasts",
                                  robust              = FALSE)
    expect_s3_class(dgeObj_output, "DGEobj")

    # testing assert statements
    ## dgeobj
    expect_error(runContrasts(),
                 regexp = "dgeObj must be specified and should be of class 'DGEobj'.")
    expect_error(runContrasts(dgeObj = NULL),
                 regexp = "dgeObj must be specified and should be of class 'DGEobj'.")
    ## designMatrixName
    msg <- "designMatrixName must be a signular character value and one of dgeobj names."
    expect_error(runContrasts(dgeObj = t_obj1),
                 regexp = msg)
    expect_error(runContrasts(dgeObj = t_obj1, designMatrixName = NULL),
                 regexp = msg)
    expect_error(runContrasts(dgeObj = t_obj1, designMatrixName = 123),
                 regexp = msg)
    expect_error(runContrasts(dgeObj = t_obj1, designMatrixName = "abc"),
                 regexp = msg)
    expect_error(runContrasts(dgeObj = t_obj1, designMatrixName = c("ReplicateGroup", "ReplicateGroup")),
                 regexp = msg)
    expect_error(runContrasts(dgeObj           = t_obj1,
                              designMatrixName = "ReplicateGroupDesign",
                              contrastList     = "XYZ"),
                 regexp = "contrastList must specified and must be a named list.")
    expect_error(runContrasts(dgeObj              = t_obj1,
                              designMatrixName    = "ReplicateGroupDesign",
                              contrastList        = contrastList,
                              foldChangeThreshold = -1),
                 regexp = "foldChangeThreshold must be greater than or equal to 0.")
    expect_error(runContrasts(dgeObj              = t_obj1,
                              designMatrixName    = "ReplicateGroupDesign",
                              contrastList        = contrastList,
                              runTopTable         = FALSE,
                              runTopTreat         = FALSE),
                 regexp = "One of runTopTable or runTopTreat must be TRUE.")
    ## runEBayes
    msg <- "runEBayes must be a singular logical value. Assigning default value TRUE"
    if (requireNamespace("statmod", quietly = TRUE)) {
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    runEBayes        = NULL),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    runEBayes        = "FALSE"),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    runEBayes        = c(FALSE, FALSE)),
                       regexp = msg)
        ## robust
        msg <- "robust must be a singular logical value. Assigning default value TRUE"
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    robust           = NULL),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    robust           = "FALSE"),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    robust           = c(FALSE, FALSE)),
                       regexp = msg)
        ## contrastSetName
        msg <- "contrastSetName must be a character value. Assigning default value: ReplicateGroupDesign_fit"

        expect_warning(
            expect_error(runContrasts(dgeObj           = t_obj1,
                                      designMatrixName = "ReplicateGroupDesign",
                                      contrastList     = contrastList),
                         regexp = "The contrastSetName already exists in dgeObj."),
            msg)

        test <- DGEobj::rmItem(t_obj1, 'ReplicateGroupDesign_fit_cm')
        test <- DGEobj::rmItem(test,   'ReplicateGroupDesign_fit_cf')

        expect_warning(runContrasts(dgeObj           = test,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = test,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  =  c("ReplicateGroup_Contrasts", "ReplicateGroup_Contrasts")),
                       regexp = msg)
        ## proportion
        msg <- "proportion must be a singular numeric value. Assigning default value 0.01"
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    proportion       = NULL),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    proportion       = "0.1"),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    proportion       = c(0.1, 0.01)),
                       regexp = msg)
        ## qValue
        msg <- "qValue must be a singular logical value. Assigning default value FALSE"
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    qValue           = NULL),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    qValue           = "FALSE"),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    qValue           = c(FALSE, FALSE)),
                       regexp = msg)
        ## IHW
        msg <- "IHW must be a singular logical value. Assigning default value FALSE"
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    IHW              = NULL),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    IHW              = "FALSE"),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    IHW              = c(FALSE, FALSE)),
                       regexp = msg)
        ## verbose
        msg <- "verbose must be a singular logical value. Assigning default value FALSE"
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    verbose          = NULL),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    verbose          = "FALSE"),
                       regexp = msg)
        expect_warning(runContrasts(dgeObj           = t_obj1,
                                    designMatrixName = "ReplicateGroupDesign",
                                    contrastList     = contrastList,
                                    contrastSetName  = "ReplicateGroup_Contrasts",
                                    verbose          = c(FALSE, FALSE)),
                       regexp = msg)
    }
})
