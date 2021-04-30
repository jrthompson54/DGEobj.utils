context("DGEobj.utils - tests for runSVA.R functions")
skip_on_cran()


test_that("runSVA.R: runSVA()", {
    # method = leek
    dgeObj_sva <- runSVA(dgeObj = t_obj1, designMatrixName = "ReplicateGroupDesign")
    expect_s3_class(dgeObj_sva, "DGEobj")

    # method = be
    dgeObj_sva <- runSVA(dgeObj = t_obj1, designMatrixName = "ReplicateGroupDesign", method = "be")
    expect_s3_class(dgeObj_sva, "DGEobj")

    # method = be and custom na.sv
    dgeObj_sva <- runSVA(dgeObj = t_obj1, designMatrixName = "ReplicateGroupDesign", n.sv = 10, method = "be")
    expect_s3_class(dgeObj_sva, "DGEobj")

    # method = be and large custom na.sv
    expect_message(dgeObj_sva <- runSVA(dgeObj = t_obj1,
                                        designMatrixName = "ReplicateGroupDesign",
                                        n.sv = 10000, method = "be"),
                   regexp = "runSVA failed due to:")
    expect_s3_class(dgeObj_sva, "DGEobj")

    ## dge
    msg <- "dgeObj must be specified, be of class 'DGEobj', and should have a 'design' attribute."
    expect_error(runSVA(),
                 regexp = msg)
    expect_error(runSVA(NULL),
                 regexp = msg)
    ## designMatrixName
    msg <- "designMatrixName must be specified, should be of class 'character', and must exist as an attribute on the dgeObj."
    expect_error(runSVA(dgeObj = t_obj1, designMatrixName = NULL),
                 regexp = msg)
    expect_error(runSVA(dgeObj = t_obj1),
                 regexp = msg)
    expect_error(runSVA(dgeObj = t_obj1, designMatrixName = 123),
                 regexp = msg)
    expect_error(runSVA(dgeObj = t_obj1, designMatrixName = c("ReplicateGroupDesign", "ReplicateGroupDesign")),
                 regexp = msg)
    expect_error(runSVA(dgeObj = t_obj1, designMatrixName = "ReplicateGroupDesign", method =  "xyz"),
                 regexp = "method must be one of 'leek' or 'be'.")
})
