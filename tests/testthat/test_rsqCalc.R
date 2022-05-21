context("DGEobj.utils - tests for rsqCalc.R functions")
skip_if(setup_failed)


test_that("rsqCalc.R: rsqCalc()", {
    dgelist <- getItem(t_obj1, "DGEList")
    log2cpm <- cpm(dgelist, log = TRUE)
    rsq <- rsqCalc(log2cpm, t_obj1$ReplicateGroupDesign_fit)
    expect_type(rsq, "double")

    ## normMatrix
    msg <- "normMatrix must be of class 'data.frame' or 'matrix'."
    expect_error(rsqCalc(), regexp = msg)
    expect_error(rsqCalc(NULL), regexp = msg)
    expect_error(rsqCalc(1:10), regexp = msg)
    ## fit
    msg <- "fit must be of class 'MArrayLM'."
    expect_error(rsqCalc(log2cpm), regexp = msg)
    expect_error(rsqCalc(log2cpm, NULL), regexp = msg)
    expect_error(rsqCalc(log2cpm, 1:10), regexp = msg)
    expect_error(rsqCalc(as.matrix(LETTERS), t_obj1$ReplicateGroupDesign_fit), regexp = "All of the entries in normMatrix must be numeric.")
})
