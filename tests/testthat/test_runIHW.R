context("DGEobj.utils - tests for runIHW.R functions")
skip_on_cran()


test_that('runIHW: runIHW()', {
    runIHW_contrastList <- getType(t_obj1, "topTable")[1:2]

    runIHW_test_one <- runIHW(runIHW_contrastList)

    expect_true(is.list(runIHW_test_one))
    expect_equal(length(runIHW_test_one), 2)
    expect_equal(names(runIHW_test_one), c("contrasts", "ihwObj"))
    expect_s4_class(runIHW_test_one$ihwObj[[1]], "ihwResult")
    expect_s4_class(runIHW_test_one$ihwObj[[2]], "ihwResult")

    expect_error(runIHW(runIWH_test_two),
                 regexp = "object 'runIWH_test_two' not found")
    ## alpha
    msg <- "alpha must be a singular numeric value between 0 and 1. Assigning default value 0.1"
    expect_warning(runIHW(runIHW_contrastList,
                          alpha = NULL),
                   regexp = msg)
    expect_warning(runIHW(runIHW_contrastList,
                          alpha  = "0.1"),
                   regexp = msg)
    expect_warning(runIHW(runIHW_contrastList,
                          alpha = c(0.1, 0.01)),
                   regexp = msg)
    expect_warning(runIHW(runIHW_contrastList,
                          alpha = -1),
                   regexp = msg)
    expect_warning(runIHW(runIHW_contrastList,
                          alpha = 2),
                   regexp = msg)
    ## FDRthreshold
    msg <- "FDRthreshold must be a singular numeric value between 0 and 1. Assigning default value 0.1"
    expect_warning(runIHW(runIHW_contrastList,
                          FDRthreshold = NULL),
                   regexp = msg)
    expect_warning(runIHW(runIHW_contrastList,
                          FDRthreshold  = "0.1"),
                   regexp = msg)
    expect_warning(runIHW(runIHW_contrastList,
                          FDRthreshold = c(0.1, 0.01)),
                   regexp = msg)
    expect_warning(runIHW(runIHW_contrastList,
                          FDRthreshold = 2),
                   regexp = msg)
    expect_warning(runIHW(runIHW_contrastList,
                          FDRthreshold = -1),
                   regexp = msg)
})


test_that('runIHW: incorrect usage', {
    expect_error(runIHW(),
                 regexp = "contrastList must be specified and should be of class 'List'")
    expect_error(runIHW(NULL),
                 regexp = "contrastList must be specified and should be of class 'List'")
})
