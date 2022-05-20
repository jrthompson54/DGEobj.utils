context("DGEobj.utils - tests for runQvalue.R functions")
skip_if(setup_failed)


test_that("runQvalue.R: runQvalue()", {
    contrast_list <- getType(t_obj1, "topTable")
    contrast_list_with_qvalue <- runQvalue(contrastList = contrast_list)
    expect_type(contrast_list_with_qvalue, "list")
    expect_equal(length(contrast_list_with_qvalue), length(contrast_list))

    contrasts_names <- names(contrast_list_with_qvalue)
    for (contrast in contrasts_names) {
        expect_true(all(c("Qvalue", "qvalue.lfdr") %in% names(contrast_list_with_qvalue[[contrast]])))
    }
    msg <- "contrastList must be of class 'list'."
    expect_error(runQvalue(),
                 regexp = msg)
    expect_error(runQvalue(contrastList = NULL),
                 regexp = msg)
    expect_error(runQvalue(contrastList = contrast_list, pvalField = "xyz"),
                 regexp = "pvalField must exist as an item in contrastList.")
})
