context("DGEobj.utils - tests for extractCol.R functions")
skip_if(setup_failed)


test_that('extractCol: extractCol()', {
    extractCol_contrastList <- getType(t_obj1, "topTable")[1:2]
    extractCol_one_test <- extractCol(extractCol_contrastList, colName = "P.Value")

    expect_true(is.data.frame(extractCol_one_test))
    expect_equal(nrow(extractCol_one_test), 951)
    expect_equal(ncol(extractCol_one_test), 2)
    expect_equal(names(extractCol_one_test), c("BDL_vs_Sham", "EXT1024_vs_BDL"))

    extractCol_two_test <- extractCol(extractCol_contrastList, colName = "P.Value", robust = FALSE)

    expect_true(is.data.frame(extractCol_two_test))
    expect_equal(nrow(extractCol_two_test), 951)
    expect_equal(ncol(extractCol_two_test), 2)
    expect_equal(names(extractCol_two_test), c("BDL_vs_Sham", "EXT1024_vs_BDL"))

    expect_error(extractCol(extractCol_three_test),
                 regexp = "object 'extractCol_three_test' not found")
    # Testing assert
    ## contrastList
    msg <- "contrastList must be a list of data.frames which all have the same colnames and same row counts."
    dataframe1 <- getType(t_obj1, "topTable")[[3]][1:10,]
    extractCol_contrastList2 <- extractCol_contrastList
    extractCol_contrastList2$dataframe1 <- dataframe1
    dataframe2 <- getType(t_obj1, "topTable")[[3]]
    extractCol_contrastList3 <- extractCol_contrastList
    colnames(dataframe2) <- 1:length(dataframe2)
    extractCol_contrastList3$dataframe2 <- dataframe2
    expect_error(extractCol(),
                   regexp = msg)
    expect_error(extractCol(contrastList = NULL),
                   regexp = msg)
    expect_error(extractCol(contrastList = extractCol_contrastList2),
                   regexp = msg)
    expect_error(extractCol(contrastList = extractCol_contrastList3),
                   regexp = msg)
    ## robust
    msg <- "robust must be a singular logical value. Assigning default value TRUE."
    expect_warning(extractCol(extractCol_contrastList, colName = "P.Value",
                              robust = NULL),
                   regexp = msg)
    expect_warning(extractCol(extractCol_contrastList, colName = "P.Value",
                              robust = "FALSE"),
                   regexp = msg)
    expect_warning(extractCol(extractCol_contrastList, colName = "P.Value",
                              robust = c(FALSE, FALSE)),
                   regexp = msg)
})
