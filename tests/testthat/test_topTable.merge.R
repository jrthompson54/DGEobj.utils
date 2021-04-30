context("DGEobj.utils - tests for topTable.merge.R functions")


test_that("topTable.merge.R: topTable.merge()", {
    # creating toptables list
    contrastList   <- getType(t_obj1, "topTable")
    contrast_table <- topTable.merge(contrastList = contrastList, digits = 2)
    expect_setequal(object = colnames(contrast_table),
                    expected = apply(X        = expand.grid(c("logFC", "AveExpr", "P.Value", "adj.P.Val"), names(contrastList)),
                                     MARGIN   =  1,
                                     FUN      =  paste,
                                     collapse = "_"))

    # testing assert statements
    ## contrastList
    msg <- "contrastList must be specified, be of class 'list' and be a named list specifically, and include items of class 'data.frame'."
    dataframe1 <- getType(t_obj1, "topTable")[[3]][1:10,]
    contrastList2 <- contrastList
    contrastList2$dataframe1 <- dataframe1
    dataframe2 <- getType(t_obj1, "topTable")[[3]]
    contrastList3 <- contrastList
    colnames(dataframe2) <- 1:length(dataframe2)
    contrastList3$dataframe2 <- dataframe2
    expect_error(topTable.merge(),
                 regexp = msg)
    expect_error(topTable.merge(contrastList = NULL),
                 regexp = msg)
    expect_error(topTable.merge(contrastList = contrastList2),
                 regexp = msg)
    expect_error(topTable.merge(contrastList = contrastList3),
                 regexp = msg)
    expect_error(topTable.merge(contrastList = contrastList, digits = 1:5),
                 regexp = "digits must be either of length 1 or the same length as colNames.")
})
