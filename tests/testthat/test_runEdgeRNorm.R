context("DGEobj.utils - tests for runEdgeRNorm.R functions")
skip_if(setup_failed)
skip_if_not_installed("ggplot2")


test_that('runEdgeRNorm: runEdgeRNorm()', {
    require(ggplot2)

    # data setup
    dgeobj <- t_obj1
    dgeobj <- rmItem(dgeobj, "DGEList")
    plot_labels <- function(n = 50) {
        a <- do.call("paste0", replicate(5, sample(LETTERS, n, TRUE), FALSE))
        paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
    }

    # no plots
    ## default
    runEdgeRNorm_one_test <- runEdgeRNorm(dgeobj)
    runEdgeRNorm_one_test_DGEList <- getType(runEdgeRNorm_one_test, "DGEList")
    expect_s3_class(runEdgeRNorm_one_test, "DGEobj")
    expect_true(is.list(runEdgeRNorm_one_test_DGEList))
    expect_equal(length(runEdgeRNorm_one_test$DGEList), 2)
    expect_equal(names(runEdgeRNorm_one_test$DGEList), c("counts", "samples"))
    ## explicit FALSE
    runEdgeRNorm_one_test <- runEdgeRNorm(dgeobj, includePlot = FALSE)
    expect_s3_class(runEdgeRNorm_one_test, "DGEobj")
    ## NULL
    runEdgeRNorm_one_test <- runEdgeRNorm(dgeobj, includePlot = NULL)
    expect_s3_class(runEdgeRNorm_one_test, "DGEobj")
    ## New more list
    runEdgeRNorm_one_test <- runEdgeRNorm(runEdgeRNorm_one_test, itemName = "newList")
    expect_s3_class(runEdgeRNorm_one_test, "DGEobj")
    runEdgeRNorm_one_test_newList <- getType(runEdgeRNorm_one_test, "DGEList")
    expect_true(is.list(runEdgeRNorm_one_test_newList))
    expect_equal(names(runEdgeRNorm_one_test_newList), c("DGEList", "newList"))
    expect_equal(length(runEdgeRNorm_one_test$newList), 2)
    expect_equal(names(runEdgeRNorm_one_test$newList), c("counts", "samples"))
    # canvasXpress
    ## with samples
    runEdgeRNorm_two_test <- runEdgeRNorm(dgeobj, normMethod = "RLE", includePlot = "canvasXpress",
                                          plotLabels = plot_labels(ncol(dgeobj)))
    expect_true(is.list(runEdgeRNorm_two_test))
    expect_equal(length(runEdgeRNorm_two_test), 2)
    expect_equal(length(runEdgeRNorm_two_test[[1]]$DGEList), 2)
    expect_s3_class(runEdgeRNorm_two_test[[1]], "DGEobj")
    expect_equal(names(runEdgeRNorm_two_test[[1]]$DGEList), c("counts", "samples"))
    expect_s3_class(runEdgeRNorm_two_test[[2]], c("canvasXpress", "htmlwidget"))
    ## with wrong number of samples
    expect_warning(runEdgeRNorm_two_test <- runEdgeRNorm(dgeobj,
                                                         normMethod = "upperquartile",
                                                         includePlot = "canvasXpress",
                                                         plotLabels = plot_labels(ncol(dgeobj) - 1)),
                   regexp = "plotLabels must be a character vector with length equal to the number of columns in dgeObj.  Assigning default values.")
    expect_true(is.list(runEdgeRNorm_two_test))
    expect_equal(length(runEdgeRNorm_two_test), 2)
    expect_equal(length(runEdgeRNorm_two_test[[1]]$DGEList), 2)
    expect_s3_class(runEdgeRNorm_two_test[[1]], "DGEobj")
    expect_equal(names(runEdgeRNorm_two_test[[1]]$DGEList), c("counts", "samples"))
    expect_s3_class(runEdgeRNorm_two_test[[2]], c("canvasXpress", "htmlwidget"))
    ## with no samples
    runEdgeRNorm_two_test <- runEdgeRNorm(dgeobj, normMethod = "RLE", includePlot = TRUE)
    runEdgeRNorm_two_test_DGEList <- getType(runEdgeRNorm_two_test[[1]], "DGEList")
    expect_true(is.list(runEdgeRNorm_two_test))
    expect_equal(length(runEdgeRNorm_two_test), 2)
    expect_true(is.list(runEdgeRNorm_two_test_DGEList))
    expect_equal(length(runEdgeRNorm_two_test[[1]]$DGEList), 2)
    expect_s3_class(runEdgeRNorm_two_test[[1]], "DGEobj")
    expect_equal(names(runEdgeRNorm_two_test[[1]]$DGEList), c("counts", "samples"))
    expect_s3_class(runEdgeRNorm_two_test[[2]], c("canvasXpress", "htmlwidget"))

    # ggplot
    ## with samples
    runEdgeRNorm_two_test <- runEdgeRNorm(dgeobj, normMethod = "RLE", includePlot = "ggplot",
                                          plotLabels = plot_labels(ncol(dgeobj)))
    expect_true(is.list(runEdgeRNorm_two_test))
    expect_equal(length(runEdgeRNorm_two_test), 2)
    expect_equal(length(runEdgeRNorm_two_test[[1]]$DGEList), 2)
    expect_s3_class(runEdgeRNorm_two_test[[1]], "DGEobj")
    expect_equal(names(runEdgeRNorm_two_test[[1]]$DGEList), c("counts", "samples"))
    expect_s3_class(runEdgeRNorm_two_test[[2]], c("gg", "ggplot"))

    ## with no samples
    runEdgeRNorm_two_test <- runEdgeRNorm(dgeobj, normMethod = "RLE", includePlot = "ggplot")
    runEdgeRNorm_two_test_DGEList <- getType(runEdgeRNorm_two_test[[1]], "DGEList")
    expect_true(is.list(runEdgeRNorm_two_test))
    expect_equal(length(runEdgeRNorm_two_test), 2)
    expect_true(is.list(runEdgeRNorm_two_test_DGEList))
    expect_equal(length(runEdgeRNorm_two_test[[1]]$DGEList), 2)
    expect_s3_class(runEdgeRNorm_two_test[[1]], "DGEobj")
    expect_equal(names(runEdgeRNorm_two_test[[1]]$DGEList), c("counts", "samples"))
    expect_s3_class(runEdgeRNorm_two_test[[2]], c("gg", "ggplot"))

    # Testing asserts
    ## DGEobj
    msg <- "dgeObj must be of class 'DGEobj'."
    expect_error(runEdgeRNorm(),
                 regexp = msg)
    expect_error(runEdgeRNorm(NULL),
                 regexp = msg)
    expect_error(runEdgeRNorm(list()),
                 regexp = msg)
    ## normMethod
    msg <- "normMethod must be only one of the following values 'TMM', 'RLE', 'upperquartile', 'none'."
    expect_error(runEdgeRNorm(dgeobj, normMethod = NULL),
                 regexp = msg)
    expect_error(runEdgeRNorm(dgeobj, normMethod = "abc"),
                 regexp = msg)
    expect_error(runEdgeRNorm(dgeobj, normMethod = 123),
                 regexp = msg)
    expect_error(runEdgeRNorm(dgeobj, normMethod = c("TMM", "RLE")),
                 regexp = msg)
    expect_error(runEdgeRNorm(runEdgeRNorm_test),
                 regexp = "object 'runEdgeRNorm_test' not found")
    ## itemName
    msg <-  "itemName must be a singular, unique and not NULL character value."
    ### null name
    expect_error(runEdgeRNorm(dgeobj, itemName = NULL),
                 regexp = msg)
    ### duplicate defult name
    expect_error(runEdgeRNorm(t_obj1),
                 regexp = msg)
    ### duplicate new name
    dgeobj <- runEdgeRNorm(dgeobj, itemName = "newList")
    expect_error(runEdgeRNorm(dgeobj, itemName = "newList"),
                 regexp = msg)
    ## includePlot
    msg <- "includePlot must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE."
    expect_warning(runEdgeRNorm(dgeobj, normMethod = "RLE", includePlot = c(TRUE, FALSE)),
                   regexp = msg)
    expect_warning(runEdgeRNorm(dgeobj, normMethod = "RLE", includePlot = "abc"),
                   regexp = msg)
    expect_warning(runEdgeRNorm(dgeobj, normMethod = "RLE", includePlot = c("canvasXpress", "ggplot")),
                   regexp = msg)
})
