context("DGEobj.utils - tests for runPower.R functions")
skip_on_cran()
skip_if_not_installed("statmod")
skip_if_not_installed("RNASeqPower")

test_that("runPower.R: runPower()", {
    # data setup
    designMatrix <- model.matrix(~ 0 + ReplicateGroup, getItem(t_obj1, "design"))

    # no plots
    ## default value
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix)
    expect_s3_class(power_plot, "data.frame")
    ## NULL value
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix, includePlots = NULL)
    expect_s3_class(power_plot, "data.frame")
    ## FALSE value
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix, includePlots = FALSE)
    expect_s3_class(power_plot, "data.frame")
    #expect_equal(dim(power_plot), c())

    # with plots
    ## canvasXpress
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix, includePlots = "canvasXpress")
    expect_type(power_plot, "list")
    expect_s3_class(power_plot$ROC, c("canvasXpress", "htmlwidget"))
    expect_s3_class(power_plot$NvP, c("canvasXpress", "htmlwidget"))
    expect_s3_class(power_plot$PowerData, "data.frame")
    ## TRUE
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix, includePlots = TRUE)
    expect_type(power_plot, "list")
    expect_s3_class(power_plot$ROC, c("canvasXpress", "htmlwidget"))
    expect_s3_class(power_plot$NvP, c("canvasXpress", "htmlwidget"))
    expect_s3_class(power_plot$PowerData, "data.frame")
    ## ggplot
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix, includePlots = "ggplot")
    expect_type(power_plot, "list")
    expect_s3_class(power_plot$ROC, c("gg", "ggplot"))
    expect_s3_class(power_plot$NvP, c("gg", "ggplot"))
    expect_s3_class(power_plot$PowerData, "data.frame")

    # Testing asserts
    ## countsMatrix
    ### Wrong data
    expect_error(runPower(counts = t_obj1),
                 regexp = "countsMatrix must be specified and must be of class matrix or dataframe.")
    ### missing data
    expect_error(runPower(),
                 regexp = "countsMatrix must be specified and must be of class matrix or dataframe.")
    ### NULL data
    expect_error(runPower(counts = NULL),
                 regexp = "countsMatrix must be specified and must be of class matrix or dataframe.")
    ## designMatrix
    ### Wrong data
    expect_error(runPower(counts = t_obj1$counts, designMatrix = t_obj1),
                 regexp = "designMatrix must be specified and must be of class matrix or dataframe.")
    ### missing data
    expect_error(runPower(counts = t_obj1$counts),
                 regexp = "designMatrix must be specified and must be of class matrix or dataframe.")
    ### NULL data
    expect_error(runPower(counts = t_obj1$counts, designMatrix = NULL),
                 regexp = "designMatrix must be specified and must be of class matrix or dataframe.")
    ## includePlots
    msg <- "includePlots must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE."
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, includePlots = c(TRUE, FALSE)),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, includePlots = "abc"),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, includePlots = c("canvasXpress", "ggplot")),
                   regexp = msg)
    ## depth
    msg <- "depth must be a vector of 3 integer values. Assigning default values 10, 100, 1000."
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, depth = c("10", "100", "1000")),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, depth = 100),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, depth = NULL),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, depth = c(10, 100, 1000, 10000)),
                   regexp = msg)
    ## N
    msg <- "N must be a vector of 4 integer values. Assigning default values 3, 6, 10, 20."
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, N = c("3", "6", "10", "20")),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, N = 3),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, N = NULL),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, N = c(3, 6, 10, 20, 30)),
                   regexp = msg)
    ## FDR
    msg <- "FDR must be a vector of 2 integer values. Assigning default values 0.05, 0.1."
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, FDR = c("0.05", "0.1")),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, FDR = 0.05),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, FDR = NULL),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, FDR = c(0.05, 0.1, 0.55)),
                   regexp = msg)
    ## effectiveSize
    msg <- "effectiveSize must be a vector of 3 integer values. Assigning default values 1.2, 1.5, 2."
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, effectSize = c("1.2", "1.5", "2")),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, effectSize = 1.2),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, effectSize = NULL),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, effectSize = c(1.2, 1.5, 2, 4)),
                   regexp = msg)
})
