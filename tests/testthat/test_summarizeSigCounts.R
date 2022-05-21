context("DGEobj.utils - tests for summarizeSigCounts.R functions")
skip_if(setup_failed)


test_that('summarizeSigCounts.R: summarizeSigCounts()', {
    myTopTables <- getType(t_obj1, "topTable")
    summarizedSigCounts <- summarizeSigCounts(myTopTables)

    expect_true(is.matrix(summarizedSigCounts))
    expect_equal(nrow(summarizedSigCounts), 4)
    expect_equal(ncol(summarizedSigCounts), 5)
    expect_equal(rownames(summarizedSigCounts), c("BDL_vs_Sham", "EXT1024_vs_BDL",
                                                  "Nint_vs_BDL","Sora_vs_BDL"))
    expect_equal(colnames(summarizedSigCounts), c("P.Value", "adj.P.Val", "Qvalue",
                                                  "qvalue.lfdr", "ihw.adj_pvalue"))

    # specify fcThreshold
    summarizedSigCounts_two <- summarizeSigCounts(myTopTables,
                                                  fcThreshold = 0.1)

    expect_true(is.matrix(summarizedSigCounts_two))
    expect_equal(nrow(summarizedSigCounts_two), 4)
    expect_equal(ncol(summarizedSigCounts_two), 5)
    expect_equal(rownames(summarizedSigCounts_two), c("BDL_vs_Sham", "EXT1024_vs_BDL",
                                                      "Nint_vs_BDL","Sora_vs_BDL"))
    expect_equal(colnames(summarizedSigCounts_two), c("P.Value", "adj.P.Val", "Qvalue",
                                                      "qvalue.lfdr", "ihw.adj_pvalue"))

    expect_error(summarizeSigCounts(myTopTable),
                 regexp = "object 'myTopTable' not found")
})


test_that('summarizeSigCounts.R: incorrect usage', {
    expect_error(summarizeSigCounts(),
                 regexp = "argument \"contrastList\" is missing, with no default")
})
