context("DGEobj.utils - tests for lowIntFilter.R functions")
skip_if_not_installed('GenomicRanges')

test_that('lowIntFilter: lowIntFilter()', {
    lowIntFilter_one_test <- lowIntFilter(t_obj1, countThreshold = 10)
    expect_s3_class(lowIntFilter_one_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_one_test$counts), 951)

    #verbose is enabled
    if (requireNamespace("zFPKM", quietly = TRUE)) {
        lowIntFilter_one_test <- lowIntFilter(t_obj1, verbose = TRUE, zfpkmThreshold = -3.0)
        expect_s3_class(lowIntFilter_one_test, "DGEobj")
        expect_equal(nrow(lowIntFilter_one_test$counts), 933)
    }

    # NULL gene length
    lowIntFilter_one_test <- lowIntFilter(t_obj1, tpmThreshold = 1, verbose = TRUE)
    expect_error(lowIntFilter(t_obj1,
                              tpmThreshold = 1,
                              verbose = TRUE,
                              geneLength = t_obj1$geneData$ExonLength[1:100]),
                 regexp = "geneLength must be specified and should be the same length as the number of rows in counts.")
    expect_s3_class(lowIntFilter_one_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_one_test$counts), 951)
    # test TPM threshold
    lowIntFilter_one_test <- lowIntFilter(t_obj1,
                                          countThreshold = 10,
                                          verbose = TRUE)
    expect_s3_class(lowIntFilter_one_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_one_test$counts), 951)

    if (requireNamespace("zFPKM", quietly = TRUE)) {
        lowIntFilter_two_test <- lowIntFilter(t_obj1, zfpkmThreshold = -3.0)
        expect_s3_class(lowIntFilter_two_test, "DGEobj")
        expect_equal(nrow(lowIntFilter_two_test$counts), 933)
    }

    lowIntFilter_three_test <- lowIntFilter(t_obj1,
                                            fpkThreshold = 5,
                                            verbose = TRUE)

    expect_s3_class(lowIntFilter_three_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_three_test$counts), 933)

    if (requireNamespace("zFPKM", quietly = TRUE)) {
        lowIntFilter_four_test <- lowIntFilter(t_obj1, countThreshold = 10, zfpkmThreshold = -3.0, sampleFraction = 0.75)
        expect_s3_class(lowIntFilter_four_test, "DGEobj")
        expect_equal(nrow(lowIntFilter_four_test$counts), 881)
    }

    expect_error(lowIntFilter(lowIntFilter_five_test),
                 regexp = "object 'lowIntFilter_five_test' not found")
    expect_error(lowIntFilter(NULL),
                 regexp = "dgeObj must be of class 'DGEobj'.")
    expect_error(lowIntFilter(),
                 regexp = "dgeObj must be of class 'DGEobj'.")
    # Testing assert
    ## sampleFraction
    msg <- "sampleFraction must be a singular numeric value. Assigning default value 0.5"
    expect_warning(lowIntFilter(t_obj1,
                                fpkThreshold   = 5,
                                sampleFraction = NULL),
                   regexp = msg)
    expect_warning(lowIntFilter(t_obj1,
                                fpkThreshold    = 5,
                                sampleFraction  = "0.5"),
                   regexp = msg)
    expect_warning(lowIntFilter(t_obj1,
                                fpkThreshold   = 5,
                                sampleFraction = c(0.5, 0.5)),
                   regexp = msg)
    ## verbose
    msg <- "verbose must be a singular logical value. Assigning default value FALSE"
    expect_warning(lowIntFilter(t_obj1,
                                fpkThreshold = 5,
                                verbose      = NULL),
                   regexp = msg)
    expect_warning(lowIntFilter(t_obj1,
                                fpkThreshold = 5,
                                verbose      = "FALSE"),
                   regexp = msg)
    expect_warning(lowIntFilter(t_obj1,
                                fpkThreshold = 5,
                                verbose         = c(FALSE, FALSE)),
                   regexp = msg)
})


test_that('lowIntFilter: incorrect usage', {
    testthat::skip_if_not_installed("zFPKM")
    expect_error(lowIntFilter(t_obj1, zfpkmThreshold = 3.0, tpmThreshold = 1),
                 regexp = "Must use zfpkmThreshold or tpmThreshold, but not both.")
})
