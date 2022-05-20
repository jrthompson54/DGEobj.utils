require(testthat)
require(stats)
require(DGEobj)
require(DGEobj.utils)

setup_failed <- TRUE

# NOTE: if the Bioconductor packages are missing from the platform we will NOT perform tests
if (requireNamespace('zFPKM',       quietly = TRUE) &&
    requireNamespace('sva',         quietly = TRUE) &&
    requireNamespace('RNASeqPower', quietly = TRUE) &&
    requireNamespace('qvalue',      quietly = TRUE) &&
    requireNamespace('limma',       quietly = TRUE) &&
    requireNamespace('IHW',         quietly = TRUE) &&
    requireNamespace('edgeR',       quietly = TRUE) &&
    requireNamespace('biomaRt',     quietly = TRUE)) {

    require(statmod)
    require(RNASeqPower)
    require(edgeR)

    t_obj1 <- readRDS(system.file("exampleObj.RDS", package = "DGEobj", mustWork = TRUE))

    setup_failed <- FALSE
}

if (setup_failed) {
    message('Test Setup Failed - likely due to missing suggested packages.  Tests will be skipped.')
}
