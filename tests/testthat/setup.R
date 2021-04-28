require(testthat)
require(stats)

require(DGEobj)
require(DGEobj.utils)
require(ggplot2)


t_obj1 <- readRDS(system.file("exampleObj.RDS", package = "DGEobj", mustWork = TRUE))
