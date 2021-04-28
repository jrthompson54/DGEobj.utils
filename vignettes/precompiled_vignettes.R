# Precompiled vignettes

library(knitr)

knit("vignettes/DGEobj.utils_Workflow.Rmd.orig", "vignettes/DGEobj.utils_Workflow.Rmd")

# remove file path such that vignettes will build with figures
replace <- readLines("vignettes/DGEobj.utils_Workflow.Rmd")
replace <- gsub("](vignettes/figure-", "](figure-", replace, fixed = T)
fileConn <- file("vignettes/DGEobj.utils_Workflow.Rmd")
writeLines(replace, fileConn)
close(fileConn)
