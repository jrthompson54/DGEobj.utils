## Comments from Maintainer

* Resolved CRAN check errors, including the noSuggests
* Some systems do show one NOTE however, as the suggested Bioconductor packages are not available for checking
* Made all Bioconductor packages in suggests required for testing, will not run tests if any are missing
* Also updated our process to check ALL RHub platforms, to so hopefully we do not run into unexpected issues

---  

## Test environments

RStudio Server Pro (ubuntu 18.04.2)  

* R 3.6.3
* R 4.0.5
* R 4.1.3
* R 4.2.0

Circle-CI

* R 4.0.5
* R 4.1.3
* rocker/verse:latest

WinBuilder

* devtools::check_win_devel()  
* devtools::check_win_release()  

RHub

* devtools::check_rhub(interactive = F,
                       platforms   = c(rhub::platforms()$name),
                       env_vars    = c("_R_CHECK_DEPENDS_ONLY_"   = "true",
                                       "_R_CHECK_FORCE_SUGGESTS_" = "false"))
---  

## R CMD check results


```
devtools::check()  

0 errors ✓ | 0 warnings ✓ | 0 notes ✓
```

---  

## Reverse dependencies


**NONE**

```
tools::package_dependencies(packages = c('DGEobj.utils'),
                            db       = available.packages(), 
                            reverse  = TRUE)

$DGEobj.utils
character(0)
```
