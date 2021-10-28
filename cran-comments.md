## Comments from Maintainer

* Resolved CRAN check notes and errors
* Made some analysis packages suggested/optional and only required when using methods requiring them

---  

## Test environments

RStudio Server Pro (ubuntu 18.04.2)  

* R 3.6.3
* R 4.0.5
* R 4.1.1

Travis-CI (ubuntu 16.04.6)

* R 3.6.3
* R 4.0.2
* R devel (2021-09-29 r80990)

WinBuilder

* devtools::check_win_devel()  
* devtools::check_win_release()  

RHub

* devtools::check_rhub(interactive = F,
                       env_vars    = c("_R_CHECK_FORCE_SUGGESTS_" = "false"))

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
