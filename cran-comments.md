## Comments from Maintainer

Resubmission Comments:
* Analysis tests (runXYZ function testing) skipped on CRAN to reduce testing time for cran
* Description rewritten to pass as many spelling checks as possible w/o industry-specific lingo
* DOI link fixed

Initial Submission Comments:
* This is a new package to be added to CRAN, it contains utilities for running DE analysis 
  in concert with the DGEobj package
* Code Coverage is 99%

---  

## Test environments

RStudio Server Pro (ubuntu 18.04.2)  

* R 3.6.3
* R 4.0.4

Travis-CI (ubuntu 16.04.6)

* R 3.6.3
* R 4.0.2
* R devel (2021-04-18 r80182)

WinBuilder

* devtools::check_win_devel()  
* devtools::check_win_release()  

RHub

* devtools::check_rhub(interactive = F)

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
revdepcheck::cran_revdeps('DGEobj.utils', bioc = T)

character(0)
```
