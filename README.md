# DGEobj.utils: A toolkit facilitating a limma/voom workflow Differential Gene Expression analysis

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/DGEobj.utils?color=9bc2cf)](https://cran.r-project.org/package=DGEobj.utils) 
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/grand-total/DGEobj.utils?color=9bc2cf)](https://cran.r-project.org/package=DGEobj.utils) 
[![Travis build status](https://app.travis-ci.com/cb4ds/DGEobj.utils.svg?branch=master)](https://app.travis-ci.com/cb4ds/DGEobj.utils?branch=master)
[![Codecov test coverage](https://app.codecov.io/gh/cb4ds/DGEobj.utils/branch/master/graph/badge.svg)](https://app.codecov.io/gh/cb4ds/DGEobj.utils?branch=master)
<!-- badges: end -->

This package implements a set of utility functions to enable a limma/voom workflow capturing
the results in the DGEobj data structure. Aside from implementing a well developed and popular
workflow in DGEobj format, the run* functions in the package illustrate how to wrap the
individual processing steps in a workflow in functions that capture important metadata,
processing parameters, and intermediate data items in the DGEobj data structure. This function-
based approach to utilizing the DGEobj data structure insures consistency among a collection of
projects processed by these methods and thus facilitates downstream automated meta-analysis.

### Functionality includes: 

#### Analysis

* **runContrasts**: Build contrast matrix and calculate contrast fits
* **runEdgeRNorm**: Run edgeR normalization on DGEobj
* **runIHW**: Apply Independent Hypothesis Weighting (IHW) to a list of topTable dataframes
* **runPower**: Run a power analysis on counts and design matrix
* **runQvalue**: Calculate and add q-value and lFDR to dataframe
* **runSVA**: Test for surrogate variables
* **runVoom**: Run functions in a typical voom/lmFit workflow

#### Utilities

* **convertCounts**: Convert count matrix to CPM, FPKM, FPK, or TPM
* **extractCol**: Extract a named column from a series of df or matrices
* **lowIntFilter**: Apply low intensity filters to a DGEobj
* **rsqCalc**: Calculate R-squared for each gene fit
* **summarizeSigCounts**: Summarize a contrast list
* **topTable.merge**: Merge specified topTable df cols
* **tpm.direct**: Convert countsMatrix and geneLength to TPM units
* **tpm.on.subset**: Calculate TPM for a subsetted DGEobj
