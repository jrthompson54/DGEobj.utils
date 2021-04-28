#' DGEobj.utils Package Overview
#'
#' This package implements a set of utility functions to enable a limma/voom workflow capturing
#' the results in DGEobj data structure. Aside from implementing a well developed and popular
#' workflow in DGEobj format, the run* functions in the package illustrate how to wrap the
#' individual processing steps in a workflow in functions that capture important metadata,
#' processing parameters, and intermediate data items in the DGEobj data structure. This function-
#' based approach to utilizing the DGEobj data structure insures consistency among a collection of
#' projects processed by these methods and thus facilitates downstream automated meta-analysis.
#'
#' @section More Information:
#' \code{browseVignettes(package = 'DGEobj.utils')}
#'
#' @docType package
#' @name DGEobj.utils-package

NULL


.onLoad <- function(libname, pkgname) {
}

.onAttach <- function(libname, pkgname) {
}
