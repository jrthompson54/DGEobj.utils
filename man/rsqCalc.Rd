% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsqCalc.R
\name{rsqCalc}
\alias{rsqCalc}
\title{Calculate R-squared for each gene fit}
\usage{
rsqCalc(normMatrix, fit)
}
\arguments{
\item{normMatrix}{A normalized log2cpm matrix.}

\item{fit}{A MArrayLM object from limma::lmFit.}
}
\value{
A vector of R-squared values for each gene fit.
}
\description{
Takes a Log2CPM numeric matrix and MArrayLM fit object from limma::lmFit
and calculates R-squared for each gene fit.
}
\examples{
   dgeObj    <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
   log2cpm   <- convertCounts(dgeObj$counts, unit = "cpm", log=TRUE, normalize = "tmm")
   fitObject <- dgeObj$ReplicateGroupDesign_fit
   rsq       <- rsqCalc (log2cpm, fitObject)

}
