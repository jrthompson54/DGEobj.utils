#' Run edgeR normalization on DGEobj
#'
#' Takes a DGEobj and adds a normalized DGEList object representing the result of
#' edgeR normalization (calcNormFactors).
#'
#' @param dgeObj A DGEobj containing counts, design data, and gene annotation.
#' @param normMethod One of "TMM", "RLE", "upperquartile", or "none". (Default = "TMM")
#' @param itemName optional string represents the name of the new DGEList. It must be unique and not exist
#' in the passed DGEobj (Default = "DGEList")
#' @param includePlot Enable returning a "canvasXpress" or "ggplot" bar plot of the norm.factors
#' produced (Default = FALSE). Possible values to pass:
#'  \itemize{
#'   \item \strong{FALSE or NULL}: Disable plot
#'   \item \strong{TRUE or "canvasXpress"}: returns "canvasXpress" bar plot of the norm.factors produced.
#'   \item \strong{"ggplot"}: returns "ggplot" bar plot of the norm.factors produced.
#' }
#' @param plotLabels Sample text labels for the plot. Length must equal the number of
#'   samples. (Default = NULL; sample number will be displayed)
#'
#' @return A DGEobj with a normalized DGEList added or a list containing the normalized DGEobj and a plot
#'
#' @examples
#' \dontrun{
#'    # NOTE: Requires the edgeR package
#'
#'    myDGEobj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'    myDGEobj <- DGEobj::resetDGEobj(myDGEobj)
#'
#'    # Default TMM normalization
#'    myDGEobj <- runEdgeRNorm(myDGEobj)
#'
#'    # With some options and plot output
#'    require(canvasXpress)
#'    myDGEobj <- DGEobj::resetDGEobj(myDGEobj)
#'    obj_plus_plot <- runEdgeRNorm(myDGEobj,
#'                                  normMethod = "upperquartile",
#'                                  includePlot = TRUE)
#'    myDGEobj <- obj_plus_plot[[1]]
#'    obj_plus_plot[[2]]
#' }
#'
#' @importFrom dplyr %>%
#' @importFrom DGEobj addItem getItem
#' @importFrom assertthat assert_that
#'
#' @export
runEdgeRNorm <- function(dgeObj,
                         normMethod  = "TMM",
                         itemName    = "DGEList",
                         includePlot = FALSE,
                         plotLabels  = NULL) {
    assertthat::assert_that(requireNamespace("edgeR", quietly = TRUE),
                            msg = "edgeR package is required to apply edgeR normalization to the given DGEobj")

    assertthat::assert_that(!missing(dgeObj),
                            !is.null(dgeObj),
                            class(dgeObj) == "DGEobj",
                            msg = "dgeObj must be of class 'DGEobj'.")
    assertthat::assert_that(!is.null(normMethod),
                            is.character(normMethod),
                            length(normMethod) == 1,
                            tolower(normMethod) %in% c("tmm", "rle", "upperquartile", "none"),
                            msg = "normMethod must be only one of the following values 'TMM', 'RLE', 'upperquartile', 'none'.")
    assertthat::assert_that(!is.null(itemName),
                            !itemName %in% names(dgeObj),
                            length(itemName) == 1,
                            msg = "itemName must be a singular, unique and not NULL character value.")

    funArgs <- match.call()
    do.call("require", list("edgeR"))

    if (is.null(includePlot)) {
        plot_type <- "none"
    } else if (is.logical(includePlot) && length(includePlot) == 1) {
        plot_type <- ifelse(includePlot, "canvasxpress", "none")
    } else if (is.character(includePlot) && length(includePlot) == 1) {
        if (tolower(includePlot) %in% c("canvasxpress", "ggplot")) {
            plot_type <- tolower(includePlot)
        } else {
            warning("includePlot must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE.")
            plot_type <- "none"
        }
    } else {
        warning("includePlot must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE.")
        plot_type <- "none"
    }

    MyDGElist <- tryCatch({
        do.call("calcNormFactors",
                list(object = do.call("DGEList",
                                      list(counts = as.matrix(DGEobj::getItem(dgeObj, "counts")))),
                     method = normMethod))
    },
    error = function(e) {
        message("Unexpected error: ", e$message, " happened during edgeR normalization of counts")
        return(NULL)
    })

    # Capture the DGEList
    itemAttr <- list(normalization = normMethod)

    dgeObj   <- DGEobj::addItem(dgeObj,
                                item     = MyDGElist,
                                itemName = itemName,
                                itemType = "DGEList",
                                funArgs  = funArgs,
                                itemAttr = itemAttr,
                                parent   = "counts")
    if (plot_type != "none") {
        if (!is.null(plotLabels) && length(plotLabels) == ncol(dgeObj)) {
            labels <-  plotLabels
            angle  <-  45
        } else {
            if (!is.null(plotLabels) && length(plotLabels) != ncol(dgeObj)) {
                warning(paste("plotLabels must be a character vector with length equal to",
                              "the number of columns in dgeObj.  Assigning default values."))
            }
            labels <- 1:ncol(dgeObj)
            angle  <- ifelse(plot_type == "canvasxpress", 90, 0)
        }
        plot_data <- data.frame(row.names = factor(labels),
                                Norm.Factors = MyDGElist$samples$norm.factors)
    }

    plot <- NULL

    if (plot_type == "canvasxpress") {
        if ("canvasXpress" %in% .packages(all.available = T)) {
            do.call("require", list("canvasXpress"))

            plot <- do.call("canvasXpress",
                            list(data             = as.data.frame(t(plot_data)),
                                 graphOrientation = "vertical",
                                 graphType        = "Bar",
                                 showLegend       = FALSE,
                                 smpLabelRotate   = angle,
                                 smpTitle         = "Samples",
                                 theme            = "CanvasXpress",
                                 widthFactor      = 1.5,
                                 title            = "Normalization Factors",
                                 xAxisTitle       = "Norm Factors",
                                 color            = "dodgerblue3",
                                 afterRender      = list(list("sortSamples",
                                                              list(sortDir = "ascending"))),
                                 decorations      = list(line = list(list(value = 1,
                                                                          width = 2,
                                                                          color = "rgb(255,0,0)"))),
                                 setMinX          = 0))
        } else {
            message('The canvasXpress package is not available, unable to create plot.')
        }
    } else if (plot_type == "ggplot") {
        if ("ggplot2" %in% .packages(all.available = T)) {
            do.call("require", list("ggplot2"))

            # resolve notes in the package related to NSE: 'no visible global function definition'
            ggplot = aes = geom_bar = geom_hline = xlab = ylab = ggtitle = theme_bw = theme = element_text <- NULL

            Norm.Factors <- NULL
            plot <- ggplot(plot_data, aes(x = labels, y = Norm.Factors)) +
                geom_bar(stat  = "identity",
                         color = "dodgerblue4",
                         fill  = "dodgerblue3",
                         width = 0.7) +
                geom_hline(yintercept = 1.0, color = "red") +
                xlab("Samples") +
                ylab("Norm Factors") +
                ggtitle("Normalization Factors") +
                theme_bw(12) +
                theme(axis.text.x = element_text(angle = angle, hjust = 1.0))
        } else {
            message('The ggplot2 package is not available, unable to create plot.')
        }
    }

    if (plot_type != "none") {
        list(dgeObj = dgeObj, plot = plot)
    } else {
        dgeObj
    }
}

