#' Run a power analysis on counts and design matrix
#'
#' Take a counts matrix and design matrix and return a power analysis using
#' the RNASeqPower package. The counts matrix should be pre-filtered to remove
#' non-expressed genes using an appropriate filtering criteria. The design matrix
#' should describe the major sources of variation so the procedure can dial
#' out those known effects for the power calculations.
#'
#' Note, both 'RNASeqPower' and 'statmod' packages are required to run this function as follow:
#' \itemize{
#' \item {'RNASeqPower' package is required to run power analysis on the given counts matrix and design matrix.}
#' \item {'statmod' package is required to run estimate dispersion calculations}
#' }
#'
#' If includePlots = FALSE (the default) or NULL, the function will return a tall skinny dataframe
#' of power calculations for various requested combinations of N and significance
#' thresholds.
#'
#' If includePlots = TRUE, "canvasXpress" or "ggplot", a list is returned with an additional two
#' "canvasXpress" or ggplots (plots) to the dataframe.
#'
#' @param countsMatrix A counts matrix or dataframe of numeric data. (Required)
#' @param designMatrix A design matrix or dataframe of numeric data. (Required)
#' @param depth A set of depth to use in the calculations.  The default depths of
#'        c(10, 100, 1000) respectively represent a detection limit, below average
#'        expression, and median expression levels, expressed in read count units.
#' @param N A set of N value to report power for. (Default = c(3, 6, 10, 20))
#' @param FDR FDR thresholds to filter for for FDR vs. Power graph. (Default = c(0.05, 0.1))
#' @param effectSize A set of fold change values to test. (Default = c(1.2, 1.5, 2))
#' @param includePlots controls adding tow plots to the returned dataframe (Default = FALSE).
#'        The two plots are; a ROC curve (FDR vs. Power) and a plot of N vs. Power.
#'        Possible values to pass:
#'        \itemize{
#'          \item \strong{FALSE or NULL}: Disable plots
#'          \item \strong{TRUE or "canvasXpress"}: returns "canvasXpress" plots.
#'          \item \strong{"ggplot"}: returns "ggplot" plots.}
#'
#' @return a dataframe of power calculations or a list of the dataframe and defined plots as defined by the "includePlots" argument.
#'
#' @examples
#' if (requireNamespace("RNASeqPower", quietly = TRUE) &&
#'     requireNamespace("statmod", quietly = TRUE)) {
#'
#'     dgeObj <- readRDS(system.file("exampleObj.RDS", package = "DGEobj"))
#'     counts <- dgeObj$counts
#'     dm     <- DGEobj::getType(dgeObj, type = "designMatrix")[[1]]
#'
#'     resultList <- runPower(countsMatrix = counts,
#'                            designMatrix = dm,
#'                            includePlots = TRUE)
#'
#'     head(resultList[[1]]) # dataframe
#'     resultList[[2]]       # ROC Curves Plot
#'     resultList[[3]]       # N vs Power Plot
#' }
#'
#' @importFrom edgeR estimateDisp DGEList calcNormFactors aveLogCPM
#' @importFrom dplyr filter arrange select %>%
#' @importFrom stats approx power
#' @importFrom assertthat assert_that
#'
#' @export
runPower <- function(countsMatrix,
                     designMatrix,
                     depth = c(10, 100, 1000),
                     N = c(3, 6, 10, 20),
                     FDR = c(0.05, 0.1),
                     effectSize = c(1.2, 1.5, 2),
                     includePlots = FALSE) {
    assertthat::assert_that(requireNamespace("RNASeqPower", quietly = TRUE),
                            msg = "RNASeqPower package is required to run power analysis on the given counts matrix and design matrix.")
    assertthat::assert_that(requireNamespace("statmod", quietly = TRUE),
                            msg = "'statmod' package is required to run estimate dispersion calculations")
    assertthat::assert_that(!missing(countsMatrix),
                            !is.null(countsMatrix),
                            class(countsMatrix)[[1]] %in% c("matrix","data.frame"),
                            msg = "countsMatrix must be specified and must be of class matrix or dataframe.")
    assertthat::assert_that(!missing(designMatrix),
                            !is.null(designMatrix),
                            class(designMatrix)[[1]] %in% c("matrix","data.frame"),
                            msg = "designMatrix must be specified and must be of class matrix or dataframe.")
    if (any(is.null(depth),
            !is.numeric(depth),
            length(depth)  != 3)) {
        warning("depth must be a vector of 3 integer values. Assigning default values 10, 100, 1000.")
        depth  <-  c(10, 100, 1000)
    }

    if (any(is.null(N),
            !is.numeric(N),
            length(N)  != 4)) {
        warning("N must be a vector of 4 integer values. Assigning default values 3, 6, 10, 20.")
        N  <-  c(3, 6, 10, 20)
    }

    if (any(is.null(FDR),
            !is.numeric(FDR),
            length(FDR)  != 2)) {
        warning("FDR must be a vector of 2 integer values. Assigning default values 0.05, 0.1.")
        FDR  <-  c(0.05, 0.1)
    }

    if (any(is.null(effectSize),
            !is.numeric(effectSize),
            length(effectSize)  != 3)) {
        warning("effectiveSize must be a vector of 3 integer values. Assigning default values 1.2, 1.5, 2.")
        effectSize  <-  c(1.2, 1.5, 2)
    }
    # Fit the BCV data and define the BCV for each depth requested.
    # Estimate dispersion
    dgelist <- countsMatrix %>%
        as.matrix() %>%
        edgeR::DGEList() %>%
        edgeR::calcNormFactors() %>%
        edgeR::estimateDisp(design = designMatrix, robust = TRUE)

    # Get a fitted CV values for each input value of depth
    # BCV is the sqrt of Dispersion
    GeoMeanLibSize  <- dgelist$samples$lib.size %>% log %>% mean %>% exp
    depth_avelogcpm <- edgeR::aveLogCPM(depth, GeoMeanLibSize)
    depthBCV        <- sqrt(approx(dgelist$AveLogCPM, dgelist$trended.dispersion,
                            xout = depth_avelogcpm, rule = 2, ties = mean)$y)

    n <- seq(min(N),max(N),1)   # For the N vs P plot
    alpha <- seq(0.05, 0.9, 0.05)  # Alpha is FDR levels

    # Initialize a dataframe for the results table
    pdat <- data.frame(depth = double(),
                       n = double(),
                       effect = double(),
                       alpha = double(),
                       powerVal = double(),
                       stringsAsFactors = FALSE)
    for (D in depth) {
        cv <- depthBCV[D == depth]
        for (Nf in n) {
            for (E in effectSize) {
                for (A in alpha) {
                    do.call("require", list("RNASeqPower"))
                    P    <- do.call("rnapower", list(depth = D, n = Nf, cv = cv, effect = E, alpha = A))
                    pdat <- rbind(pdat, c(depth = D, n = Nf, effect = E, alpha = A, powerVal = P))
                }
            }
        }
    }
    colnames(pdat) <- c("depth", "n", "effect", "alpha", "power")
    if (is.null(includePlots)) {
        plot_type <- "none"
    } else if (is.logical(includePlots) && length(includePlots) == 1) {
        plot_type <- ifelse(includePlots, "canvasxpress", "none")
    } else if (is.character(includePlots) && length(includePlots) == 1) {
        if (tolower(includePlots) %in% c("canvasxpress", "ggplot")) {
            plot_type <- tolower(includePlots)
        } else {
            warning("includePlots must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE.")
            plot_type <- "none"
        }
    } else {
        warning("includePlots must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE.")
        plot_type <- "none"
    }
    rocdat <- dplyr::filter(pdat, n %in% N)
    rocdat$depth <- as.factor(rocdat$depth)
    # N vs Power
    # Filter to just a few FDR thresholds
    ndat <- dplyr::filter(pdat, alpha %in% FDR)
    ndat$depth <- as.factor(ndat$depth)
    ndat$FDR <- ndat$alpha

    result <- pdat

    if (plot_type == "canvasxpress") {
        if ("canvasXpress" %in% .packages(all.available = T)) {
            do.call("require", list("canvasXpress"))
            do.call("require", list("htmlwidgets"))

            rocdat   <- rocdat %>%
                dplyr::arrange(alpha)
            cx_data  <- rocdat %>%
                dplyr::select(alpha, power)
            var_data <- rocdat %>%
                dplyr::select(depth, n, effect)
            var_data$n      <- paste0("n:", var_data$n)
            var_data$effect <- paste0("effect: ", var_data$effect)
            events <- do.call("JS",
                              list("{'mousemove' : function(o, e, t) {
                                                       if (o != null && o != false) {
                                                           t.showInfoSpan(e, '<b>Alpha</b>: ' + o.y.data[0][0] +
                                                                             '<br><b>Power</b>: ' + o.y.data[0][1]);
                                                        };}}"))
            roc <- do.call("canvasXpress", list(data                 = cx_data,
                                              varAnnot             = var_data,
                                              segregateVariablesBy = list("effect", "n"),
                                              layoutType           = "rows",
                                              dataPointSize        = 5,
                                              spiderBy             = "depth",
                                              shapeBy              = "depth",
                                              colorBy              = "depth",
                                              title                = "ROC curves",
                                              xAxisTitle           = "FDR",
                                              yAxisTitle           = "Power",
                                              events               = events,
                                              afterRender          = list(list("switchNumericToString",
                                                                               list("depth",FALSE)))))
            ndat     <- ndat %>%
                dplyr::arrange(n)
            cx_data  <- ndat %>%
                dplyr::select(n, power)
            var_data <- ndat %>%
                dplyr::select(depth, FDR, effect)
            var_data$FDR    <- paste0("FDR:", var_data$FDR)
            var_data$effect <- paste0("effect: ", var_data$effect)
            events <- do.call("JS",
                              list("{'mousemove' : function(o, e, t) {
                                                       if (o != null && o != false) {
                                                           t.showInfoSpan(e, '<b>N</b>: ' + o.y.data[0][0] +
                                                                          '<br><b>Power</b>: ' + o.y.data[0][1]);
                                                        };}}"))
            NvP <- do.call("canvasXpress", list(data                 = cx_data,
                                              varAnnot             = var_data,
                                              segregateVariablesBy = list("FDR", "effect"),
                                              layoutType           = "rows",
                                              dataPointSize        = 5,
                                              spiderBy             = "depth",
                                              shapeBy              = "depth",
                                              colorBy              = "depth",
                                              title                = "N vs Power",
                                              xAxisTitle           = "N",
                                              yAxisTitle           = "Power",
                                              events               = events,
                                              afterRender          = list(list("switchNumericToString",
                                                                               list("depth",FALSE)))))

            result <- list(PowerData = pdat, ROC = roc, NvP = NvP)
        } else {
            message('The canvasXpress package is not available, unable to create plots.')
        }
    } else if (plot_type == "ggplot") {
        if ("ggplot2" %in% .packages(all.available = T)) {
            do.call("require", list("ggplot2"))

            # resolve notes in the package related to NSE: 'no visible global function definition'
            ggplot = aes = geom_line = scale_x_continuous = scale_y_continuous = facet_grid <- NULL
            label_both = ggtitle = xlab = ylab = expand_limits = theme = element_text = theme_gray <- NULL

            effect <- NULL
            roc <- ggplot(rocdat, aes(x = alpha, y = power, fill = depth, shape = depth, color = depth)) +
                geom_line(size = 1) +
                scale_x_continuous(breaks = seq(0, 1, 0.2)) +
                scale_y_continuous(breaks = seq(0, 1, 0.2)) +
                facet_grid(effect ~ n, labeller = label_both) +
                ggtitle("ROC curves") +
                xlab("\nFDR") +
                ylab("Power") +
                expand_limits(x = 0, y = 0) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                theme_gray(18)
            NvP <- ggplot(ndat, aes(x = n, y = power, fill = depth, shape = depth, color = depth)) +
                geom_line(size = 1) +
                scale_y_continuous(breaks = seq(0, 1, 0.2)) +
                facet_grid(FDR ~ effect, labeller = label_both) +
                ggtitle("N vs Power") +
                xlab("\nN") +
                ylab("Power") +
                expand_limits(x = 0, y = 0) +
                theme_gray()

            result <- list(PowerData = pdat, ROC = roc, NvP = NvP)
        } else {
            message('The canvasXpress package is not available, unable to create plots.')
        }
    }
    result
}
