#' The 'CSTWAS' package.
#'
#' @docType package
#' @name package
#'
#' @import tidyverse
#' @import ggplot2
#' @import ggrepel
#' @import RColorBrewer
#' @import progress
#' @import dplyr
#' @import ggVennDiagram
#' @import IRanges
#' @import EnsDb.Hsapiens.v79
#' @importFrom MASS mvrnorm
#' @importFrom dplyr select
#' @importFrom dplyr rename
#'
#' @description Transcriptome-wide association study (TWAS) is introduced to identify
#' significant expression-trait associations through imputations. It has been widely
#' used to analyze tissue-specific associations with the reference expression
#' quantitative trait loci (eQTL) panel. To increase the statistical power of TWAS results,
#' meta-analysis methods aggregating TWAS results across multiple tissues are developed.
#' However, most existing meta-analysis methods lose interpretation of disease etiology and
#' have limited power to identify weaker associations when only a few tissues are weakly
#' activated. Therefore, we developed the cross-tissue subset-based meta-analysis method,
#' so called CSTWAS.In this package, we aggregate the TWAS results across tissues and
#' perform meta-analysis through the subset-based test. R functions are provided for researchers
#' to integrate TWAS results across multiple tissues and visualize the results.
#'
#'
NULL
