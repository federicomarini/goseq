#' Androgen stimulation of prostate cancer Cell lines.
#' 
#' This data set gives the RNA-seq data from an experiment measuring the
#' effects of androgen stimulation on prostate cancer.  Information is given
#' about all (ENSEMBL) genes for which there was at least one mapping read in
#' either the treated or untreated RNA-seq experiment.  The \code{edgeR}
#' package was used to determine which genes were differentially expressed.
#' The details of the analysis can be found in the \code{goseq} vignette.
#' 
#' 
#' @name genes
#' @docType data
#' @format A named vector of ENSEMBL genes, with 1 representing differential
#' expression.
#' @references Li, H., Lovci, M. T., Kwon, Y. S., Rosenfeld, M. G., Fu, X. D.,
#' Yeo, G. W. (2008) \emph{Determination of tag density required for digital
#' transcriptome analysis: application to an androgen-sensitive prostate cancer
#' model} Proceedings of the National Academy of Sciences of the United States
#' of America Date: Dec 23 Vol: 105 Issue: 51 Pages: 20179-84
#' @source Determination of tag density required for digital transcriptome
#' analysis: application to an androgen-sensitive prostate cancer model, 2008,
#' Li et. al.
#' @keywords datasets
#' @examples
#' 
#' data(genes)
#' head(genes)
#' 
NULL
