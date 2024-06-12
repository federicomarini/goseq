#' goseq 
#'
#' Detects Gene Ontology and/or other user defined categories which are 
#' over/under represented in RNA-seq data.
#'
#' @importFrom mgcv gam smoothCon s mono.con pcls Predict.matrix
#' @importFrom graphics lines plot
#' @importFrom stats approx binomial dhyper median phyper runif
#' @importFrom utils data installed.packages tail
#' @importFrom AnnotationDbi as.list select
#' @importFrom GO.db GO.db
#' @importFrom GenomicFeatures transcriptLengths
#' @importFrom GenomeInfoDb `genome<-`
#' @importFrom methods show
#' @importFrom rtracklayer browserSession ucscGenomes ucscTableQuery getTable
#' @importFrom BiocGenerics relist toTable
#' @importFrom BiasedUrn dWNCHypergeo pWNCHypergeo
#' @import geneLenDataBase
#'
#' @name goseq-pkg
#' @docType package
#' @keywords internal
"_PACKAGE"


