#############################################################################
# Description: Generates the weighting curve that is used to generate the length corrected NULL distribution.
# A spline is fit to obtain a functional relationship between gene length and likelihood of differential exrpession
# Notes: By default genome and id are used to fetch length data from GeneLenDataBase, but the length of each gene can be supplied with bias.data
# Author: Matthew Young
# Date Modified: 20/12/2010




#' Probability Weighting Function
#' 
#' Calculates a Probability Weighting Function for a set of genes based on a
#' given set of biased data (usually gene length) and each genes status as
#' differentially expressed or not.
#' 
#' It is essential that the entire analysis pipeline, from summarizing raw
#' reads through to using \code{goseq} be done in just one gene identifier
#' format.  If your data is in a different format you will need to obtain the
#' gene lengths and supply them to the \code{nullp} function using the
#' \code{bias.data} argument.  Converting to a supported format from another
#' format should be avoided whenever possible as this will almost always result
#' in data loss.
#' 
#' \code{NA}s are allowed in the bias.data vector if you do not have
#' information about a certain gene.  Setting a gene to \code{NA} is preferable
#' to removing it from the analysis.
#' 
#' If \code{bias.data} is left as \code{NULL}, \code{nullp} attempts to use
#' \code{\link{getlength}} to fetch GO category to gene identifier mappings.
#' 
#' It is recommended you review the fit produced by the \code{nullp} function
#' before proceeding by leaving \code{plot.fit} as \code{TRUE}.
#' 
#' @param DEgenes A named binary vector where 1 represents DE, 0 not DE and the
#' names are gene IDs.
#' @param genome A string identifying the genome that \code{genes} refer to.
#' For a list of supported organisms run \code{\link{supportedGenomes}}.
#' @param id A string identifying the gene identifier used by \code{genes}.
#' For a list of supported gene IDs run \code{\link{supportedGeneIDs}}.
#' @param bias.data A numeric vector containing the data on which the DE may
#' depend.  Usually this is the median transcript length of each gene in bp.
#' If set to \code{NULL} \code{nullp} will attempt to fetch length using
#' \code{\link{getlength}}.
#' @param plot.fit Plot the PWF or not?  Calls \code{\link{plotPWF}} with
#' default values if \code{TRUE}.
#' @return A data frame with 3 columns, named "DEgenes", "bias.data" and "pwf"
#' with the rownames set to the gene names.  Each row corresponds to a gene
#' with the DEgenes column specifying if the gene is DE (1 for DE, 0 for not
#' DE), the bias.data column giving the numeric value of the DE bias being
#' accounted for (usually the gene length or number of counts) and the pwf
#' column giving the genes value on the probability weighting function. This
#' object is usually passed to \code{goseq} to calculate enriched categories or
#' \code{plotPWF} for further plotting.
#' 
#' @author Matthew D. Young \email{myoung@@wehi.edu.au}
#' 
#' @seealso \code{\link{supportedGenomes}}, \code{\link{supportedGeneIDs}},
#' \code{\link{goseq}}, \code{\link{getlength}}
#' 
#' @export
#' 
#' @references Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010)
#' \emph{Gene ontology analysis for RNA-seq: accounting for selection bias}
#' Genome Biology Date: Feb 2010 Vol: 11 Issue: 2 Pages: R14
#' 
#' @examples
#' data(genes)
#' pwf <- nullp(genes, 'hg19', 'ensGene')
#' 
nullp <- function(DEgenes, genome, id, bias.data = NULL, plot.fit = TRUE) {
  # Input Checking
  if (!is.null(bias.data) & length(bias.data) != length(DEgenes)) {
    stop("bias.data vector must have the same length as DEgenes vector!")
  }
  # Factors cause strange things to happen, remove them if they exist
  bias.data <- unfactor(bias.data)
  DEgenes <- unfactor(DEgenes)

  # Fetch length data from geneLenDataBase
  if (is.null(bias.data)) {
    bias.data <- getlength(names(DEgenes), genome, id)
  }

  # Fit a spline to the binary vector of DE calls vs gene length
  # May not have bias data for some of the entries, return NA at those positions
  pwf <- rep(NA, length(DEgenes))
  w <- !is.na(bias.data)
  pwf[w] <- makespline(bias.data[w], DEgenes[w])

  # Make a data frame which contains all the data used to make the fit and the fit itself
  out <- data.frame(DEgenes = DEgenes, bias.data = bias.data, pwf = pwf, stringsAsFactors = FALSE)
  rownames(out) <- names(DEgenes)

  # Plot the PWF if the arument has been specified
  if (plot.fit) {
    plotPWF(out)
  }
  out
}
