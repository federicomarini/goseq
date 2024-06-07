#############################################################################
# Description: Generates the weighting curve that is used to generate the length corrected NULL distribution.
# A spline is fit to obtain a functional relationship between gene length and likelihood of differential exrpession
# Notes: By default genome and id are used to fetch length data from GeneLenDataBase, but the length of each gene can be supplied with bias.data
# Author: Matthew Young
# Date Modified: 20/12/2010


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
