#############################################################################
# Description: Plot the probability weighting function
# Notes: "auto" binsize tries to determine the best binsize to display the relationship b/w bias.data and %DE
# The "pwf" input can be either a list with the data.frame as an entry, or just the data.frame
# Author: Matthew Young
# Date Modified: 13/12/2010

plotPWF <- function(pwf, binsize = "auto", pwf_col = 3, pwf_lwd = 2, xlab = "Biased Data in <binsize> gene bins.", ylab = "Proportion DE", ...) {
  # We shouldn't try and plot NAs obviously...
  w <- !is.na(pwf$bias.data)
  o <- order(pwf$bias.data[w])
  # What is the total range in the fit?
  rang <- max(pwf$pwf, na.rm = TRUE) - min(pwf$pwf, na.rm = TRUE)
  if (rang == 0 & binsize == "auto") binsize <- 1000
  if (binsize == "auto") {
    # A low number of starting genes to bin, usually 100
    binsize <- max(1, min(100, floor(sum(w) * .08)))
    resid <- rang
    # Turn off warnings till we've worked out what we're doing
    oldwarn <- options()$warn
    options(warn = -1)
    # Keep increasing the number of genes in each bin until the scatter around the lines reaches the cutoff.
    # Stop if we reach only 10 bins for the entire plot
    while (binsize <= floor(sum(w) * .1) & resid / rang > .001) {
      binsize <- binsize + 100
      # Assign each gene a "bin number"
      splitter <- ceiling(1:length(pwf$DEgenes[w][o]) / binsize)
      # Determine the percentage DE in each bin
      de <- sapply(split(pwf$DEgenes[w][o], splitter), mean)
      # Determine the average length in each bin
      binlen <- sapply(split(as.numeric(pwf$bias.data[w][o]), splitter), mean)
      # Calculate the residuals, how much the binned data deviates from the PWF
      resid <- sum((de - approx(pwf$bias.data[w][o], pwf$pwf[w][o], binlen)$y)^2) / length(binlen)
    }
    options(warn = oldwarn)
  } else {
    # Assign each gene a "bin number"
    splitter <- ceiling(1:length(pwf$DEgenes[w][o]) / binsize)
    # Determine the percentage DE in each bin
    de <- sapply(split(pwf$DEgenes[w][o], splitter), mean)
    # Determine the average length in each bin
    binlen <- sapply(split(as.numeric(pwf$bias.data[w][o]), splitter), mean)
  }
  # Now we've settled on a binsize, plot it
  # Did the user specify the labels? If so we can't put in the defaults or they'll be used twice and errors result.
  xlab <- gsub("<binsize>", as.character(binsize), xlab)
  if ("xlab" %in% names(list(...))) {
    if ("ylab" %in% names(list(...))) {
      plot(binlen, de, ...)
    } else {
      plot(binlen, de, ylab = ylab, ...)
    }
  } else if ("ylab" %in% names(list(...))) {
    plot(binlen, de, xlab = xlab, ...)
  } else {
    plot(binlen, de, xlab = xlab, ylab = ylab, ...)
  }
  # Add the PWF
  lines(pwf$bias.data[w][o], pwf$pwf[w][o], col = pwf_col, lwd = pwf_lwd)
}
