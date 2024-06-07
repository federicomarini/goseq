#############################################################################
# Description: Reverses the direction of a bidirectional mapping, removes an entries pointing to NULL
# Notes:
# Author: Matthew Young
# Date Modified: 13/12/2010

reversemapping <- function(map) {
  tmp <- unlist(map, use.names = FALSE)
  names(tmp) <- rep(names(map), times = as.numeric(summary(map)[, 1]))
  return(split(names(tmp), as.vector(tmp)))
}
