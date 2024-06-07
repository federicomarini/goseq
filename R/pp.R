pp <- function(total, count, i = i) {
  if (missing(count)) {
    count <- evalq(i, envir = parent.frame())
  }
  if (missing(total)) {
    total <- evalq(stop, envir = parent.frame())
  }
  message(round(100 * (count / total)), "%   \r")
}
