#' Read a whole file.
read_text_file <- function(path) {
  readChar(path, file.info(path)$size)
}

#' Read a Csnippet
read_Csnippet <- function(path) {
  library(pomp)
  Csnippet(read_text_file(path))
}
