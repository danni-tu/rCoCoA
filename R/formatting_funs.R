#' Print new line in HTML files
#'
#' This function prints a new line in HTML files, which is useful when printing other things such as headers or plots.
#'
#' @examples
#' nl()
#'
#' @export
nl <- function(){
  cat("  \n")
  cat("  \n")
}

#' Round and format numbers for table output
#'
#' Converts a numeric value to a character value with the specified number of digits after the decimal.
#'
#' @param x a numeric value
#' @param digits an integer indicating the number of digits to display
#'
#' @return The formatted number (as a character).
#'
#' @examples
#' round(0.100, digits = 2)
#' round2(0.100, digits = 2)
#'
#' @keywords internal
#'

round2 <- function(x, digits){
  format(round(x, digits), nsmall = digits)
}

#' Scale vector using mean and standard deviation
#'
#' Linearly rescale a vector to have mean = 0 and sd = 1.
#'
#' @param x a numeric vector
#' @param na.rm logical value to be passed to mean() and sd()
#'
#' @return The scaled vector.
#'
#' @examples
#' scale2(1:10, na.rm = TRUE)
#'
#' @keywords internal
#'
scale2 <- function(x, na.rm = TRUE){
  (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
}

#' Output data.frame as pandoc.table()
#'
#' This is a wrapper around pandoc.table() that ensures all output will be printed on the same line.
#'
#' @param x a data.frame
#'
#' @examples
#' ptab(cars)
#'
#' @keywords internal
#'
ptab <- function(x){
  pander::pandoc.table(x, split.table = Inf)
}
