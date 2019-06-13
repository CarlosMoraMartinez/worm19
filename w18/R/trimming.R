
#' trim
#' 
#' Trims whitespaces from the end and the start of a string.
#' @keywords trim, trimming
#' @export

trim <- function (x) gsub("^\\s+|\\s+$", "", x)


#' trim2
#' 
#' Trims whitespaces from the end and the start of a string.
#' @keywords trim, trim2, trimming
#' @export

trim2 <- function (x) gsub("^_+|_+$", "", x)


#' trimextcsv
#' 
#' Trims ".csv" from the end of .csv files.
#' @keywords trim, trimextcsv, trimming
#' @export

trimextcsv <- function (x) gsub(".csv", "", x)