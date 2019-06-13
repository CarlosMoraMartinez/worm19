
#' getGCProp
#' 
#' Calculates the GC content of a character sequence.
#' @param seq an input character sequence
#' @return The proportion of GCs (range 0-1).
#' @keywords GC, strings
#' @import Biostrings
#' @export

getGCProp <- function(seq){
  seq <- DNAString(seq)
  l <- length(seq)
  alp <- alphabetFrequency(seq)
  GC <- alp["G"] + alp["C"]
  return(GC/l)
}