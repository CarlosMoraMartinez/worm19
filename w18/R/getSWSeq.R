
#' getSWSeq
#' 
#' Given a sequence and 2 coordinates, returns the subsequence within 
#' the two coordinates, masked or unmasked, depending on the parameter "mask".
#' 
#' @param start interval start
#' @param end interval end
#' @param mask "original" if not masked, otherwise the masked sequence 
#' will be used
#' @param originalSeq raw DNA sequence (character string)
#' @param maskedSeq masked DNA sequence (usually repeats have been 
#' replaced by '+') (character string)
#' @return A character string with the subsequence found.
#' @keywords superwindows, strings
#' @export

getSWSeq <- function(start, end, mask, originalseq, maskedseq){
  if (as.character(mask) == "original"){
    seq <- originalseq
  } else {
    seq <- maskedseq
  }
  return(substr(seq, start, end))
}