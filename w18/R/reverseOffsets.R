
#' reverseOffsets
#' 
#' Calculates the offsets, starts and ends in the reverse complement
#' strand for the seqs of a given strand.
#' 
#' @param data dataframe with start, end, offsets, width and strand
#' @param len length of the original sequence where the subsequences in 
#' argument data were found
#' @param strand used to control it reverses only the negative strands
#' @return A data frame with the same dimensions as the original. Start 
#' keeps being before end.
#' @keywords strand, reverse, complement
#' @import Biostrings
#' @export

reverseOffsets <- function(data, len, strand = "-"){  
  data[data$strand == strand, 
       c("start", "end", "offsets")] <- sapply(data[data$strand == strand, 
                                                    c("end","start","offsets")], 
                                               FUN = function(x, len){
                                                 return(len-x+1)
                                                 }, 
                                               len = len)
  return(data)
}
