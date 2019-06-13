
#' filterseqsbyDist
#' 
#' Takes a dataframe with all the sequences and reduces them to a maximum 
#' distance from the first exon start. It also allows to filter out sequences
#' with length <= min_width.
#' 
#' @param seqs input sequences dataframe
#' @param dist integer vector of length 2; the first one is the distance
#' to the first exon allowed for the promoter (<= 0) and the second one 
#' is the maximum distance of introns (>= 0) 
#' @param min_width Minimum sequence length allowed (int). Defaults to 0.
#' @return A dataframe object equivalent to the input seqs.
#' @seealso trimSequencetoDist
#' @keywords filters, seqs, dist
#' @import stringr
#' @export

filterseqsbyDist <- function(seqs, dist, min_width = 0){
  seqs <- seqs[seqs$dist_to_first_exon > dist[1] & 
                 seqs$dist_to_first_exon < dist[2] & seqs$width >= min_width, ]
  
  if(nrow(seqs) > 0){
    for(i in 1:nrow(seqs)) {
      seqs[i, ] <- trimSequencetoDist(seqs[i, ], dist)
      }
    seqs$width <- seqs$end - seqs$start + 1
    seqs$masked_bp <- seqs$width - str_count(seqs$masked.seq, pattern="[ACTG]")
  }
  return(seqs)
}