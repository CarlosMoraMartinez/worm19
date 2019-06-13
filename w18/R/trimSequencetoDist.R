
#' trimSequencetoDist
#' 
#' Gets the region of a sequence that is closer to the first exon than a 
#' given limit. Calculates start, end, and sequence (masked and original).
#' 
#' @param seq a row of a dataframe with the sequence data
#' @param dist integer vector of length 2. The first one number is the distance
#' to the first exon allowed for the promoter (<= 0) and the second one is 
#' the maximum distance of introns (>= 0).
#' @return A dataframe row with the same format as seq, with trimmed data.
#' @keywords trimming, trim, sequence, distance
#' @import stringr
#' @export

trimSequencetoDist <- function(seq, dist){
  ## Calculate start and end
  if(seq$region == "upstream"){
    if(seq$width + abs(seq$dist_to_first_exon) <= abs(dist[1])){
      return(seq)
    } else {
      if(seq$strand == "+"){
        seq$start <- seq$end + (dist[1] - seq$dist_to_first_exon)
        seq$original.seq <- substr(seq$original.seq, 
                                   nchar(seq$original.seq) - (seq$end-seq$start), 
                                   nchar(seq$original.seq))
        seq$masked.seq <- substr(seq$masked.seq, 
                                 nchar(seq$original.seq) - (seq$end-seq$start), 
                                 nchar(seq$original.seq))
      } else if (seq$strand == "-"){
        seq$end <- seq$start + (abs(dist[1]) + seq$dist_to_first_exon)
        seq$original.seq <- substr(seq$original.seq, 1, 
                                   (seq$end - seq$start + 1))
        seq$masked.seq <- substr(seq$original.seq, 1, 
                                 (seq$end - seq$start + 1))
      }
    }
  } else if (seq$region == "intron"){
    if(seq$width + seq$dist_to_first_exon <= dist[2]) return(seq)
    else{
      if(seq$strand == "+"){
        seq$end = seq$start + (dist[2] - seq$dist_to_first_exon)
        seq$original.seq <- substr(seq$original.seq, 1, 
                                   (seq$end - seq$start + 1))
        seq$masked.seq <- substr(seq$original.seq, 1, 
                                 (seq$end - seq$start + 1))
      }else if(seq$strand == "-"){
        seq$start = seq$end - (dist[2] - seq$dist_to_first_exon)
        seq$original.seq <- substr(seq$original.seq, 
                                   nchar(seq$original.seq) - (seq$end-seq$start), 
                                   nchar(seq$original.seq))
        seq$masked.seq <- substr(seq$masked.seq, 
                                 nchar(seq$original.seq) - (seq$end-seq$start), 
                                 nchar(seq$original.seq))
      }
    }
  }
  return(seq)
}