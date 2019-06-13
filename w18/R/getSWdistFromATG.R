
#' getSWdistFromATG
#' 
#' Calculates superwindows distance to first exon (not necessarily ATG) 
#' counting from the end of the superwindow (not the beginning!). 
#' 
#' @param start interval start within original sequence 
#' (not chromosomal coordinate)
#' @param end interval end within original sequence
#' @param strand interval strand (equal to sequence strand)
#' @param region "upstream" or "inton"
#' @param seqdist distance of original sequence to first exon start
#' @param seqwidth original sequence width
#' @return An integer number with the distance (in bp).
#' @keywords superwindows, distance
#' @export

getSWdistFromATG <- function(start, end, strand, region, seqdist, seqwidth){
  dist = NA
  if(strand == "+" & region == "upstream"){
    dist = seqdist - (seqwidth - start)
  } else if (strand == "-" & region == "upstream"){
    dist = seqdist - end + 1
  } else if (strand == "+" & region == "intron"){
    dist = seqdist + end - 1
  } else if (strand == "-" & region == "intron"){
    dist = seqdist + (seqwidth - start)
  }
  return (dist)
}