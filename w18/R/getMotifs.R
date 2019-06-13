
#' getMotifs
#' 
#' Searches motifs in a string of characters, filters the ones that 
#' don't have the core (myfilt).
#' 
#' @param seqstr the sequence where the motifs are to be aligned
#' @param mymotifs a list of PWMatrixes (the ones that are to be aligned 
#' to seq), each with a name
#' @param myfilt a character vector with the motif core pattern 
#' (e.g "[AT]CAA") to filter the matrix alignments
#' @param minscore is the min.score argument for Biostrings::matchPWM function
#' @return A dataframe with all the matches present in both the + and - strand
#' @seealso filterMatches, reverseOffsets
#' @keywords motifs, filters, strings 
#' @import Biostrings
#' @export

getMotifs <- function(seqstr, mymotifs, myfilt, minscore){
  seq <- DNAString(seqstr)
  seqrv <- reverseComplement(seq)
  result <- data.frame()
  
  for (i in 1:length(mymotifs)){  
    mot_name <- names(mymotifs)[i]
    mat <- mymotifs[[i]]
    
    # Align PWM to the sequence as it is (forward strand) 
    # and introduce matches in data.frame
    fwmatch <- filterMatches(matchPWM(pwm=mat, 
                                      subject=seq, 
                                      min.score=minscore, 
                                      with.score=TRUE), 
                             myfilt[[mot_name]])
    l <- nrow(fwmatch)
    fwmatch <- data.frame(
      fwmatch, 
      strand = rep("+", l), 
      motif_name = rep(mot_name, l),
      minscore = rep(minscore, l), 
      stringsAsFactors = FALSE)
    
    # Align PWM to the Reverse Complement of the sequence
    # and introduce matches in data.frame
    rvmatch <- filterMatches(matchPWM(pwm=mat, 
                                      subject=seqrv, 
                                      min.score=minscore, 
                                      with.score=TRUE), 
                             myfilt[[mot_name]])
    l <- nrow(rvmatch)
    rvmatch <- data.frame(
      rvmatch,
      strand = rep("-", l), 
      motif_name = rep(mot_name, l),
      minscore = rep(minscore, l), 
      stringsAsFactors = FALSE)
    
    result <- rbind(result, rbind(fwmatch, rvmatch))
  }
  
  # Recalculate the within-sequence coordinates (indexes) of the matches 
  # that are in the Reverse Complement sequence
  result <- reverseOffsets(result, nchar(seq), "-")
  return (result)   
}
