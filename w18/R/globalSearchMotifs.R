
#' globalSearchMotifs
#' 
#' Searches a set of PWMs (including filtering) in a set of sequences 
#' and returns a list of lists, each one of them with a dataframe containing 
#' the matches in the original seq or the masked seq (masked repeats or 
#' other features). 
#' Be careful! before starting the search it will check whether a file exists 
#' with the same output name. If TRUE, it will load and won't perform
#' any search.
#' 
#' @param seqs dataframe with sequences. Must have "original.seq" and
#' "masked.seq" columns as strings and "seqID" with unique identifiers
#' @param PWMsubset vector with names of PWMs to use in search. 
#' @param mymotifs list of PWMs (matrixes). Must contain the PWMs in PWMsubset.
#' @param myfilt a string vector with regular expressions 
#' used to filter out PWM matches.
#' @param seqfile name of the set of sequences (the file from which seqs have
#' been read, without the extension). Will be used to write output files.
#' 
#' @return A list with an element for every row of the originl sequences 
#' dataframe, each one containing all the motif matches in both the original 
#' and masked sequence. It also saves the list in a .RData file.
#' @seealso getMotifs
#' @keywords motifs, search, PWMs
#' @import Biostrings
#' @export

globalSearchMotifs <- function(seqs, PWMsubset, mymotifs, myfilt, seqfile){
  motifList <- list()
  PWMsubsetname <- paste(unlist(sapply(PWMsubset, substr, 1, 1)), 
                         sep=",", collapse="")
  motifFile <- paste(seqfile, PWMset, PWMsubsetname,"mots.R", 
                     sep="_", collapse="_")
  
  if(!file.exists(motifFile)){      # If it not exists...
    for(i in 1:nrow(seqs)){
      this_seq_motifs <- list()
      this_seq_motifs[["original"]] <- getMotifs(seqs[i, ]$original.seq, 
                                                 mymotifs, myfilt, minscore)
      this_seq_motifs[["masked"]] <- getMotifs(seqs[i, ]$masked.seq, 
                                               mymotifs, myfilt, minscore)  
      motifList[[seqs[i, ]$seqID]] <- this_seq_motifs  
    }
    save(motifList, file = motifFile)
  } else {
    load(motifFile)     # If it exists, it just will be loaded
  }
  return(motifList)
}
