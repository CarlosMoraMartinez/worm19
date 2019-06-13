
#' globalMakeWindowsPar
#' 
#' Splits the data in batches and calls globalMakeWindows in parallel.
#' 
#' @param seqs data frame with sequences. Must have "original.seq" and 
#' "masked.seq" fields with strings and "seqID" with unique identifiers.
#' @param seqfile name of the set of sequences (the file from which seqs 
#' have been read without extension). Will be used to write output files.
#' @param motifList a list with data frames containing motif matches, as 
#' returned by globalSearchMotifs().
#' @param PWMset list of PWMs
#' @param PWMsubsetname a string with the name of the PWMsubset 
#' (used to name files), usually tie first letter of each motif name
#' @param max_width maximum window span (usually 700 bp)
#' @param nmotifs number of different motifs needed to make a window (6 
#' by default)
#' @param PWMsubset vector with names of PWMs to use in search.
#' @param n_partitions number of batches to split the data (24 by default)
#' 
#' @return A dataframe with all the information regarding windows. It also
#' writes a .csv and a .bed file with the genomic intervals.
#' @seealso getPartitions, globalMakeWindows, globalSearchMotifs
#' @keywords windows, motifs, search, PWMs, parallel
#' @import Biostrings, doParallel, foreach
#' @export

globalMakeWindowsPar <- function(seqs, motifList, seqfile, 
                                  PWMset,PWMsubsetname, max_width, 
                                  nmotifs = 6, PWMsubset = NULL, 
                                  n_partitions = 24){
  packs <- search()
  packs <- packs[grep(packs, pattern="package:")]
  packs <- sapply(packs, FUN=function(s){unlist(strsplit(s, split=":"))[2]})
  
  partitions <- getPartitions(length(motifList), n_partitions)
  windows <-foreach(numfile = 1:n_partitions, 
                    .packages = packs, .combine = rbind, 
                    .export=ls(), .verbose=TRUE) %dopar% {
    part_names <- names(motifList)[partitions[[numfile]]]
    motifListPart <- motifList[part_names]
    seqsPart <- seqs[seqs$seqID %in% part_names, ]
    seqfilePart <- paste(seqfile, '_part', numfile, sep="", collapse="")
    win <- globalMakeWindows(seqsPart, motifListPart, seqfilePart, 
                             PWMset, PWMsubsetname, max_width, 
                             min(window_min_motifs))
    win
  }
  write.table(windows, file = paste(seqfile, PWMset, 
                                    PWMsubsetname, max_width,
                                    "bp_winCoords_ALL.csv", 
                                    sep="_", collapse="_"), 
              row.names=FALSE, sep="\t")
  return(windows)
}
