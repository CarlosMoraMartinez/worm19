
#' globalMakeWindows
#'
#' Produces a data frame with a row for every region of max_width bearing 
#' nmotifs different motifs. As input it uses motifList, which is the output
#' of globalSearchMotifs() function.
#'
#' @param seqs dataframe with sequences. Must have "original.seq" and
#' "masked.seq" columns as strings and "seqID" with unique identifiers.
#' @param seqfile name of the set of sequences (the file from which seqs
#' have been read, without extension). Will be used to write output files.
#' @param motifList a list with data frames containing motif matches, as
#' returned by globalSearchMotifs().
#' @param PWMset list of PWMs
#' @param PWMsubsetname a string with the name of PWMsubset (used to
#' name files). Usually is the first letter of each motif name, tied
#' @param max_width maximum window span (usually 700 bp)
#' @param nmotifs number of different motifs needed to make a window
#' @param PWMsubset vector with names of PWMs to use in search.
#'
#' @return A dataframe with all the information regarding windows. It also
#' writes both a .csv and a .bed file with the genomic intervals.
#' @seealso globalSearchMotifs
#' @keywords windows, motifs, search, PWMs
#' @import Biostrings
#' @export

globalMakeWindows <- function(seqs, motifList, seqfile, PWMset, PWMsubsetname,
                            max_width, nmotifs = 6, PWMsubset = NULL){
  windows <- data.frame()
  # This outer loop is for genomic regions (upstream or introns)
  for(i in names(motifList)){
    region <- motifList[[i]]
    # This other loop is for the masked or unmasked
    for (j in names(region)){
      if(nrow(region[[j]]) > 0){
        this_seq_windows <- getWindows(region[[j]],
                                       max_width, nmotifs, PWMsubset)
        windows <- rbind(windows,
                        data.frame(seqID=rep(i, nrow(this_seq_windows)),
                                   this_seq_windows,
                                   seq.mask = rep(j, nrow(this_seq_windows))))
      }
    }
  } # windows calculated

  ## Calculate extra columns for windows
  if(nrow(windows) > 0){
    # windows<-read.table("all_windowsCoords1.csv",
    #                     header=TRUE, sep="\t", stringsAsFactors=FALSE)

    for (i in 1:nrow(seqs)){
      windows[windows$seqID == seqs[i, ]$seqID,
              "chrom_start"] <- windows[windows$seqID == seqs[i, ]
                                        $seqID,]$start + seqs[i, ]$start - 1

      windows[windows$seqID == seqs[i, ]$seqID,
              "chrom_end"] <- windows[windows$seqID == seqs[i, ]
                                      $seqID, ]$end + seqs[i, ]$start - 1

      windows[windows$seqID == seqs[i, ]$seqID,
              "strand"] <- rep(seqs[i, ]$strand,
                               nrow(windows[windows$seqID == seqs[i, ]$seqID, ]))

      windows[windows$seqID == seqs[i, ]$seqID,
              "chrom"] <- rep(seqs[i, ]$seqnames,
                              nrow(windows[windows$seqID == seqs[i, ]$seqID, ]))
    }

    windows[, "window_width"] <- windows$end - windows$start + 1
    windows[, "windowID"] <- paste(windows$seqID,
                                   1:nrow(windows), windows$seq.mask, sep="_")
    write.table(windows, file = paste(seqfile, PWMset, PWMsubsetname,
                                      max_width, "bp_winCoords.csv",
                                      sep="_", collapse="_"),
                row.names=FALSE, sep="\t")

    winfilename <- paste(seqfile, PWMset,PWMsubsetname, max_width,
                         "bp_win.bed", sep="_", collapse="_")
    windows_bed <- data.frame(chrom = windows$chrom,
                              chromStart = windows$chrom_start,
                              chromEnd = windows$chrom_end,
                              name = windows$windowID,
                              score = 0, strand = windows$strand)
    write.table(windows_bed, file = winfilename, sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
  }
  return(windows)
}
