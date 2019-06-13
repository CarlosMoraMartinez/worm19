
#' readSeqFile
#' 
#' Reads a file with sequences and fills missing data to prevent 
#' further crashing.
#' @param seqfile the input sequences file
#' @return The same file without possible missing data.
#' @keywords sequences, strings
#' @export

readSeqFile <- function(seqfile){
  seqs <- read.table(seqfile, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  
  ## Eliminate sequences that have been manually detected to be wrong 
  ## (due to liftOver inconsistencies, etc.)
  # seqs<-seqs[!seqs$seqID %in% wrongseq, ]  
  
  seqs <- seqs[!duplicated(seqs[, names(seqs)!="seqID"]), ]
  if(!any(grepl("masked.seq", names(seqs)))){
    seqs$masked.seq <- "AAA"
  }
  if(!any(grepl("dist_to_first_exon", names(seqs)))){
    seqs$dist_to_first_exon <- 0
  }
  if(!any(grepl("region", names(seqs)))){
    seqs$region <- "intron"
  }
  if(!any(grepl("^gene$", names(seqs)))){
    seqs$gene <- seqs$seqID
  }
  return(seqs)
}