
#' makeSuperwindowsPar
#' 
#' Reads a windows file and calls all the functions that are required 
#' to make superwindows (in parallel). As an output, writes superwindows 
#' files, both with and without sequences.
#' 
#' @param filename name of the windows file
#' @param seqs dataframe with the original sequences analysed
#' @param n_partitions number of threads to use in parallel 
#' @return Nothing, writes superwindows files directly.
#' @seealso cleanFilename, getSuperwindows, completeSW, reduceWin, writeBedFile
#' @keywords superwindows, windows, bed
#' @import data.table, GenomicRanges, doParallel, foreach
#' @export

makeSuperwindowsPar <- function(filename= "", seqs, 
                                outfolder = "./superwindows/", win = NULL, 
                                n_partitions){
  packs <- search()
  packs <- packs[grep(packs, pattern="package:")]
  packs <- sapply(packs, FUN = function(s){unlist(strsplit(s, split=":"))[2]})
  seqnames <- unique(seqs$seqID)
  partitions <- getPartitions(length(seqnames), n_partitions)
  swfiles <- foreach(numfile = 1:n_partitions, 
                     .packages = packs, .combine = list, .export = ls(), 
                     .verbose=TRUE) %dopar% {
    part_names <- seqnames[partitions[[numfile]]]
    winPart <- win[win$seqID %in% part_names, ]
    seqsPart <- seqs[seqs$seqID %in% part_names, ]
    filePart <- paste(filename, '_part', numfile, sep = "", collapse = "")
    swfile <- makeSuperwindows(filePart, seqs = seqs, 
                               outfolder = outfolder, win = winPart)
    swfile
  }
  return (swfiles)
}
