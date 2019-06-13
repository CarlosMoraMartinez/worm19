
#' makeSuperwindows
#' 
#' Reads a windows file and calls all the functions that are required 
#' to make superwindows. As an output, writes superwindows files, both with 
#' and without sequences.
#' 
#' @param filename name of the windows file
#' @param seqs dataframe with the original sequences analysed
#' @return Nothing, writes superwindows files directly.
#' @seealso cleanFilename, getSuperwindows, completeSW, reduceWin, writeBedFile
#' @keywords superwindows, windows, bed
#' @import data.table, GenomicRanges
#' @export

makeSuperwindows <- function(filename ="", seqs, 
                             outfolder = "./superwindows/", win = NULL){
  if (!file.exists(gsub("\\/$" , "", outfolder))){
    dir.create(outfolder)
  }
  outname <- paste(outfolder, "superwin_", cleanFilename(filename), 
                   sep="", collapse="")
  if(is.null(win)){
    win <- read.table(filename, header = TRUE, 
                      sep = "\t", stringsAsFactors = FALSE)
  }
  winoriginal <- win[win$seq.mask == "original", ]
  winmasked <- win[win$seq.mask == "masked", ]
  
  ## Make superwindows (only coordinates)
  sworiginal <- getSuperwindows(winoriginal)
  if(nrow(winmasked) > 0){
    swmasked <- getSuperwindows(winmasked)
  }
  
  ## Get the rest of the window data
  sworiginal2 <- completeSW(sworiginal, winoriginal)
  if(nrow(winmasked) > 0){
    swmasked2 <- completeSW(swmasked, winmasked)
  }
  
  ## Join superwindows from masked and original seqs. Create superwindows IDs
  if(nrow(winmasked) > 0){
    allsuperwindows <- rbind(sworiginal2, swmasked2)
    } else {
      allsuperwindows <- sworiginal2
    }
  allsuperwindows$seqID <- as.character(allsuperwindows$seqID)
  allsuperwindows[, "swID"] <- paste(allsuperwindows$seqID,
                                     1:nrow(allsuperwindows),
                                     "sw", allsuperwindows$seq.mask, sep="_")
  
  ## Write .bed file and superwindows file
  write.table(allsuperwindows, file = paste(outname, ".csv", 
                                            sep="", collapse=""), 
              row.names=FALSE, quote=FALSE, sep="\t")
  writeBedFile(allsuperwindows, allsuperwindows$swID, outname)
  
  ## Get sequences, get dist to first exon and calculate number 
  ## of non-masked nucleotides
  seqIndexes <- sapply(allsuperwindows$seqID, 
                       FUN = function(x, y){which(y %in% x)}, y = seqs$seqID)
  allsuperwindows[, "windows.seq"] <- mapply(FUN = getSWSeq,
                                           allsuperwindows$start, 
                                           allsuperwindows$end,
                                           allsuperwindows$seq.mask, 
                                           seqs[seqIndexes, ]$original.seq,
                                           seqs[seqIndexes, ]$masked.seq)
  
  allsuperwindows[, "dist_to_first_exon"] <- mapply(FUN = getSWdistFromATG,
                                                   allsuperwindows$start, 
                                                   allsuperwindows$end, 
                                                   allsuperwindows$strand,
                                                   seqs[seqIndexes, ]$region, 
                                                   seqs[seqIndexes, ]$dist_to_first_exon,
                                                   seqs[seqIndexes, ]$width)
  
  allsuperwindows[, "width_withoutreps"]<- str_count(allsuperwindows$windows.seq,
                                                     pattern="[ACGT]")
  
  ## Write superwindows file with all data
  swfilename <- paste(outname, "_withSeqs.csv", sep="", collapse="")
  write.table(allsuperwindows, file = swfilename, row.names=FALSE, 
              quote=FALSE, sep="\t")
  return(swfilename)
}