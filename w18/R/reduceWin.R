
#' reduceWin
#' 
#' Takes data of several overlapping windows and calculates data for 
#' the superwindow: start, end, chromosome, seq.mask and a string with 
#' pasted windowIDs.
#' 
#' @param sw the superwindow dataframe with superwindows sharing the same seqID
#' @param win the windows dataframe with windows sharing the same seqID
#' @return A dataframe with all the superwindow data.
#' @keywords windows, superwindows
#' @export 

reduceWin <- function(sw, win){
  swdata <- data.frame()
  ## For each of the superwindows within a sequence,
  ## calculate data from windows
  for(i in 1:nrow(sw)){
    this_win <- win[win$chrom_start >= sw[i, ]$chrom_start
                    & win$chrom_end <= sw[i, ]$chrom_end, ]
    mergedmotifs <- mergeMotifs(this_win)
    for(j in 1:length(mergedmotifs)){
      swdata[i, names(mergedmotifs)[j]] <- mergedmotifs[j]
    }
    swdata[i, "start"] <- min(this_win$start)
    swdata[i, "end"] <- max(this_win$end)
    swdata[i, "seq.mask"] <- paste(unique(this_win$seq.mask), 
                                   sep="", collapse="")
    
    ## This will be slower but will ensure that we are doing it right
    swdata[i, "chrom"] <- paste(unique(this_win$chrom), sep="", collapse="")
    swdata[i, "windowIDs"] <- paste(this_win$windowID, sep=",", collapse=",")
  }
  return(swdata)
}