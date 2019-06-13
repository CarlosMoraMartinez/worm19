
#' getSuperwindows
#' 
#' Given a data frame object with window data, it creates a GenomicRanges 
#' object and uses the "reduce" function to get non-redundant intervals
#' (i.e., superwindows). Then, it converts these data to dataframe again.
#' 
#' @param win dataframe with window data (must have at least "seqID", 
#' "chrom_start", "chrom_end", "strand" columns)
#' @return A dataframe with the superwindow coordinates.
#' @seealso grangesFromDF
#' @keywords windows, superwindows
#' @import GenomicRanges, data.table
#' @export

getSuperwindows <- function(win){
  # Make GRanges object
  grwin <- grangesFromDF(win)
  # Use reduce to calculate superwindows
  sw <- list()
  for(i in unique(seqnames(grwin))){
    sw[[i]] <- reduce(grwin[seqnames(grwin) == i])
  }
  
  ## Bind all in a new data frame
  # rbindlist in library(data.table) seems to be the most efficient 
  # solution for binding data frames
  swdf <- rbindlist(lapply(sw, as.data.frame))
  names(swdf) <- c("seqID", "chrom_start", "chrom_end", "width", "strand")
  return(swdf)
}