
#' mergeMotifs
#' 
#' Merges all the motif data from several windows into one row, avoiding
#' redundancies.
#' 
#' @param win window dataframe with the following fields:
#' "width", "signature", "offsets", "strands", "seqs", "starts", "ends".
#' @return A dataframe with 1 row, containing the mentioned variables. 
#' Each field contains the all the motifs pasted into a single string, 
#' using "_" as separators.
#' @keywords merge, motifs
#' @import data.table
#' @export

mergeMotifs <- function(win){
  win2 <- win[, c("width", "signature", "offsets", "strands", "seqs",
                  "starts", "ends")]
  winlist <- apply(win2, 1, FUN = function(x){
    as.data.frame(sapply(x, strsplit, split="_"))
  })
  windf <- rbindlist(winlist) # requires data.table package
  windf <- windf[!duplicated(windf), ]
  row <- sapply(windf, paste, sep="_", collapse="_")
  names(row) <- c("widths", "signature", "offsets",
                 "strands", "seqs", "starts", "ends")
  return(row)
}
