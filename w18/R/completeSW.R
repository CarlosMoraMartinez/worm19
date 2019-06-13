
#' completeSW
#' 
#' Calls the function reduceWin, passing as arguments the superwindows 
#' and windows within 1 sequence each time, so that superwindow data 
#' are filled. Returns a data frame with all superwindows with complete data
#' (start, end, chrom, mask, pasted string with windowIDs).
#'  
#' @param sw dataframe with all superwindows
#' @param win dataframe with all windows
#' @seealso reduceWin
#' @keywords windows, superwindows
#' @export

completeSW <- function(sw, win){
  fullsw <- sw
  for(i in unique(sw$seqID)){
    this_sw <- sw[sw$seqID == i, ]
    this_w <- win[win$seqID == i, ]
    reducedwin <- reduceWin(this_sw, this_w)
    for(j in 1:ncol(reducedwin))
      fullsw[sw$seqID == i, names(reducedwin)[j]] <- reducedwin[, j]
  }
  return(fullsw)
}