
#' winToDataframe
#' 
#' Given a set of motifs that form a window, it returns a dataframe with 
#' the data for the window (i.e, pastes motif names, offsets, strands, 
#' starts, ends in a string and calculates start and end) of the window.
#' 
#' @param window_motifs a dataframe with motif data (start, end, 
#' offsets, motif_name, strand, width)
#' @keywords windows, motifs
#' @export

winToDataframe <- function(window_motifs){
  window <- data.frame(width = paste(window_motifs$width, 
                                     sep="_", collapse="_"), 
                      signature = paste(window_motifs$motif_name, 
                                        sep="_", collapse="_"), 
                      offsets = paste(window_motifs$offsets, 
                                      sep="_", collapse="_"), 
                      strands = paste(window_motifs$strand, 
                                      sep="_", collapse="_"),
                      seqs = paste(window_motifs$seq, 
                                   sep="_", collapse="_"), 
                      starts = paste(window_motifs$start, 
                                     sep="_", collapse="_"), 
                      ends = paste(window_motifs$end, 
                                   sep="_", collapse="_"), 
                      start = min(window_motifs$start),
                      end = max(window_motifs$end),
                      different_motifs = length(unique(window_motifs$motif_name)),
                      stringsAsFactors = FALSE)
  
  return(window)
}