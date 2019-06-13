
#' getWindows 
#'
#' Given a dataframe with within-sequence coordinates for elements
#' of different kinds, it gets all the windows of a maximum length
#' maxwidth with at least nkinds different kinds of elements.
#'
#' @param matches data.frame with start, end, width, motif_names,
#' strand, offsets
#' @param maxwidth maximum width of the window
#' @param nkinds number of different kinds of motif_names that need
#' to be in a window
#' @return A dataframe with the window data (start, end, and
#' pasted motif information)
#' @keywords windows, motifs
#' @export

getWindows <- function(matches, maxwidth = 800,
                       nkinds = 8, PWMsubset = NULL){
  ## Do search centered in each motif
  windows <- data.frame()
  winstart <- min(matches$offsets)
  winend <- min(max(matches$offsets), winstart + maxwidth)
  window_motifs <- matches[matches$offsets >= winstart &
                             matches$offsets <= winend, ]

  if(getWindowCond(window_motifs, nkinds, PWMsubset)){
    windows <- rbind(windows, winToDataframe(window_motifs))
  }
  while((winend - winstart) == maxwidth){
    matches <- matches[matches$offsets != winstart, ]
    winstart <- min(matches$offsets)
    winend <- min(max(matches$offsets), winstart + maxwidth)
    window_motifs <- matches[matches$offsets >= winstart &
                               matches$offsets <= winend, ]
    if(getWindowCond(window_motifs, nkinds, PWMsubset)){
      windows <- rbind(windows, winToDataframe(window_motifs))
    }
  }
  return (windows)
}
