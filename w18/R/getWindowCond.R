
#' getWindowCond
#'
#' Given a set of motifs (within a genomic range), it calculates whether
#' they are enough to form a window or not.

#' @param window_motifs data frame with motifs. Column "motif_name"
#' is compared to PWMsubset
#' @param nkinds minimum number of different motifs that a window needs
#' to have (in case PWMsubset is null). Defaults to 6.
#' @param PWMsubset character vector with motifs that are necessary to
#' form a window
#' @return If aa set of motifs is given (PWMsubset), returns TRUE if all
#' the motifs names in PWMsubset are found in window_motifs. If, otherwise,
#' PWMsubset is NULL, it returns TRUE if there are at least nkinds
#' kinds of motifs.
#' @keywords windows, motifs
#' @export

getWindowCond <- function(window_motifs, nkinds=6, PWMsubset){
  if(is.null(PWMsubset)){
    cond <- length(unique(window_motifs$motif_name)) >= nkinds
    
    ## This is changed in the CIPF's cluster copy to == nkinds
    ## to make the windows of only a specific number of motifs!!!

  } else {
    cond <- all(PWMsubset %in% unique(window_motifs$motif_name))
  }
  return(cond)
}
