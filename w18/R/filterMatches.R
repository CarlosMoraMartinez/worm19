
#' filterMatches
#' 
#' Filters sequences that do not match a pattern.
#' 
#' @param views an object of the class XStringViews
#' @param pattern a string of characters that contains the motif core
#' pattern (e.g.: "[AT]CAA") to filter the matrix alignments. Accepts "" 
#' as a pattern.
#' @return A dataframe with the matches (all the 'views' objects in it).
#' @keywords filters, sequences, matches
#' @import Biostrings
#' @export

filterMatches <- function(views, pattern = ""){
  matches <- data.frame(start = start(views), 
                        end = end(views), 
                        width = width(views), 
                        offsets = (start(views) + width(views)/2), 
                        seq = as.character(views))
  filtered_matches <- matches[grep(pattern=pattern, x=matches$seq), ]
  return(filtered_matches)
}
