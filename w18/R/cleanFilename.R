#' cleanFilename
#' 
#' Removes the extension and path in a file name f.
#' @param f Complete file path
#' @return Clean file name without "/" and the file extension.
#' @keywords trimming
#' @seealso trimextcsv
#' @export

cleanFilename <- function(f){
  f <- unlist(strsplit(f, split="/"))
  f <- f[length(f)]
  f <- trimextcsv(f)
  return(f)
}