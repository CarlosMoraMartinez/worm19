
#' getFiltFromPWMs
#' 
#' Given a set of PWMs, produces a set of strings representing each PWM 
#' (a regular expression allowing at each position all the nucleotides 
#' with probability > 0).
#' 
#' @param PWMs list of matrixes, with names
#' @return A vector of strings, with names.
#' @keywords filters, PWMs, strings
#' @export

getFiltFromPWMs <- function(PWMs){
  filt <- list()
  for(i in 1:length(PWMs)){
    aux <- character()
    for(j in 1:ncol(PWMs[[i]])){
      if(length(which(PWMs[[i]][, j] >0)) >1){ 
        aux[j] <- paste("[", paste(names(which(PWMs[[i]][, j] >0)), 
                                   sep="", collapse=""), "]", 
                        sep="", collapse="")  
      } else {
        aux[j] <- names(which(PWMs[[i]][, j] >0))
      }
    }
    filt[[i]] <- paste(aux, sep="", collapse="")
  } 
  
  filt <- unlist(filt)
  names(filt) <- names(PWMs)
  
  return(filt)
}