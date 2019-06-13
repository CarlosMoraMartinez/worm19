
#' getRandomPWMsAndFilts
#' 
#' Given a set of PWMs and filter (eg. [AT]TTG[CG]A), returns n lists of PWMs
#' with shuffled columns and its corresponding filters.
#' 
#' @param PWMs list of matrixes, with names
#' @param filt list of filters (regular expression-like, eg. [AT]TTG[CG]A) 
#' @param n number of permutated lists (defaults to n=10).
#' @return n lists of PWMs with shuffled columns (1 shuffled PWM per 
#' originalPWM) and its corresponding filters (a list of lists, with n 
#' permutations of each PWM.)
#' @keywords PWMs, filters
#' @export

getRandomPWMsAndFilts <- function(PWMs, filt, n=10){
  newfilt <- filt
  mats <- PWMs
  ranBlock <- list()
  ranMatFilt <- list()
  for(i in 1:n){
    for(m in names(PWMs)){
      newOrder <- sample(1:ncol(PWMs[[m]]), ncol(PWMs[[m]]), replace=FALSE)
      mats[[m]] <- PWMs[[m]][, newOrder] 
      
      filtAux <- unlist(strsplit((filt[[m]]), split="[][]+"))
      filtAux <- filtAux[filtAux != ""]
      filtAux <- filtAux[newOrder] 
      filtAux <- paste("[", filtAux[1], "][", 
                       paste(filtAux[2:(length(filtAux)-1)], sep="][", 
                             collapse="][")
                       , "][", filtAux[length(filtAux)], "]", 
                       sep="", collapse="")
      newfilt[[m]] <- filtAux   
    }
    ranBlock[["PWMs"]] <- mats
    ranBlock[["filters"]] <- newfilt
    ranMatFilt[[paste("randommats1.", i, sep="", collapse="")]] <- ranBlock
  }  
  return(ranMatFilt)
}