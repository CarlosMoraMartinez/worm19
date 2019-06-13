
#' getRandomPWMs
#' 
#' Given a set of PWMs, returns n lists of PWMs with shuffled columns.
#' 
#' @param PWMs list of matrixes, with names
#' @param n number of permutated lists (defaults to n=10).
#' @return A list of lists, with n permutations of each PWM (1 shuffled PWM
#' per originl PWM)
#' @keywords PWMs
#' @export

getRandomPWMs <- function(PWMs, n=10){
  mats <- PWMs
  ranPWMs <-list()
  for(i in 1:n){
    for(m in names(PWMs)) mats[[m]] <- PWMs[[m]][, sample(1:ncol(PWMs[[m]]), 
                                                          ncol(PWMs[[m]]), 
                                                          replace=FALSE)] 
    ranPWMs[[paste("randommats1.", i, sep="", collapse="")]] <- mats
  }  
  return(ranPWMs)
}