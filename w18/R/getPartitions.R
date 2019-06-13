
#' getPartitions
#' 
#' Divides n tasks into k groups. Used to split datasets and run in parallel.
#' It returns a list in which each element is a vector with the indexes (from 1
#' to n) in each partition
#' 
#' @param n number of tasks
#' @param k number of groups
#' @keywords partitions, parallel
#' @import doParallel, foreach
#' @export

getPartitions <- function(n, k){
  tareas <- list()
  num <- floor(n/k)
  for(i in 1:k){
    if(i < k){
      tareas[[i]] <- c((1 + num * (i-1)):(num * i))
    }
    if(i == k){
      tareas[[i]] <- c((1 + num * (i-1)):n)
    }
  }
  return(tareas)
}
