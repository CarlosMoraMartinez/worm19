
#' joinSWFiles
#' 
#' When sequences have been split into many files, this function rejoins 
#' the superwindows files sorting them by condition.
#' 
#' @param swfiles the superwindows files to be rejoined
#' @param types different regular expression to get the conditions
#' @param rmark how different splits are marked
#' @param mark used to label the output file with all splits together (will 
#' replace rmark in the name)
#' @return A vector with the new swfiles.
#' @keywords join, superwindows
#' @export

joinSWFiles <- function(swfiles, types = list(maxlen = "[2785]000*_", 
                                              tfgroup = "EGhIS[pl]+"),
                        rmark = "split[0-9]+", mark  = "splitjoin"){
  dfsw <- data.frame(swfiles = swfiles, 
                     stringsAsFactors=FALSE)
  for(t in names(types)){
    dfsw[,t] <- unlist(regmatches(swfiles, 
                                  regexpr(swfiles, pattern = types[[t]])))
  }
  dfsep <- expand.grid(sapply(dfsw[, names(types)], unique))
  dfsep <- sapply(dfsep, as.character)
  
  ## Generate the new files
  newfiles <- NULL
  for(sep in 1:nrow(dfsep)){
    splits <- dfsw[apply(as.matrix(dfsw[, names(types)]), 1,
                         FUN = function(x, y) all(unlist(x)==y), 
                         unlist(dfsep[sep, ])), "swfiles"]
    splits <- as.character(unlist(splits))
    
    if(length(splits) > 0){
      fname <- gsub(rmark, mark, splits[1])
      newfiles[length(newfiles) + 1] <- fname
      alldata <- data.frame()
      for(f in splits){
        alldata <- rbind(alldata, read.table(f, sep = "\t", 
                                             header = TRUE, 
                                             stringsAsFactors = FALSE))
      }
      write.table(alldata, file = fname, sep = "\t", 
                  row.names = FALSE, quote = FALSE)
    }
  }
  return(newfiles)
}