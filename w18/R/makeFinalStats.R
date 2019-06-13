
#' makeFinalStats
#' 
#' Makes stats from superwindows. Reads superwindows file and:
#' 1) calculates stats per gene
#' 2) calculates mean and sum per gene
#' 
#' @param swfiles files with superwindows (do not need to include group)
#' @param seqs dataframe with the original sequences analysed
#' @param dataset name of condition of the swfiles
#' @param makeplots if TRUE, makes boxplots of all variables (defaults to
#' FALSE)
#' @param dist_intervals calculates stats for intervals centered around TSS
#' @param doDist if TRUE, only considers 1 dist interval (-Inf, Inf). Defaults
#' to TRUE.
#' 
#' @return A dataframe with the final stats.
#' @seealso filterSWbyDist, filterseqsbyDist, getGeneStats, makeBoxplots
#' @keywords stats, boxplots, graphs
#' @import stringr
#' @export

makeFinalStats <- function(swfiles, seqs, dataset, 
                           makeplots = FALSE, dist_intervals = NULL,
                           doDist = TRUE){
  
  if(is.null(dist_intervals) & doDist){
    dist_intervals <- expand.grid(c(-50001, -20000, -15000, -10000,
                                    -5000, -3000, -2000, -1000, -500, 0),
                                  c(10000000, 50001, 2000, 15000, 10000, 
                                    5000, 2000, 1000, 500, 0))
  } else if(!doDist){
    dist_intervals <- expand.grid(-Inf, Inf)
  }
  seqsReserva <- seqs
  finalstats <- data.frame()
  
  for(swfile in swfiles){
    allsuperwindows_full <- read.table(swfile, header = TRUE, 
                                       sep="\t", stringsAsFactors = FALSE)

    ## Select superwindows within given distances from first exon
    for(i in 1:nrow(dist_intervals)){
      dist <- as.numeric(dist_intervals[i, ])
      outdata <- paste(swfile, paste(dist, sep = ":",
                                    collapse = ":"), 
                      sep = "_", collapse = "_")
      outnameaux <- unlist(strsplit(outdata, "[._]"))
      outnameaux <- outnameaux[outnameaux != "csv" & outnameaux != ""]
      outnameaux <- paste(outnameaux, sep = "_", collapse = "_")
                    # paste(outnameaux[c(6:9, 11)], sep="_", collapse="_")
      
      allsuperwindows <- filterSWbyDist(allsuperwindows_full, dist)
      seqs <- filterseqsbyDist(seqsReserva, dist)
      
      if(nrow(allsuperwindows) != 0 & nrow(seqs) !=0){
        ## Get table with data per gene and save it
        genestats <- getGeneStats(seqs, allsuperwindows)
        write.table(genestats, file = paste("genestats", outnameaux, 
                                            sep = "_", collapse = "_"),
                    row.names = FALSE, quote = FALSE, sep="\t")
        
        ## Get means by group and add them to the finalstats dataframe
        means <- sapply(genestats[, 2:(ncol(genestats)-1)], 
                        FUN = function(x, group){tapply(x, group, mean, 
                                                        na.rm = TRUE)}, 
                        group = genestats$group)
        sums<- sapply(genestats[, 2:(ncol(genestats)-1)], 
                      FUN = function(x, group){tapply(x, group, sum,
                                                      na.rm = TRUE)}, 
                      group = genestats$group)
        
        names(sums) <- paste("sum", names(sums), sep="_")
        finalstats <- rbind(finalstats, data.frame(dataset = outdata,
                                                   group = rownames(means), 
                                                   means, sums, 
                                                   stringsAsFactors=FALSE))
        
        ## Plot the gene data
        if(makeplots) makeBoxplots(genestats, outnameaux)
        
      }
    } # end different dists
  } # end read swfiles
  write.table(finalstats, file = paste("finalstats_byIntervals_", dataset,
                                       ".csv", sep = "", collapse = ""),
              row.names = FALSE, quote = FALSE, sep = "\t", dec = ",")
  
  return(finalstats)
}
