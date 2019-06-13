
#' makeFinalStats2
#' 
#' Makes stats from superwindows. Reads superwindows file and:
#' 1) calculates stats per gene
#' 2) calculates mean and sum per gene
#' 
#' Works better for genes which sequences are not divided in introns and
#' exons (it does not take into account dist_from_first_exon). Also, it
#' allows for absence of gene names.
#' 
#' @param swfiles files with superwindows (do not need to include group)
#' @param seqs dataframe with the original sequences analysed
#' @param dataset name of condition of the swfiles
#' @param makeplots if TRUE, makes boxplots of all variables (defaults to
#' FALSE)
#' @param typicalName whether the gene names are in the format 'gene_numseq' 
#' or not. If FALSE, "seqID" becomes "geneID" (defaults to FALSE).
#' 
#' @return A dataframe with the final stats.
#' @seealso getGeneStats, makeBoxplots
#' @keywords stats, boxplots, graphs
#' @export

makeFinalStats2 <- function(swfiles, seqs, dataset, 
                            makeplots = FALSE, typicalName = FALSE){
  # seqsReserva <- seqs
  finalstats <- data.frame()
  
  for(swfile in swfiles){
    allsuperwindows_full <- read.table(swfile, header = TRUE, 
                                       sep="\t", stringsAsFactors = FALSE)
    outdata <- cleanFilename(swfile)
    outnameaux <- outdata
    allsuperwindows <- allsuperwindows_full[allsuperwindows_full$seqID
                                            %in% seqs$seqID, ]
    
    if(nrow(allsuperwindows) != 0 & nrow(seqs) != 0) {
      ## Get table with data per gene and save it
      genestats <- getGeneStats(seqs, allsuperwindows, typicalName)
      write.table(genestats, file = paste("genestats", 
                                          outnameaux, sep="_", 
                                          collapse="_"),
                  row.names=FALSE, quote=FALSE, sep="\t")
      
      ## Get means by group and add them to the finalstats dataframe
      means <- sapply(genestats[, 2:(ncol(genestats)-1)], 
                      FUN = function(x, group){tapply(x, group, mean, 
                                                      na.rm = TRUE)}, 
                      group = genestats$group)
      sums <- sapply(genestats[, 2:(ncol(genestats)-1)], 
                     FUN = function(x, group){tapply(x, group, sum, 
                                                     na.rm = TRUE)}, 
                     group = genestats$group)
      
      if(length(unique(genestats$group)) > 1){
        colnames(sums) <- paste("sum", colnames(sums), sep="_")
        } else {
          names(sums) <- paste("sum", names(sums), sep="_")
          }
      
      if(length(unique(genestats$group)) == 1){
        finalstats <- rbind(finalstats, 
                            data.frame(dataset = outdata,
                                       group = unique(genestats$group),
                                       as.data.frame(t(means)),
                                       as.data.frame(t(sums)),
                                       stringsAsFactors = FALSE))
      } else { 
        finalstats <- rbind(finalstats, data.frame(dataset = outdata, 
                                                   group=rownames(means), 
                                                   means, sums, 
                                                   stringsAsFactors = FALSE))}
      
      if(any(names(seqs) == "GCs")){
        aux <- aggregate(GCs ~group.1, seqs, mean)
        names(aux)[1] <- "group"
        finalstats <- merge(finalstats, aux)
      }
      ## Plot the gene data
      if(makeplots & length(unique(genestats$group)) > 1){
        makeBoxplots(genestats, outnameaux)
      }
    }
  } # end read swfiles
  
  write.table(finalstats, 
              file = paste("finalstats_byIntervalsENCODE_", 
                           dataset, ".csv", sep = "", collapse = ""),
              row.names = FALSE, quote = FALSE, sep="\t", dec=",")
  
  return(finalstats)
}