
#' getGeneStats
#' 
#' This function gets basic stats for each gene.
#' 
#' @param seqs dataframe with original sequences data
#' @param allsuperwindows dataframe with superwindow data. Must have
#' "seqID", "width", "width_withoutreps"... columns (gene column no needed)
#' @param typicalName whether the seqIDs are in the form 'genename_seqnumber'
#' or not. If FALSE, the "seqID" column is used to split the data instead 
#' of the "gene" column (defaults to TRUE).
#' @return A dataframe with basic stats calculated for each gene.
#' @keywords gene, stats
#' @export

getGeneStats <- function(seqs, allsuperwindows, typicalName = TRUE){
  ## Create dataframe to store results
  if (typicalName) {
    genestats <- data.frame(gene = unique(seqs$gene),
                            stringsAsFactors = FALSE)
  } else {
    genestats <- data.frame(gene = unique(seqs$seqID), 
                            stringsAsFactors = FALSE)}
  
  ## Get superwindows gene names
  if (typicalName) {
    sw_genes <- sapply(allsuperwindows$seqID, 
                       FUN = function(x){a <- strsplit(x,
                                                       split="_")[[1]][1]}) 
  } else {
    sw_genes <- allsuperwindows$seqID
    }
  
  ## Get results for each row (1 row per gene)
  for (i in 1:nrow(genestats)){
    geneseqs <- seqs[seqs$gene == genestats[i, "gene"], ]
    
    # Get superwindows for a given gene
    sw <- allsuperwindows[which(sw_genes == genestats[i, "gene"]), ]
    
    # Separate superwindows in two categories:
    #   1) from original sequence 
    #   2) from masked sequence
    swmasked <- sw[sw$seq.mask == "masked", ]
    sworiginal <- sw[sw$seq.mask == "original", ]
    
    # Fill the results dataframe
    genestats[i, "num_seqs"] <- nrow(geneseqs)
    genestats[i, "totalLength"] <- sum(geneseqs$width)
    genestats[i, "unmaskedLength"] <- sum(geneseqs$width) - sum(geneseqs$masked_bp)
    genestats[i, "prop_masked"] <- sum(geneseqs$masked_bp)/sum(geneseqs$width)
    
    genestats[i, "num_superw"] <- nrow(sworiginal)
    genestats[i, "num_superw_Notmasked"] <- nrow(swmasked)
    genestats[i, "num_superw_perkb"] <- 1000*genestats[i, "num_superw"]/genestats[i, "totalLength"]
    genestats[i, "num_superw_masked_perkb"] <- 1000*genestats[i, "num_superw_Notmasked"]/ genestats[i, "unmaskedLength"]
    
    genestats[i, "coverage_superw"] <- sum(sworiginal$width)/sum(geneseqs$width)
    genestats[i, "coverage_superw_masked"] <- sum(swmasked$width_withoutreps)/genestats[i, "unmaskedLength"]
    genestats[i, "any_superwin"]<- as.numeric(genestats[i, "num_superw"] > 0)
    genestats[i, "any_superwin_masked"] <- as.numeric(genestats[i, "num_superw_Notmasked"] > 0)
    genestats[i, "group"] <- unique(geneseqs$group.1)
  }
  return (genestats)
}