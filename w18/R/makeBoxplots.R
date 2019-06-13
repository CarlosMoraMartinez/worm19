
#' makeBoxplots
#' 
#' Makes boxplots from GeneStats. Creates a new directory ./boxplots
#' and saves the boxplots generated in PDF format.
#' @param genestats stats generated from function getGeneStats
#' @param outnameaux name for the output PDF file
#' @return Does not return anything; writes file with finalstats directly.
#' @seealso getGeneStats
#' @keywords stats, boxplot
#' @import grDevices
#' @export

makeBoxplots <- function(genestats, outnameaux){
  if (!file.exists("./boxplots")){
    dir.create(file.path(getwd(), "./boxplots/"))
  }
  pdf(paste("./boxplots/", outnameaux, ".pdf", sep = "", collapse = ""))
  par(mfrow = c(1, 2))
  for(j in 2:(ncol(genestats) -1)){
    boxplot(genestats[, j]~genestats$group, col = c("red", "blue"))
    title(names(genestats)[j])
    differentgroups <- unique(genestats$group)
    pval <- wilcox.test(genestats[genestats$group == 
                                    differentgroups[1], j],
                        genestats[genestats$group == 
                                    differentgroups[2], j])$p.value
    pval <- round(pval, 4)
    legend("topright", 
           legend = paste("wilcox p-val=", pval, sep="", collapse=""))
  }
  dev.off()
}