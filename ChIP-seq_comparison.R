######################################################################
##        Analysis of various ChIP-seq experiments from ENCODE      ##
######################################################################
## Remove previous environment
rm(list = ls())

## Load packages and set working directory
require(GenomicFeatures)
require(rtracklayer)
require(Biostrings)
require(ChIPseeker)
require(clusterProfiler) 
require(ReactomePA)
require(org.Ce.eg.db)
require(TxDb.Celegans.UCSC.ce11.refGene)
require(ggplot2)

if (dir.exists("/media/erick/OS/Bioinformatics/TFM/DA/NGS")){
  parentDir <- "/media/erick/OS/Bioinformatics/TFM/DA/NGS"
} else { 
  parentDir <- getwd() 
}
setwd(parentDir)
if (!(dir.exists("plots"))){ system("mkdir plots") }

## We also need to create a TxDb object for ce11 genome:
txdb <- TxDb.Celegans.UCSC.ce11.refGene

# From now on, we will follow mostly the ChipSeeker pipeline:
# http://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html


################################
##   Load experimental data   ##
################################
## Create a list with all our ChIP-seq experiments:
list.files(pattern = ".bed", recursive = T)
files <- list(
  ceh43_emb = file.path(parentDir, "ceh-43/ceh-43_embryo.bed"),
  mef2_L1 = file.path(parentDir, "mef-2/mef-2_L1.bed"), 
  unc55_emb = file.path(parentDir, "unc-55/unc-55_embryo.bed"),
  unc55_L2 = file.path(parentDir, "unc-55/unc-55_L2.bed"),
  unc62_emb = file.path(parentDir, "unc-62/unc-62_embryo.bed"),
  unc62_L1 = file.path(parentDir, "unc-62/unc-62_L1.bed"),
  unc62_L2 = file.path(parentDir, "unc-62/unc-62_L2.bed"),              
  unc62_L3 = file.path(parentDir, "unc-62/unc-62_L3.bed"),              
  unc62_L4 = file.path(parentDir, "unc-62/unc-62_L4.bed"),             
  unc62_YA = file.path(parentDir, "unc-62/unc-62_YA.bed"),
  flh1_emb = file.path(parentDir, "flh-1/flh-1_embryo.bed"), 
  hlh13_emb = file.path(parentDir, "hlh-13/hlh-13_embryo.bed"),
  pat9_emb = file.path(parentDir, "pat-9/pat-9_embryo.bed"),
  sem2_emb = file.path(parentDir, "sem-2/sem-2_embryo.bed")
)


##############################################
##   Individual analysis of each BED file   ##
##############################################
for (file in 1:length(files)){
  ## Set name of the experiment and read file:
  peak <- as.character(files[file])
  name <- strsplit(peak, "/")[[1]][length(strsplit(peak, "/")[[1]])]
  name <- strsplit(name, "[.]")[[1]][1]
  dir <- strsplit(peak, "/")[[1]][length(strsplit(peak, "/")[[1]]) - 1]
  
  ## Load BED file:
  cat(sprintf("> Processing %s file...\n", name))
  
  extraCols_narrowPeak <- c(signalValue = "numeric", 
                            pValue = "numeric",
                            qValue = "numeric", 
                            peak = "integer")
  peak <- import(peak, format = "BED",
                 extraCols = extraCols_narrowPeak)
  
  ## Make coverage plot of all ChIP peaks:
  plot <- covplot(peak, weightCol = "signalValue")
  if(!(dir.exists(sprintf("%s/%s/plots/", parentDir, dir)))){
    dir.create(sprintf("%s/%s/plots/", parentDir, dir))
  }
  plot <- plot + ggtitle(sprintf("ChIP Peaks of %s over chromosomes", 
                                 paste0(strsplit(name, "_")[[1]][1], 
                                        " (", strsplit(name, "_")[[1]][2], ")" )))
  ggsave(plot, filename = sprintf("%s/%s/plots/%s_peaks_coverage_plot.jpg", 
                                  parentDir, dir, name))
  
  ## Profile of ChIP peaks binding to TSS regions:
  # First of all, for calculating the profile of ChIP peaks binding to 
  # TSS regions, we should prepare the TSS regions, which are defined as 
  # the flanking sequence of the TSS sites. Then, we align the peaks that
  # are mapping to these regions, and generate the tagMatrix.
  if (!exists("promoter")){
    promoter <- getPromoters(TxDb = txdb,
                             upstream = 5000,
                             downstream = 5000)
  }
  # tagMatrix <- getTagMatrix(peak, windows = promoter)
  # 
  # ## Average profile of ChIP binding to TSS regions
  # plotAvgProf <- plotAvgProf(tagMatrix, xlim = c(-5000, 5000),
  #                            xlab = "Genomic Region (5'->3')",
  #                            ylab = "Read Count Frequency")
  # plotAvgProf <- plotAvgProf + ggtitle(sprintf("Average profile of %s ChIP-seq binding to TSS regions", 
  #                                              paste0(strsplit(name, "_")[[1]][1], 
  #                                                     " (", strsplit(name, "_")[[1]][2], ")" )))
  # ggsave(plotAvgProf, filename = sprintf("%s/%s/plots/%s_average_profile_TSS_binding.jpg", 
  #                                        parentDir, dir, name))
  
  ## Make annotation of peaks:
  # Defaults to nearest gene!!
  # But we can add the info from flanking, more distant genes 
  peakAnno <- annotatePeak(peak, 
                           tssRegion = c(-5000, 5000),
                           TxDb = txdb,
                           addFlankGeneInfo = F)
  ## Visualize peak annotation:
  # Pie-like plot
  # plotAnnoPie(peakAnno) 
  # Bar-plot
  png(sprintf("%s/%s/plots/%s_plot_peakAnno_bar.png", 
              parentDir, dir, name))
  plotAnnoBar(peakAnno,
              title = sprintf("Feature distribution of %s peaks", 
                              paste0(strsplit(name, "_")[[1]][1], 
                                     " (", strsplit(name, "_")[[1]][2], ")" )))
  while (!is.null(dev.list())) dev.off()
  
  ## Visualize distribution of TF-binding loci relative to TSS:
  # The distance from the peak (binding site) to the TSS of the 
  # nearest gene is calculated by annotatePeak and reported in the 
  # output. plotDistToTSS calculate the percentage of binding sites 
  # upstream and downstream from the TSS of the nearest gene, and
  # visualize the distribution.
  plot <- plotDistToTSS(peakAnno)
  plot <- plot + ggtitle(sprintf("Distribution of TF-binding loci relative to TSS for %s peaks", 
                                 paste0(strsplit(name, "_")[[1]][1], 
                                        " (", strsplit(name, "_")[[1]][2], ")" )))
  ggsave(plot, filename = sprintf("%s/%s/plots/%s_distribution_TFBS_plot.jpg", 
                                  parentDir, dir, name))
  
  ## Functional enrichment analysis:
  pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId,
                            organism = "celegans")
  plot <- dotplot(pathway1, showCategory = 10)
  plot <- plot +  ggtitle(sprintf("Pathway enrichment analysis of %s peaks",
                                  paste0(strsplit(name, "_")[[1]][1],
                                         " (", strsplit(name, "_")[[1]][2], ")" )))
  ggsave(plot, filename = sprintf("%s/%s/plots/%s_pathway_all_genes.jpg",
                                  parentDir, dir, name),
         width = 12, height = 5, units = "in")
  # We could also do a barplot with the same information:
  # barplot(pathway1)
  
  
  ## Make a dataframe of peakAnno
  peakAnno.df <- as.data.frame(peakAnno)
  # Select those peaks which have a distance of max -5000 from TSS
  peaks.promoter <- peakAnno.df[which(-5000 <= peakAnno.df$distanceToTSS &
                                  peakAnno.df$distanceToTSS <= 5000), ]
  peaks.promoter
  write.table(peaks.promoter, file = sprintf("%s/%s/%s_promoter_peaks_5kb.bed",
                                             parentDir, dir, name), 
              sep = "\t", col.names = F, row.names = F, quote = F)
  # Convert gene IDs to ENTREZID
  keytypes(org.Ce.eg.db)
  promoter.genes <- bitr(peaks.promoter$geneId,
                         fromType = "ENTREZID", 
                         toType = "ALIAS",
                         OrgDb = "org.Ce.eg.db")
  write.table(promoter.genes$ALIAS, file = sprintf("%s/%s/%s_promoter_genes_5kb.txt",
                                             parentDir, dir, name), 
              sep = "\t", col.names = F, row.names = F, quote = F)
}


## We can also look into a region in more detail: for example, cat-2
filescat2 <- list(ceh43.emb = file.path(parentDir, "ceh-43/ceh-43_embryo.bed"),
                  mef2.L1 = file.path(parentDir, "mef-2/mef-2_L1.bed"),
                  unc55.emb = file.path(parentDir, "unc-55/unc-55_embryo.bed"),
                  unc62.emb = file.path(parentDir, "unc-62/unc-62_embryo.bed"))
# scales::show_col(hue_pal()(length(filescat2)))
p <- covplot(filescat2, weightCol = "V7", chr = "chrII",
             xlim=c(2.55e5, 2.57e5))
p <- p + scale_colour_manual(values = alpha(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
                                            0.5),
                             name = "")
p <- p + scale_fill_manual(values = alpha(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
                                          0.5),
                           name = "")
p <- p + ggtitle("ChIP Peaks in cat-2 CRR") +
  theme(axis.text.x = element_text(color = "#000000",
                                   size = 10),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(color = "#000000",
                                   size = 10),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5,
                                  size = 17),
        legend.title = element_text(face = "bold"),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.position = "top",
        panel.background = element_rect(fill = "white")) 
ggsave(p, file = "./plots/peaks_CRR_cat-2.jpg",
       width = 400, height = 200, units = 'mm')



###########################################################
##   Analysis of various ChIP-seq experiments together   ##
###########################################################
## Define promoter regions
if (!exists("promoter")){
  promoter <- getPromoters(TxDb = txdb,
                           upstream = 5000,
                           downstream = 5000)
}

## All experiments:
if (!(file.exists("all_chips_matrix.Rda"))){
  all.files <- files[c(1, 2, 3, 5, 11:14)]
  all.tagMatrix <- lapply(all.files, getTagMatrix, windows = promoter)
  save(all.tagMatrix, 
       file = "all_chips_matrix.Rda")
} else {
  load ("all_chips_matrix.Rda")
}
# Average profile binding
all.plotAvgProf <- plotAvgProf(all.tagMatrix, 
                               xlim = c(-5000, 5000),
                               xlab = "Genomic Region (5'->3')", 
                               ylab = "Read Count Frequency")
all.plotAvgProf <- all.plotAvgProf + 
  ggtitle("Average profile binding of all ChIP-seq experiments to TSS regions") +
  theme(axis.text.x = element_text(color = "#000000",
                                   size = 10),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(color = "#000000",
                                   size = 10),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5,
                                  size = 17),
        # legend.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.position = "top",
        panel.background = element_rect(fill = "white")) 
ggsave(all.plotAvgProf,
       filename = sprintf("%s/plots/average_profile_TSS_binding_all.jpg",
                          parentDir),
       width = 500, height = 300, units = 'mm')
# We can use plotAnnoBar to compare their genomic annotation:
all.files <- files[c(1, 2, 3, 5, 11:14)]
peakAnnoList <- lapply(all.files, annotatePeak, TxDb = txdb,
                       tssRegion = c(-5000, 5000),
                       verbose = FALSE)
plot <- plotAnnoBar(peakAnnoList) 
plot <- plot + ggtitle ("Feature distribution of embryonic and control experiments") +
  theme(axis.text.x = element_text(color = "#000000",
                                   size = 10),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(color = "#000000",
                                   size = 10),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5,
                                  size = 17),
        # legend.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.position = "right",
        panel.background = element_rect(fill = "white")) 
ggsave(plot,
       filename = sprintf("%s/plots/feature_distribution_all.jpg",
                          parentDir),
       width = 500, height = 300, units = 'mm')


## Control experiments:
if (!(file.exists("control_chips_matrix.Rda"))){
  control.files <- files[11:14]
  control.tagMatrix <- lapply(control.files, getTagMatrix, windows = promoter)
  save(control.tagMatrix, 
       file = "control_chips_matrix.Rda")
} else {
  load ("control_chips_matrix.Rda")
}
# Average profile binding
control.plotAvgProf <- plotAvgProf(control.tagMatrix, 
                                   xlim = c(-5000, 5000),
                                   xlab = "Genomic Region (5'->3')", 
                                   ylab = "Read Count Frequency")
control.plotAvgProf <- control.plotAvgProf + 
  ggtitle("Average profile binding of control ChIP-seq experiments to TSS regions") +
  theme(axis.text.x = element_text(color = "#000000",
                                   size = 10),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(color = "#000000",
                                   size = 10),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5,
                                  size = 17),
        # legend.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.position = "top",
        panel.background = element_rect(fill = "white"))
ggsave(control.plotAvgProf,
       filename = sprintf("%s/plots/average_profile_TSS_binding_controls.jpg",
                          parentDir),
       width = 500, height = 300, units = 'mm')

# Venn diagram of overlapping peaks
control.files <- files[11:14]
peakAnnoList <- lapply(control.files, annotatePeak, TxDb = txdb,
                       tssRegion = c(-5000, 5000),
                       verbose = FALSE)
plot <- plotAnnoBar(peakAnnoList) 
plot <- plot + ggtitle ("Feature distribution of control experiments") +
  theme(axis.text.x = element_text(color = "#000000",
                                   size = 10),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(color = "#000000",
                                   size = 10),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5,
                                  size = 17),
        # legend.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.position = "right",
        panel.background = element_rect(fill = "white")) 
ggsave(plot,
       filename = sprintf("%s/plots/feature_distribution_control.jpg",
                          parentDir),
       width = 500, height = 300, units = 'mm')
control.overlap <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(control.overlap)


## Emb/L1 experiments:
if (!(file.exists("embryonic_chips_matrix.Rda"))){
  embryonic.files <- files[c(1, 2, 3, 5)]
  embryonic.tagMatrix <- lapply(embryonic.files, getTagMatrix, windows = promoter)
  save(embryonic.tagMatrix, 
       file = "embryonic_chips_matrix.Rda")
} else {
  load ("embryonic_chips_matrix.Rda")
}
# Average profile binding
embryonic.plotAvgProf <- plotAvgProf(embryonic.tagMatrix, 
                                     xlim = c(-5000, 5000),
                                     xlab = "Genomic Region (5'->3')", 
                                     ylab = "Read Count Frequency")
embryonic.plotAvgProf <- embryonic.plotAvgProf + 
  ggtitle("Average profile binding of embryonic ChIP-seq experiments to TSS regions") +
  theme(axis.text.x = element_text(color = "#000000",
                                   size = 10),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(color = "#000000",
                                   size = 10),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5,
                                  size = 17),
        # legend.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.position = "top",
        panel.background = element_rect(fill = "white"))
ggsave(embryonic.plotAvgProf,
       filename = sprintf("%s/plots/average_profile_TSS_binding_embryonic.jpg",
                          parentDir),
       width = 500, height = 300, units = 'mm')

# Venn diagram of overlapping peaks
embryonic.files <- files[c(1, 2, 3, 5)]
peakAnnoList <- lapply(embryonic.files, annotatePeak, TxDb = txdb,
                       tssRegion = c(-5000, 5000),
                       verbose = FALSE)
plot <- plotAnnoBar(peakAnnoList) 
plot <- plot + ggtitle ("Feature distribution of embryonic experiments") +
  theme(axis.text.x = element_text(color = "#000000",
                                   size = 10),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(color = "#000000",
                                   size = 10),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5,
                                  size = 17),
        # legend.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.position = "right",
        panel.background = element_rect(fill = "white")) 
ggsave(plot,
       filename = sprintf("%s/plots/feature_distribution_embryonic.jpg",
                          parentDir),
       width = 500, height = 300, units = 'mm')
embryonic.overlap <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(embryonic.overlap)
