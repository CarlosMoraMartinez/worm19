######################################################################
##        Analysis of a Single Cell dataset (Cao et al., 2017)      ##
######################################################################
## Remove previous environment 
rm(list = ls())

## Load packages and set working directory
require(monocle)
require(dplyr)
require(ggplot2)

if (dir.exists("/media/erick/OS/Bioinformatics/TFM/DA/SingleCell")){
  parentDir <- "/media/erick/OS/Bioinformatics/TFM/DA/SingleCell"
} else { 
  parentDir <- getwd() 
  }
setwd(parentDir) 


## Create directory for results and run analysis
if (!dir.exists("TFM_results")){
  dir.create("TFM_results")
  
  ## Load Single Cell dataset
  if (!file.exists("Cao_et_al_2017_vignette.RData")){
    download.file(
      "http://waterston.gs.washington.edu/sci_RNA_seq_gene_count_data/Cao_et_al_2017_vignette.RData",
      destfile = "Cao_et_al_2017_vignette.RData")
    load("Cao_et_al_2017_vignette.RData")
    save(list=ls(), file = "Cao_et_al_2017_vignette.RData")
  } else {
    load("Cao_et_al_2017_vignette.RData") 
    }
  setwd("TFM_results")  

  ## Exclude some unclassified cells that could interfere in the analysis
  dim(cds)
  cells <- which(pData(cds)$cell.type != c("Failed QC"))
  cds <- (cds)[,cells]
  cells <- which(pData(cds)$cell.type != c("Unclassified glia"))
  cds <- (cds)[,cells]
  cells <- which(pData(cds)$cell.type != c("Unclassified neurons"))
  cds <- (cds)[,cells]
  dim(cds)
  rm(cells)
  ## Exclude doublets as they could interfere in the analysis
  dim(cds.neurons)
  cells <- which(pData(cds.neurons)$neuron.type != "Doublets")
  cds.neurons <- (cds.neurons)[,cells]
  dim(cds.neurons)
  rm(cells)
  
  ## Check the markers for dopaminergic neurons and ciliated neurons
  plot.expr(cds, "dat-1", cell_size = 0.5)
  p <- plot.expr(cds.neurons, "dat-1", cell_size = 0.5)
  p <- p + 
    ggtitle("Expression of dat-1 in all the neuronal clusters") +
    theme(axis.text.x = element_text(color = "#000000",
                                     size = 10),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(color = "#000000",
                                     size = 10),
          axis.title.y = element_text(size = 13),
          plot.title = element_text(hjust = 0.5,
                                    size = 17),
          legend.position = "top",
          panel.background = element_rect(fill = "white")) 
  p
  filename <- paste0(getwd(), "/", "dat-1_neuronal_expression.png")
  ggsave(filename, plot = p,
         width = 250, height = 250, units = "mm")
  
  p <- plot.expr(cds, "che-11", cell_size = 0.5)
  p <- p + 
    ggtitle("Expression of che-11 in all the cells of the dataset") +
    theme(axis.text.x = element_text(color = "#000000",
                                     size = 10),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(color = "#000000",
                                     size = 10),
          axis.title.y = element_text(size = 13),
          plot.title = element_text(hjust = 0.5,
                                    size = 17),
          legend.position = "top",
          panel.background = element_rect(fill = "white")) 
  p
  filename <- paste0(getwd(), "/", "che-11_general_expression.png")
  ggsave(filename, plot = p,
         width = 250, height = 250, units = "mm")
  
  p <- plot.expr(cds.neurons, "che-11", cell_size = 0.5)
  p <- p + 
    ggtitle("Expression of che-11 in all the neuronal clusters") +
    theme(axis.text.x = element_text(color = "#000000",
                                     size = 10),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(color = "#000000",
                                     size = 10),
          axis.title.y = element_text(size = 13),
          plot.title = element_text(hjust = 0.5,
                                    size = 17),
          legend.position = "top",
          panel.background = element_rect(fill = "white")) 
  p
  filename <- paste0(getwd(), "/", "che-11_neuronal_expression.png")
  ggsave(filename, plot = p,
         width = 250, height = 250, units = "mm")
  
  # Alternatively we can ask the dataset
  cells <- cds
  plot.expr(cells, "dat-1")
  plot.expr(cells, "che-11")
  
  neurons <- cds.neurons
  plot.expr(neurons, "dat-1")
  plot.expr(neurons, "che-11")
  rm(cells, neurons)
  
  ######################################################################
  ########## Dopaminergic contrast with Cao et al. (2017) data #########
  ######################################################################
  ## DEG (Differential expression analysis) between dopa and ciliated
  if (file.exists("dopaminergic_genes_L2_Cao.rda")){
    load("dopaminergic_genes_L2_Cao.rda")
  } else {
    cat("Calculating DEG for DA neurons...\n")
    ciliated.vs.da.DEG = two.set.differential.gene.test(
      cds.neurons, is.cell.type(cds.neurons, "Ciliated sensory neurons"),
      is.cell.type(cds.neurons, "Dopaminergic neurons"), formal = T)  
    
    ## Filter results for q-value
    rows <- which(ciliated.vs.da.DEG$q.val <= 0.05)
    cil.vs.da.significant <- ciliated.vs.da.DEG[rows, ]
    
    ## Filter genes in dopaminergic neurons
    table(cil.vs.da.significant$higher.expr == "Set 2")
    da.rows <- which(cil.vs.da.significant$higher.expr == "Set 2")
    dopaminergic.genes <- cil.vs.da.significant[da.rows,]
    
    ## Add WormBase IDs to the dataframe
    for (i in 1:length(dopaminergic.genes[,1])){
      row <- which(fData(cds.neurons)[,2] == dopaminergic.genes[i, 1])
      if (i == 1){
        gene_id <- as.character(fData(cds.neurons)[row, 1])
      } else {
        gene_id <- c(gene_id, as.character(fData(cds.neurons)[row, 1]))
      }
    }
    if (length(gene_id) == nrow(dopaminergic.genes)){ # Check
    dopaminergic.genes <- tibble::add_column(dopaminergic.genes, 
                                             gene_id, .before = 1)
    }
    ## Change the name of a gene since is dead and merged into another
    row <- which(dopaminergic.genes$gene == "F36H2.5")
    dopaminergic.genes[row, "gene"] <- "F36H2.3"
    dopaminergic.genes[row, "gene_id"] <- "WBGene00009500"
    
    ## Write .tsv file and save filtered dataframe
    write.table(dopaminergic.genes, 
                file = "dopaminergic_genes_L2_Cao.tsv", 
                sep = "\t", row.names = F)
    save(dopaminergic.genes, file = 'dopaminergic_genes_L2_Cao.rda')
  } 
  
  
  ######################################################################
  ############ Neuron contrasts with Cao et al. (2017) data ############
  ######################################################################
  ## Compare every type of neuron against the others
  setwd(parentDir)

  ## Check if this part of the script has been already executed
  if (!dir.exists("neuronal_contrasts")){
    system("mkdir neuronal_contrasts")
    setwd("neuronal_contrasts")
    
    ## How many different types of neurons are there in the study?
    neuron.types
    
    ## Modified function from `is.neuron.type` function
    ## to include various types of neurons in the subsetting
    is.neuron.type2 <- function(cds, x) {
      with(pData(cds), !is.na(neuron.type) & neuron.type %in% x)
    }
    
    ## Make contrasts of every neuron type vs. all neurons
    for (i in neuron.types){
      # We take apart each neuron type from the rest of the neuronal
      # population in order to make the DEG contrasts: each contrast
      # will take into account a specific neuronal type (type)
      # vs. the rest of neuronal types of the dataset (all.neurons)
      type <- i
      neuron <- is.neuron.type(cds.neurons, type)
      all.neurons <- neuron.types != type
      all.neurons <- neuron.types[all.neurons]
      all.neurons <- is.neuron.type2(cds.neurons, all.neurons)
      
      cat(sprintf("Calculating DEG for %s neurons...\n", type))
      cores = detectCores() - 1
      type.vs.all = two.set.differential.gene.test(
        cds.neurons, neuron, all.neurons, formal = T, cores = cores)
      
      ## Filter results for q-value
      rows <- which(type.vs.all$q.val <= 0.05)
      type.vs.all <- type.vs.all[rows, ]
      
      ## Filter genes in type neurons
      # Select genes with higher expression in each type
      rows <- which(type.vs.all$higher.expr == "Set 1")
      type.genes <- type.vs.all[rows,]
      
      ## Add WormBase IDs to the dataframe
      for (j in 1:length(type.genes[,1])){
        row <- which(fData(cds.neurons)[,2] == type.genes[j, 1])
        if (j == 1){
          gene_id <- as.character(fData(cds.neurons)[row, 1])
        } else {
          gene_id <- c(gene_id, as.character(fData(cds.neurons)[row, 1]))
        }
      }
      if (length(gene_id) == nrow(type.genes)) { # Check
      type.genes <- tibble::add_column(type.genes, gene_id, .before = 1)
      }
      
      ## Save dataframe
      # Convert neuron type name if it has "/" caracters in it
      name <- strsplit(type, split = "/")
      for (j in 1:length(name[[1]])){
        if (j == 1){
          type <- paste(name[[1]][j])
        } else {
          type <- paste(type, name[[1]][j], sep=".")
        }
      }
      # Convert neuron type name to "xxxx.neurons" format
      filename <- strsplit(type, " ")
      filename <- paste0(filename[[1]][1:(length(filename[[1]]))], 
                         sep = ".", collapse = "")
      filename <- substr(filename, 1, nchar(filename) - 1)
      filename <- gsub('[-]','.', filename)
      filename <- gsub('[+]','.pos', filename)
      filename <- gsub('[()]','', filename)
      
      assign(paste(filename, "genes", sep = "."), type.genes)
      filename <- paste(filename, "genes", "rda", sep = ".")
      save(type.genes, file = filename)
      
      ## Save .tsv file
      filename <- gsub(".rda", ".tsv", filename)
      write.table(type.genes, file = filename,
                  row.names = F, sep = "\t")
    }
  } else {
    ## If the contrasts have already been made, load the .rda files
    setwd("neuronal_contrasts")
    files <- list.files()
    rda_files <- grep(list.files(), pattern = ".rda")
    for (i in rda_files){
      name <- strsplit(files[i], "[.]")
      name <- paste0(name[[1]][1:(length(name[[1]]) - 2)], 
                     sep = ".", collapse = "")
      name <- substr(name, 1, nchar(name) - 1)
      load(files[i])
      assign(name, type.genes)
    }
    rm(files, rda_files, i, name, type.genes)
  }
  
  
  ######################################################################
  ##########  Extract gene lists from some of these contrasts  #########
  ######################################################################
  ## We want to get the list of genes coming from:
  ## - ASE neuron (ASEL + ASER)
  ## - Touch receptor neurons 
  ## - RIA neuron
  ## - SDQ/ALN/PLN neurons
  ## - GABAergic neuron
  ## Plus the list of dopaminergic genes.
  
  ## Extract list of dopa genes
  setwd(paste0(parentDir, "/TFM_results/"))
  if (!dir.exists("gene_lists")){
    system("mkdir gene_lists")
  }
  setwd("./gene_lists")
  col.selected <- c("gene_id", "gene", "p.val", "q.val")
  dopa_genes <- dopaminergic.genes[col.selected]
  write.table(dopa_genes, file = "DA_genes.tsv",
              sep = "\t", row.names = F, quote = F)
  
  ## ASE neuron (ASEL + ASER)
  ASE <- rbind(ASEL, ASER)
  ASE_genes <- as_tibble(unique(ASE$gene))
  gene_id <- ASE$gene_id[which(ASE$gene %in% ASE_genes$value)]
  gene_id <- unique(gene_id)
  ASE_genes <- tibble::add_column(ASE_genes, gene_id, .before = 1)
  stats <- ASE[which(ASE_genes$gene_id %in% ASE$gene_id), c("p.val", "q.val")]
  ASE_genes <- cbind(ASE_genes, stats)
  colnames(ASE_genes) <- col.selected
  ASE_genes <- ASE_genes[order(ASE_genes$p.val), ]
  write.table(ASE_genes, file = "ASE_genes.tsv", 
              sep = "\t", row.names = F,  quote = F)
  
  ## Touch receptor neurons 
  touch_receptor_genes <- Touch.receptor[col.selected]
  write.table(touch_receptor_genes, file = "Touch_receptor_genes.tsv",
              sep = "\t", row.names = F, quote = F)
  
  ## RIA neuron
  RIA_genes <- RIA[col.selected]
  write.table(RIA_genes, file = "RIA_genes.tsv",
              sep = "\t", row.names = F, quote = F)
  
  ## SDQ/ALN/PLN neurons
  SDQ_ALN_PLN_genes <- SDQ.ALN.PLN[col.selected]
  write.table(SDQ_ALN_PLN_genes, file = "SDQ_ALN_PLN_genes.tsv",
              sep = "\t", row.names = F, quote = F)
  
  ## GABAergic neurons
  GABAergic_genes <- GABAergic[col.selected]
  write.table(GABAergic_genes, file = "GABAergic_genes.tsv",
              sep = "\t", row.names = F, quote = F)
} else {
  ## If the gene lists have already been made, we load the .tsv files
  setwd("TFM_results/gene_lists") 
  files <- list.files()
  for (i in files){
    name <- gsub("_genes.tsv", "", i)
    file <- read.table(i, header = T, sep = "\t")
    assign(name, file)
  }
  rm(i, file, files, name)
}
