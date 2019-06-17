################################################################
##                   MAIN WORM19 PIPELINE                     ##
################################################################
######################################################
##         1. Search all motifs and windows         ##
######################################################
rm(list = ls()) ## Clean workspace

### This script has two parts:
### 1) Get data
parentdir <- "/media/erick/OS/Bioinformatics/TFM/DA/worm19"
setwd(parentdir)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(w18, plyr, Biostrings, BSgenome, GenomicRanges, 
               doParallel, foreach, stringr, data.table)

## Paths
seqsfolder <- "../worm19/all_data/all_seqs" # Path to the sequences
newdir <- "../worm19/190522"              # Path to the output
if (!dir.exists(newdir)) {
  dir.create(newdir)
} else {
  print ("El directorio ya existe.")
}
superwindir <- "../worm19/190522/sw1/"     # Path for superwindows
if (!dir.exists(superwindir)){
  dir.create(superwindir)
} else {
  print ("El directorio sw1 ya existe.")
}

## Settings
PWMset <- "worm19"
sequenceSet <- "worm19"  # Organism 
minscore <- "70%"   # Min. score
PWMsize <- "small"  # PWM size
## Calculate windows
window_min_motifs <- c(6, 7, 8)
max_widthVec <- c(700)

## PWM names
names.mymotifs <- c("COUP", "DLX1", "DLX2", "ETS",
                    "MADS", "MEIS", "PAX", "PBX")
mymotifs <- vector("list", length(names.mymotifs))
names(mymotifs) <- names.mymotifs
nmotifs <- length(mymotifs) # 8 motifs

## PWM files
for (i in 1:length(names.mymotifs)){
  name <- paste0(names.mymotifs[i], "_PWM.txt")
  if (i == 1){
    names.PWMfiles <- name
  } else {
    names.PWMfiles <- c(names.PWMfiles, name)
  }
}

## Get PWMs if they are read from a file
for(i in 1:length(names.mymotifs)){
  temp <- read.table(file = file.path(paste(c(getwd(), 
                                              "/original_PWMs"), 
                                            sep = "", collapse = ""), 
                                      names.PWMfiles[i]), 
                     sep = "\t", header = FALSE)
  colnames(temp) <- c("A", "C", "G", "T")
  mymotifs[[names.mymotifs[i]]] <- t(as.matrix(temp))
}
allmymotifs <- mymotifs

## This is the filter for worm matrixes
## It has to be the same lenght as names.PWMfiles
## The RE has been provided by experimental data (Jimeno, Flames)
myfilt <- c("TGACC", # COUP
            "TAATT", # DLX1
            "[GAT]ATAATT[GT]", # DLX2
            "[GAC][AC]GGA[AT][GA]", # ETS
            "[AT][AT][AT][TAG]TAG", # MADS
            "[GAT]TGTC[GAT]", # MEIS (DTGTCD = HGACAH)
            "GGA[GA][GC]A", # PAX
            "GAT[ACGT][ACGT]AT") # PBX
names(myfilt) <- names(mymotifs)

## Get sequence files
seqfiles <- list.files(path = seqsfolder, full.names = TRUE)

### 2) Once settings are fixed and data is loaded,
### it looks for windows and superwindows in each sequence
### in which the motifs appear to look for enhancer regions
## Start the cluster
packs <- search()
packs <- packs[grep(packs, pattern="package:")]
packs <- sapply(packs, 
                FUN = function(s){unlist(strsplit(s, split=":"))[2]})

nCPUcores = detectCores()
if (nCPUcores < 3) {
  registerDoSEQ()
} else {
  cl = makeCluster(nCPUcores - 1, 
                   type = "FORK", outfile = "log_w19_main.txt")
  registerDoParallel(cl)
}

## For each one of the sequence files
foreach(i = 1:length(seqfiles), .packages = packs, 
        .export=ls(envir=globalenv()),
        .verbose=TRUE) %dopar% {
          ## Read sequences
          file <- seqfiles[i]
          seqs <- read.table(file, header=TRUE, sep="\t", 
                             stringsAsFactors=FALSE) 
          
          # Each chain (seqs$original.seq) has 10k bases
          # seqs<-seqs[!seqs$seqID %in% wrongseq, ]  
          seqs <- seqs[!duplicated(seqs[, names(seqs)!="seqID"]), ]
          seqs$masked.seq <-'A'
          seqs$group <- 'not_grouped'
          
          seqfile <- unlist(strsplit(file, split="/"))
          seqfile <- trimextcsv(paste(newdir, '/', 
                                      seqfile[length(seqfile)], 
                                      sep="", collapse=""))
          
          # Search motifs in each sequence in variable "seqs" and store such
          # dataframes in a list. The list can be subsetted with seqs$ID. 
          # Each element of the list (each seqs$ID) is another list 
          # with 2 data frames: [1] original sequence; and [2] masked sequence
          swfiles <- NULL
          
          ## Search motifs and make windows and superwindows
          for(PWMsubset in (list(names.mymotifs))){
            cat(paste(Sys.info()[['nodename']], 
                      Sys.getpid(), sep='-'), "now searching motifs\n")
            
            mymotifs <- allmymotifs[PWMsubset]
            PWMsubsetname <- paste(unlist(sapply(PWMsubset, substr, 1, 1)), 
                                   sep=",", collapse="")
            motifList <- globalSearchMotifs(seqs, PWMsubset, mymotifs, 
                                            myfilt, seqfile)
            
            ## Loop for different window max widths (currently 700 and 800) 
            for(max_width in max_widthVec){
              cat(paste(Sys.info()[['nodename']], 
                        Sys.getpid(), sep='-'), "now making windows of", 
                  max_width, "\n")
              windows <- globalMakeWindows(seqs, motifList, seqfile, PWMset, 
                                           PWMsubsetname, max_width, 
                                           min(window_min_motifs))
              for(min_motifs in window_min_motifs){
                cat(paste(Sys.info()[['nodename']], 
                          Sys.getpid(), sep='-'), "now making superwindows of",
                    min_motifs, "\n")
                w_aux <- windows[windows$different_motifs == min_motifs, ]
                swfiles[length(swfiles)
                        + 1] <- makeSuperwindows(filename = paste(seqfile,
                                                                  PWMset, 
                                                                  PWMsubsetname, 
                                                                  max_width, 
                                                                  min_motifs,
                                                                  'mots',
                                                                  sep="_", 
                                                                  collapse="_"), 
                                                 seqs = seqs,
                                                 outfolder = superwindir, 
                                                 win = w_aux)
              } ## for window_min_motifs
            } ## for different window widths (max_widthVec)
          } ## for all the motifs (mymotifs)
        } ## for parallel execution (each seqfile)

## Don't forget to stop the cluster when finished
stopCluster(cl)





######################################################
##         2. Build complete gene summary           ##
######################################################
####################################################################
##              Load template of complete_gene_summary            ##
####################################################################
rm(list = ls())

parentDir <- "/media/erick/OS/Bioinformatics/TFM/DA/worm19"
currentDir <- paste0(parentDir, "/random_genes")
setwd(currentDir)

## Read original all genes file
all <- read.table("template_worm19_complete_gene_summary.tsv", 
                  sep = "\t", header = T, stringsAsFactors = F)

## Select only the columns that we're gonna use
all <- all[1:11]


######################################################################
##       Check how many windows of n motifs has each species        ##
######################################################################
## Load homology table
currentDir <- paste0(parentDir, "/reference")
setwd(currentDir)
homology <- read.table("worm_homologues_of_celegans.tsv",
                       sep = "\t", header = T, stringsAsFactors = F)
colnames(homology)[1] <- "elegans_genes"

## Define iterators: species and number of motifs
n.motifs <- c(6, 7, 8)
species <- c("brenneri", "briggsae", "elegans", "japonica", "remanei")
motifs <- c("COUP", "DLX1", "DLX2", "ETS", "MADS", "MEIS", "PAX", "PBX")

## Prepare final dataframe
for (h in 1:length(n.motifs)){
  for (i in 1:length(species)){
    col <- ncol(all) + 1 ; all[,col] <- 0
    colnames(all)[col] <- sprintf("any_%d_%s", n.motifs[h], species[i])
    for (j in 1:length(motifs)){
      if (n.motifs[h] == n.motifs[length(n.motifs)]) { next }
      col <- ncol(all) + 1 ; all[,col] <- 0
      colnames(all)[col] <- sprintf("any_%d_without_%s_%s", 
                                    n.motifs[h], motifs[j], species[i])
    }
  }
}
for (h in 1:length(n.motifs)) {
  col <- ncol(all) + 1 ; all[,col] <- 0
  colnames(all)[col] <- sprintf("ortho_with_win_%d_mots", n.motifs[h])
}


## Complete final dataframe with data
for (i in 1:length(species)){
  ## Load reference genes table for each species
  currentDir <- paste0(parentDir, "/reference")
  setwd(currentDir)
  filename <- paste("gene", "reference", "c", species[i], "10k.tsv",
                    sep = "_")
  reference <- read.table(file = filename, sep="\t", header=T, 
                          stringsAsFactors=F)
  
  ## Load windows files:
  currentDir <- paste0(parentDir, "/190522")
  setwd(currentDir)
  files <- list.files()
  # Select text files of 700 pb without seqs
  files.select <- files[grep(species[i], files)] 
  files.select <- files.select[grep(".csv", files.select)]
  files.select <- files.select[grep("700", files.select)]
  file.win <- tryCatch( { 
    file.win <- read.table(files.select,
                           sep = "\t", header = T, stringsAsFactors = F)
  }, error = function(e) {
    print (sprintf("Error reading %s windows file.", species[i]))
    file.win <- NA
  })
  if (is.na(file.win)){ next }
  
  
  ## For each one of the sequence files:
  for(j in 1:length(n.motifs)){
    # Split file.win for windows of n.motifs[j] motifs or more.
    file <- file.win[which(file.win$different_motifs >= n.motifs[j]), ]
    # For all the possible genes of that species (all the genes in the
    # windows file), count how many windows have been detected and
    # write it in the `all` dataframe
    genes <- unique(file$seqID)
    genes.unique.ids <- character()
    for (k in 1:length(genes)){
      id <- strsplit(genes[k], "_")[[1]][1]
      genes.unique.ids <- c(genes.unique.ids, id)
    }
    genes.unique.ids <- unique(genes.unique.ids)
    for (k in 1:length(genes)){
      seqid <- genes[k]
      id <- strsplit(seqid, "_")[[1]][1]
      # Is this the first time we access to this gene?
      if (!(id %in% genes.unique.ids)) { next }
      wbid <- reference$gene[which(reference$seqID == seqid)]
      seqids <- reference$seqID[which(reference$gene == wbid)]
      
      # To accelerate the searches we use the homology table we loaded
      column <- paste("any", n.motifs[j],  species[i],  sep = "_")
      ref.col <- paste(species[i], "genes", sep = "_")
      row <- which(wbid == homology[ref.col])
      # If we have found the gene...
      if (length(row) > 0){
        # ...we search for the C. elegans orthologue
        w.id <- homology[row[1], "elegans_genes"]
        row <- which(all$elegans == w.id)
      }
      # If we haven't found a match in the homology table,
      # we search the gene in the `all` dataframe
      if (length(row) == 0) {
        count = 0
        while (length(row) == 0){
          homologue = F
          worm.ids <- strsplit(x = (all[, species[i]]), split = "_")
          for (l in 1:length(worm.ids)){
            row <- which(wbid %in% worm.ids[[l]])
            count = count + 1
            if (length(row) != 0){ 
              row <- l
              break 
            }
          } # end for (l)
          if (count == length(worm.ids)) { 
            row <- NA
            break 
          }
        }  # end while 
      } # end else
      
      # How many matches on the windows file have that specific gene?
      matches = 0
      for (m in 1:length(seqids)){
        match <- seqids[m] %in% file$seqID
        if (match){
          windows <- sum(file$seqID == seqids[m])
          matches <- sum (matches, windows)
        }
      } # end for (m)
      
      # If we have found where to put the result (row, column),
      # assign the number of matches to that position
      if (!(is.na(row))){
        if (matches == 0){
          # That gene has zero matches on the windows file
          all[row, column] <- 0
        } else {
          all[row, column] <- all[row,column] + matches
        }
        # Now that we have assign the total of windows for that gene, 
        # we analyse the composition of the windows:
        # We want to check how many windows there are for each gene
        # without each one of the motifs.
        for (n in 1:length(motifs)){
          column <- sprintf("any_%d_without_%s_%s", n.motifs[j], 
                            motifs[n], species[i])
          windows <- file[which(file$seqID %in% seqids),]
          for(win in 1:nrow(windows)){
            signature.split <- strsplit(windows$signature[win], "_")
            # If that motif is not on the signature, we sum +1
            if (!(motifs[n] %in% signature.split[[1]])) { 
              all[row, column] <- all[row,column] + 1
            }
          } # end win
        } # end for length(motifs) (n)
      } # end !(is.na(row))
      
      # We eliminate this gene from the list of unique ids 
      # since we have already counted all the windows for 
      # the promoter + intron sequences.
      if (id %in% genes.unique.ids){
        genes.unique.ids <- genes.unique.ids[which(genes.unique.ids != id)]
      }
    } # end for length(genes) (k)
    
    ## Informative of the state of loop
    cat("> Finalizado C.", species[i], "con", n.motifs[j], "motivos.\n")
    
  } # end for length(n.motifs) (j)
} # end for length(species) (i)


######################################################################
##      Check the number of orthologs for each C. elegans gene      ##
######################################################################
for (j in 1:length(n.motifs)){
  cols <- colnames(all)[grep(n.motifs[j], colnames(all))]
  cols <- cols[grep("without", cols, invert = T)]
  target <- which(colnames(all) == cols[length(cols)])
  # Without elegans: we search orthologs
  cols <- cols[grep("elegans", cols, invert = T)]
  cols <- cols[1:(length(cols) - 1)]
  
  # For each row of `all`, we check all cols with j number of motifs...
  for (row in 1:nrow(all)){
    count = 0
    for (m in 1:length(cols)){
      if (all[row, cols[m]] != 0){
        count = count + 1
      } else { next }
    }
    # ...and assign the number of orthologs
    all[row, target] <- count
  }
}


######################################################################
##                     Write final output file                      ##
######################################################################
currentDir <- paste0(parentDir, "/random_genes")
setwd(currentDir)
write.table(all, 
            file = "190522_worm19_complete_gene_summary_windows.tsv",
            sep = "\t", row.names = F, quote = F)





######################################################
##         3. Get random sets of all genes          ##
######################################################
###################################################################
##   Create random gene sets for all types of the gene.summary   ##
###################################################################
## Clean environment and set working directory
rm(list=ls())
setwd("/media/erick/OS/Bioinformatics/TFM/DA/worm19/random_genes")

## Load libraries
pacman::p_load(parallel, doParallel, data.table, magrittr)

## Load file
f <- list.files()
file <- f[grep("190522_worm19_complete_gene_summary_windows.tsv", f)]
all <- read.table(file, sep="\t", stringsAsFactors=F, header=T)

## Create dir for saving the results
if (!(dir.exists("all_types_windows"))){
  system("mkdir all_types_windows")
}
setwd("./all_types_windows")
if (!(dir.exists("all_random"))) { system("mkdir all_random") }

## Split data according to "type" category 
types <- sort(unique(all$type[all$type != "other" & all$type != "neuronal"]))


for (i in 1:length(types)){
  ## Select type for each iteration
  type <- all[all$type == types[i], ]
  alloth <- all[all$type != types[i],]
  
  ## Select those genes that have an upstream distance greater or equal as
  ## the minimum upstream distance in dopa genes (or the type of election) 
  ## But never less!
  othtrim <- alloth[alloth$upstream_considered >= min(type$upstream_considered),]
  
  ## Calculate Empirical Cumulative Distribution Function for type genes
  percentile <- ecdf(type$upstream_considered + type$intronic_considered)
  plot(percentile)
  summary(percentile)
  
  ## In what part of the distribution are the upstreams of the genes selected?
  prob <- percentile(othtrim$upstream_considered + othtrim$intronic_considered) 
  
  ## Settings for the random groups
  groups <- list()
  needed <- 10000
  p <- 0.05
  all_othtrim <- othtrim
  all_type <- type
  all_groups <- data.frame()
  
  ## Settings for paralelization
  ncores = detectCores()
  cl <- makeCluster(ncores - 1, type="FORK")
  registerDoParallel(cl)
  
  ## For every conservation score (from 0 to 4), we get 10000 rows with 
  ## genes that match criteria (as many genes as type genes we have)
  for (orts in 0:4){
    # Select all type genes with that number of orthologs
    type <- all_type[all_type$num_orthologs >= orts, ]
    # Select every other gene with that number of orthologs
    othtrim <- all_othtrim[all_othtrim$num_orthologs >= orts, ]
    gsize <- nrow(type)
    sample <- 1:nrow(othtrim)
    
    ## We have to decide which genes we want to select:
    this_orts <- parLapply(cl=cl, X=1:needed, fun=function(x, type,
                                                           othtrim, p, 
                                                           gsize, sample,
                                                           orts){
      a <- 1
      b <- 1
      percentile <- ecdf(type$upstream_considered + 
                           type$intronic_considered)
      prob <- percentile(othtrim$upstream_considered + 
                           othtrim$intronic_considered)
      while(TRUE){
        # We select a sample of `gsize` size with the probability given by prob
        # and we save the indexes of those genes in the variable `ind`
        ind <- sample(sample, gsize, replace=FALSE, prob = prob)
        
        # We make a Wilcox Test between the upstream/intronic region of the type genes and
        # all the other genes selected and look for a size result that is not significantly 
        # different (that is, that has a p-value greater than the one we selected (0.05)).
        # If both (the upstream and the intronic region) are similar to the type size,
        # we exit the loop.
        a <- wilcox.test(othtrim$upstream_considered[ind], 
                         type$upstream_considered)$p.value
        b <- wilcox.test(othtrim$intronic_considered[ind], 
                         type$intronic_considered)$p.value
        if((a >= p) 
           & (b >= p)
        ){ break }
      } # end while
      ## Once we have selected the `gsize` genes, we concatenate 
      ## their WormBaseID and save it as a single string in the variable `res`
      aux <- othtrim[ind, ]
      res <- data.frame(genes = paste(aux$elegans, sep="_", collapse="_"))
      
      ## For the number of motifs considered (6, 7, 8), we select the column in C. elegans 
      ## that corresponds to the window with that number of motifs (`col`) and count
      ## the rows that have more than 0 win for that gene (`var`). 
      ## Finally, we set the percentage of (those) genes that have windows
      ## in C. elegans for that number of motifs.
      mots <- paste("any", 8:6, sep="_")
      for(m in mots){
        col <- names(aux)[grepl("_elegans", names(aux)) & 
                            grepl(m, names(aux))]
        col <- col [grep ("without", col, invert = T) ]
        var <- aux[, col] > 0
        res[1, m] <- sum(var)/nrow(aux)
      }
      return(res)
    }, type=type, othtrim=othtrim, p=p, gsize=gsize, sample=sample, orts=orts) %>% rbindlist %>% as.data.frame
    this_orts[, "conservation"] <- orts
    all_groups <- rbind(all_groups, this_orts)
  } #for othologs
  
  ## Stop the cluster
  stopCluster(cl)
  
  ## Save results
  name <- paste0("./all_random/", types[i], 
                 "_all_random_groups_WithoutOverlaps.tsv")
  write.table(all_groups, file = name,
              sep="\t", row.names = F, quote = F)
}





######################################################
##       4. Get random sets of conserved genes      ##
######################################################
###################################################################
##   Create random gene sets for all types of the gene.summary   ##
##        according to the number of orthologous with win        ##
###################################################################
## Clean environment and set working directory
rm(list=ls())
setwd("/media/erick/OS/Bioinformatics/TFM/DA/worm19/random_genes")

## Load libraries
pacman::p_load(parallel, doParallel, data.table, magrittr)

## Load file
f <- list.files()
file <- f[grep("190522_worm19_complete_gene_summary_windows.tsv", f)]
all <- read.table(file, sep="\t", stringsAsFactors=F, header=T)

## Create dir for saving the results
if (!(dir.exists("all_types_windows"))){
  system("mkdir all_types_windows")
}
setwd("./all_types_windows")
if (!(dir.exists("conserved_random"))) { system("mkdir conserved_random") }

## Split data according to "type" category 
types <- unique(all$type[all$type != "other" & all$type != "neuronal"]) 


for (i in 1:length(types)){
  ## Select type for each iteration
  type <- all[all$type == types[i], ]
  alloth <- all[all$type != types[i],]
  
  ## Select those genes that have an upstream distance greater or equal as
  ## the minimum upstream distance in dopa genes (or the type of election) 
  ## But never less!
  othtrim <- alloth[alloth$upstream_considered >= 
                      min(type$upstream_considered),]
  
  ## Calculate Empirical Cumulative Distribution Function for type genes
  percentile <- ecdf(type$upstream_considered + type$intronic_considered)
  plot(percentile)
  summary(percentile)
  
  ## In what part of the distribution are the upstreams of the genes selected?
  prob <- percentile(othtrim$upstream_considered + 
                       othtrim$intronic_considered) 
  
  ## Settings for the random groups
  groups <- list()
  needed <- 10000
  p <- 0.05
  all_othtrim <- othtrim
  all_type <- type
  all_groups <- data.frame()
  
  ## Settings for paralelization
  ncores = detectCores()
  cl <- makeCluster(ncores - 1, type="FORK")
  registerDoParallel(cl)
  
  ## For every conservation score (from 0 to 4), we get 10000 rows with 
  ## genes that match criteria (as many genes as type genes we have)
  for (orts in 0:4){
    # Select all type genes with that number of orthologs
    type <- all_type[all_type$num_orthologs >= orts, ]
    # Select every other gene with that number of orthologs
    othtrim <- all_othtrim[all_othtrim$num_orthologs >= orts, ]
    gsize <- nrow(type)
    sample <- 1:nrow(othtrim)
    
    ## We have to decide which genes we want to select:
    this_orts <- parLapply(cl=cl, X=1:needed, fun=function(x, type, 
                                                           othtrim, p, 
                                                           gsize, sample,
                                                           orts){
      a <- 1
      b <- 1
      percentile <- ecdf(type$upstream_considered + 
                           type$intronic_considered)
      prob <- percentile(othtrim$upstream_considered + 
                           othtrim$intronic_considered)
      while(TRUE){
        # We select a sample of `gsize` size with the probability given by prob (less random)
        # and we save the indexes of those genes in the variable `ind`
        ind <- sample(sample, gsize, replace=FALSE, prob = prob)
        
        # We make a Wilcox Test between the upstream/intronic region of the type genes and
        # all the other genes selected and look for a size result that is not significantly 
        # different (that is, that has a p-value greater than the one we selected (0.05)).
        # If both (the upstream and the intronic region) are similar to the type size,
        # we exit the loop.
        a <- wilcox.test(othtrim$upstream_considered[ind], 
                         type$upstream_considered)$p.value
        b <- wilcox.test(othtrim$intronic_considered[ind], 
                         type$intronic_considered)$p.value
        if((a >= p) 
           & (b >= p)
        ){ break }
      } # end while
      ## Once we have selected the `gsize` genes, we concatenate their WormBaseID
      ## and save it as a single string in the variable `res`
      aux <- othtrim[ind, ]
      res <- data.frame(genes = paste(aux$elegans, sep="_", collapse="_"))
      
      ## For the number of motifs considered (6, 7, 8), we select the column in C. elegans 
      ## that corresponds to the window with that number of motifs (`col`) and count
      ## the rows that have more than 0 win for that gene (`var`). 
      ## Finally, we set the percentage of (those) genes that have windows
      ## in C. elegans for that number of motifs.
      mots <- paste(8:6, "mots", sep="_")
      rows.conserved <-  which(aux$num_orthologs == orts)
      for(m in mots){
        rows.conserved.win <- which(aux[, sprintf("ortho_with_win_%s", m)]
                                    == orts)
        rows <- rows.conserved[which(rows.conserved %in% rows.conserved.win)]
        col <- which(colnames(type) == sprintf("ortho_with_win_%s", m))
        var <- type[rows, col] > 0
        res[1, m] <- sum(var)/length(rows.conserved)
      }
      return(res)
    }, type=type, othtrim=othtrim, p=p, gsize=gsize, sample=sample, orts=orts) %>% rbindlist %>% as.data.frame
    this_orts[, "conservation"] <- orts
    all_groups <- rbind(all_groups, this_orts)
  } #for othologs
  
  ## Stop the cluster
  stopCluster(cl)
  
  ## Save results
  name <- paste0("./conserved_random/", types[i],
                 "_all_random_groups_WithoutOverlaps.tsv")
  write.table(all_groups, file = name,
              sep="\t", row.names = F, quote = F)
}





######################################################
##          5. Extract info + write tables          ##
######################################################
#######################################################
##     Extract the info of the w19 windows files     ##
#######################################################
## Clear environment
rm(list = ls())

## Load libraries
require(Biostrings)
require(GenomicFeatures)
require(ggplot2)
require(ggseqlogo)

## Set directories
parentDir <- "/media/erick/OS/Bioinformatics/TFM/DA/worm19"
referenceDir <- paste0(parentDir, "/reference")
winDir <- paste0(parentDir, "/190522")
seqsDir <- paste0(parentDir, "/all_data/all_seqs")

## Load complete gene summary file
setwd(paste0(parentDir, "/random_genes"))
f <- list.files()
file <- f[grep("190522_worm19_complete_gene_summary_windows.tsv", f)]
gene.summary <- read.table(file, sep = "\t", 
                           stringsAsFactors = F, header = T)

## Initialize some variables
species <- colnames(gene.summary)[1:5]
types <- sort(unique(gene.summary$type[gene.summary$type != "other" & 
                                         gene.summary$type != "neuronal"])) 
n.motifs <- c(6, 7, 8)
motifs <- c("COUP", "DLX1", "DLX2", "ETS", 
            "MADS", "MEIS", "PAX", "PBX")

## Create output directory and it as setwd
newDir <- paste0(parentDir, "/win_results")
if (!(dir.exists(newDir))) {
  dir.create(newDir)
} else { print ("El directorio ya existe!")}
setwd(newDir)


## Run analysis
if (!(file.exists("DA_genes.Rda"))){
  table.types <- data.frame()
  table.motifs <- data.frame()
  ## Make table.types, table.motifs and table.genes
  for (i in 1:length(types)) { 
    for (h in 1:length(species)) { 
      for (n in 1:length(n.motifs)) { 
        ## Get the index from where we start adding info to the table at each iteration
        if (h == 1 & n == 1) { index <- 2 }
        
        ## Select type genes
        type <- gene.summary[gene.summary$type == types[i],]
        assign(types[i], type)
        
        ##### 1) Make table with information of windows for each neuronal type #####
        # Annotate number of genes of each type
        table.types[1, i] <- types[i]
        table.types[2, i] <- nrow(type)
        if (length(rownames(table.types)) == 2) {
          rownames(table.types) <- c("type", "n.genes")
        }
        ## Calculate number of total windows
        index <- index + 1
        table.types[index, i] <- sum(type[, sprintf("any_%d_%s",
                                                    n.motifs[n], species[h])])
        rownames(table.types)[index] <- sprintf("any_%d_%s",
                                                n.motifs[n], species[h])
        # Calculate percentage of genes with windows of n motifs
        index <- index + 1
        table.types[index, i] <- (sum(type[, sprintf("any_%d_%s",
                                                     n.motifs[n], species[h])] > 0) 
                                  / nrow(type)) * 100
        rownames(table.types)[index] <- paste0("%_genes_", sprintf("any_%d_%s",
                                                                   n.motifs[n], species[h]))
        # Calculate percentage of genes with windows in all orthologs
        index <- index + 1
        table.types[index, i] <- (sum(type[which(type[, sprintf("ortho_with_win_%d_mots",
                                                                n.motifs[n])] == type$num_orthologs), 
                                           sprintf("ortho_with_win_%d_mots", n.motifs[n])] != 0)
                                  / sum(type$num_orthologs != 0) * 100)
        rownames(table.types)[index] <- paste0("%_conserved_genes_", 
                                               sprintf("any_%d_%s",
                                                       n.motifs[n], species[h]))
        
        
        ## Calculate iteratively same percentages without one of the motifs
        if (n.motifs[n] != n.motifs[length(n.motifs)]){
          for (m in 1:length(motifs)){
            # Number of total windows
            index <- index + 1
            table.types[index, i] <- sum(type[, sprintf("any_%d_without_%s_%s",
                                                        n.motifs[n], motifs[m], species[h])])
            rownames(table.types)[index] <- sprintf("any_%d_without_%s_%s",
                                                    n.motifs[n], motifs[m], species[h])
            # Percentage of genes with windows
            index <- index + 1
            table.types[index, i] <- (sum(type[, sprintf("any_%d_without_%s_%s",
                                                         n.motifs[n], motifs[m], species[h])] > 0) 
                                      / nrow(type)) * 100
            rownames(table.types)[index] <- paste0("%_genes_", sprintf("any_%d_without_%s_%s",
                                                                       n.motifs[n], motifs[m], species[h]))
          }
        }
        
        
        ##### 2) Make table with information of motif composition #####
        ## Load reference table
        setwd(referenceDir)
        reference <- read.table(sprintf("gene_reference_c_%s_10k.tsv", species[h]), 
                                sep = "\t", stringsAsFactors = F, header = T)
        
        ## Load each windows file
        setwd(winDir)
        gene <- character()
        common_name <- character()
        f <- list.files()
        f <- f[grep(species[h], f)]
        f <- f[grep("700", f)]
        f <- f[grep(".csv", f)] 
        file.win <- read.table(f, sep = "\t", header = T,
                               stringsAsFactors = F)
        # Finally we keep the file as a subset with (>=6, >=7 or 8) motifs
        win <- subset(file.win, file.win$different_motifs >= n.motifs[n])
        # How many windows do we have for that species and n.motifs in those genes?
        rows <- which(reference$gene %in% strsplit(type[,species[h]], "_"))
        seqs.id <- reference$seqID[rows]
        win.type <- win[which(win$seqID %in% seqs.id),]
        if (nrow(win.type) == 0) { index <- index + 2; next }
        
        ## Complete table with % of motifs for every type
        name <- sprintf("%s_%d_%s", gsub("[.]", "_", types[i]), 
                        n.motifs[n], species[h])
        signatures <- strsplit(win.type$signature, "_")
        coup = 0
        dlx1 = 0
        dlx2 = 0
        ets = 0
        mads = 0
        meis = 0
        pax = 0
        pbx = 0
        for (j in 1:length(signatures)){
          if ("COUP" %in% signatures[[j]]) { coup = coup + 1 }
          if ("DLX1" %in% signatures[[j]]) { dlx1 = dlx1 + 1 }
          if ("DLX2" %in% signatures[[j]]) { dlx2 = dlx2 + 1 }
          if ("ETS" %in% signatures[[j]]) { ets = ets + 1 }
          if ("MADS" %in% signatures[[j]]) { mads = mads + 1 }
          if ("MEIS" %in% signatures[[j]]) { meis = meis + 1 }
          if ("PAX" %in% signatures[[j]]) { pax = pax + 1 }
          if ("PBX" %in% signatures[[j]]) { pbx = pbx + 1 }
        }
        table.motifs[(nrow(table.motifs) + 1), 1] <- name
        table.motifs[nrow(table.motifs), 2] <- (coup / length(signatures)) * 100
        table.motifs[nrow(table.motifs), 3] <- (dlx1 / length(signatures)) * 100
        table.motifs[nrow(table.motifs), 4] <- (dlx2 / length(signatures)) * 100
        table.motifs[nrow(table.motifs), 5] <- (ets / length(signatures)) * 100
        table.motifs[nrow(table.motifs), 6] <- (mads / length(signatures)) * 100
        table.motifs[nrow(table.motifs), 7] <- (meis / length(signatures)) * 100
        table.motifs[nrow(table.motifs), 8] <- (pax / length(signatures)) * 100
        table.motifs[nrow(table.motifs), 9] <- (pbx / length(signatures)) * 100
        colnames(table.motifs) <- c("type", "COUP", "DLX1", "DLX2", "ETS", 
                                    "MADS", "MEIS", "PAX", "PBX")
        
        ## Add mean and median of distance to ATG to windows table !!
        distance_ATG <- numeric()
        region <- character()
        gene <- character()
        common_name <- character ()
        for (j in 1:nrow(win.type)){
          seqID <- win.type$seqID[j]
          strand <- win.type$strand[j]
          ref <- reference[which(reference$seqID == seqID),]
          if (strand == "+"){
            if (ref$region == "upstream"){
              pos.ATG <- ref$end + ref$dist_to_first_exon
              dist.ATG <- abs(pos.ATG - win.type[j, "chrom_end"])
            } else { # ref$region == "intronic"
              pos.ATG <- ref$start - ref$dist_to_first_exon
              dist.ATG <- abs(pos.ATG - win.type[j, "chrom_start"])
            }
          } else { # strand == "-"
            if (ref$region == "upstream") {
              pos.ATG <- ref$start - ref$dist_to_first_exon
              dist.ATG <- abs(pos.ATG - win.type[j, "chrom_start"])
            } else { # ref$region == "intronic"
              pos.ATG <- ref$end + ref$dist_to_first_exon
              dist.ATG <- abs(pos.ATG - win.type[j, "chrom_end"])
            }
          }
          region <- c(region, ref$region)
          distance_ATG <- c(distance_ATG, dist.ATG)
          
          ## Create a new colum with the WBID
          gene <- c(gene, reference$gene[which(win.type$seqID[j] 
                                               == reference$seqID)])
          
          ## Create a new colum with the common name of the gene
          name <- type$common_name[which(gene[j] == type[,species[h]])]
          if (name[1] == ""){
            name <- type$common_name2[which(gene[j] == type[,species[h]])]
          }
          common_name <- c(common_name, name[1])
        }
        
        # Add mean of all win
        index <- index + 1
        table.types[index, i] <- round(mean(distance_ATG))
        rownames(table.types)[index] <- paste0("ATG_mean_distance_", 
                                               sprintf("any_%d_%s",
                                                       n.motifs[n], 
                                                       species[h]))
        # Add median of all win
        index <- index + 1
        table.types[index, i] <- round(median(distance_ATG))
        rownames(table.types)[index] <- paste0("ATG_median_distance_", 
                                               sprintf("any_%d_%s",
                                                       n.motifs[n], 
                                                       species[h]))
        
        
        ##### 3) Make table with information of windows for each gene #####
        ## Add the column with the WBID
        win.type <- tibble::add_column(win.type, gene, .before = 1)
        ## Add the colum with the common name of the gene
        win.type <- tibble::add_column(win.type, common_name, .before = 2)
        
        ## Initialize table.gene as a subset of win.type
        table.gene <- win.type[, c("gene", "common_name", "chrom", 
                                   "chrom_start", "chrom_end", 
                                   "window_width")]
        # Get the number of windows for each gene
        num_win <- numeric()
        for (j in 1:length(table.gene$gene)){
          num_win <- c(num_win, length(which(table.gene$gene[j] 
                                             == win.type$gene)))
        }
        table.gene <- tibble::add_column(table.gene, num_win, 
                                         .before = "chrom")
        colnames(table.gene)[which(colnames(table.gene)
                                   == "num_win")] <- sprintf("num_win_%d_mots", 
                                                             n.motifs[n])
        # Get the frequency of motifs for each win
        signatures <- strsplit(win.type$signature, "_")
        for (j in 1:length(signatures)){
          coup = 0
          dlx1 = 0
          dlx2 = 0
          ets = 0
          mads = 0
          meis = 0
          pax = 0
          pbx = 0
          if ("COUP" %in% signatures[[j]]) { coup = sum("COUP" == signatures[[j]]) }
          if ("DLX1" %in% signatures[[j]]) { dlx1 = sum("DLX1" == signatures[[j]]) }
          if ("DLX2" %in% signatures[[j]]) { dlx2 = sum("DLX2" == signatures[[j]]) }
          if ("ETS" %in% signatures[[j]])  { ets  = sum("ETS" == signatures[[j]]) }
          if ("MADS" %in% signatures[[j]]) { mads = sum("MADS" == signatures[[j]]) }
          if ("MEIS" %in% signatures[[j]]) { meis = sum("MEIS" == signatures[[j]]) }
          if ("PAX" %in% signatures[[j]])  { pax  = sum("PAX" == signatures[[j]]) }
          if ("PBX" %in% signatures[[j]])  { pbx  = sum("PBX" == signatures[[j]]) }
          table.gene[j, "n_motifs_win"] <- length(signatures[[j]])
          table.gene[j, "%_COUP"] <- coup / length(signatures[[j]]) * 100 # % COUP
          table.gene[j, "%_DLX1"] <- dlx1 / length(signatures[[j]]) * 100
          table.gene[j, "%_DLX2"] <- dlx2 / length(signatures[[j]]) * 100
          table.gene[j, "%_ETS"] <- ets / length(signatures[[j]]) * 100
          table.gene[j, "%_MADS"] <- mads / length(signatures[[j]]) * 100
          table.gene[j, "%_MEIS"] <- meis / length(signatures[[j]]) * 100
          table.gene[j, "%_PAX"] <- pax / length(signatures[[j]]) * 100
          table.gene[j, "%_PBX"] <- pbx / length(signatures[[j]]) * 100
          if (sum(table.gene[j, 9:(9 + (length(motifs) - 1))]) == 100){ next } # check
        }
        ## Add column with distance to ATG and the type of region (upstream/intronic)
        table.gene <- tibble::add_column(table.gene, distance_ATG,
                                         .after = "window_width")
        table.gene <- tibble::add_column(table.gene, region,
                                         .after = "distance_ATG")
        
        ## Add signature column to final table.gene
        # This column is alphabetically ordered in the original windows file, so we are going
        # to reorder it using the $starts column to know the real order of the motifs in the seq
        signatures <- strsplit(win.type$signature, split = "_")
        starts <- strsplit(win.type$starts, split = "_")
        for (j in 1:nrow(win.type)){
          df <- as.data.frame(signatures[[j]], 
                              starts[[j]])
          sig <- as.character(df[sort(rownames(df)), 1])
          sig <- paste0(sig, sep = "_", collapse = "")
          sig <- substr(sig, 1, nchar(sig) - 1)
          if (j == 1){ # New column
            table.gene[j, ncol(table.gene) + 1] <- sig
          } else {
            table.gene[j, ncol(table.gene)] <- sig
          }
        }
        colnames(table.gene)[ncol(table.gene)] <- sprintf("signature_%d_mots",
                                                          n.motifs[n])
        
        ## When finished the table.genes for each type, write tables:
        filename <- paste0(newDir,
                           sprintf("/%s_genes_%s_%d_mots.tsv", 
                                   gsub("[.]", "_", types[i]), 
                                   species[h], n.motifs[n]))
        write.table(table.gene, file = filename, 
                    sep = "\t", row.names = F, quote = F)
        cat(sprintf(">> Finalizado anotación de %s en C. %s con %d motivos.\n", 
                    types[i], species[h], n.motifs[n]))
        
        
        ##### 4) Get logos from conserved C. elegans genes #####
        if (h == 1){
          cat(sprintf(">> Haciendo logo de %s con %d motivos...\n", 
                      types[i], n.motifs[n]))
          ## How many conserved genes with windows?
          conserved <- type[which(type[, sprintf("ortho_with_win_%d_mots",
                                                 n.motifs[n])] == type$num_orthologs),]
          if (nrow(conserved) == 0) { next }
          
          ## Load reference table
          setwd(referenceDir)
          reference <- read.table(sprintf("gene_reference_c_%s_10k.tsv",
                                          species[h]), 
                                  sep = "\t", 
                                  stringsAsFactors = F, 
                                  header = T)
          
          ## Load windows file
          ## Load each windows file
          setwd(winDir)
          gene <- character()
          common_name <- character()
          f <- list.files()
          f <- f[grep(species[h], f)]
          f <- f[grep("700", f)]
          f <- f[grep(".csv", f)] 
          file.win <- read.table(f, sep = "\t", header = T, 
                                 stringsAsFactors = F)
          # Finally we keep the file as a subset with (>=6, >=7, or 8) motifs
          win <- subset(file.win, file.win$different_motifs 
                        >= n.motifs[n])
          # How many windows do we have for that species and n.motifs in those genes?
          rows <- which(reference$gene %in% type[,species[h]])
          seqs.id <- reference$seqID[rows]
          win.type <- win[which(win$seqID %in% seqs.id),]
          ## Add a new colum with the WBID
          gene <- character()
          for (j in 1:nrow(win.type)){
            gene <- c(gene, reference$gene[which(win.type$seqID[j] 
                                                 == reference$seqID)])
          }
          win.type <- tibble::add_column(win.type, gene, .before = 1)
          
          ## Add a new colum with the common name of the gene
          common_name <- character()
          for (j in 1:nrow(win.type)){
            name <- type$common_name[which(win.type$gene[j]
                                           == type[,species[h]])]
            if (name == ""){
              name <- type$common_name2[which(win.type$gene[j] 
                                              == type[,species[h]])]
            }
            common_name <- c(common_name, name[1])
          }
          win.type <- tibble::add_column(win.type, common_name,
                                         .before = 2)
          
          ## Some genes have a n motif window in all 4 species but not in elegans.
          # Those genes are:
          table((conserved$elegans %in% unique(win.type$gene)))
          conserved[(!(conserved$elegans %in% unique(win.type$gene))), 
                    c("common_name", "any_6_elegans", "any_7_elegans", 
                      "any_8_elegans", "ortho_with_win_8_mots")]
          
          # We drop them from the analysis. 
          conserved <- conserved[(conserved$elegans %in%
                                    unique(win.type$gene)), ]
          
          ## Load seqs file
          setwd(seqsDir)
          f <- list.files()
          file <- f[grep("elegans", f)]
          all.seqs <- read.table(file, sep = "\t", header = T, 
                                 stringsAsFactors = F)
          
          ## Select only conserved genes
          win.type <- win.type[which(win.type$gene %in% 
                                       conserved$elegans),]
          if (nrow(win.type) == 0) { next }
          # Get seqs, signatures, starts and ends for each window
          seqs <- strsplit(win.type$seqs, "_")
          signatures <- strsplit(win.type$signature, "_")
          positions_start <- strsplit(win.type$starts, "_")
          positions_end <- strsplit(win.type$ends, "_")
          strands <- strsplit(win.type$strands, "_")
          
          ## Create a column for each motif that may be in the window
          table.pwm <- win.type[, c(1, 2, 3)]
          for (j in 1:length(motifs)){
            table.pwm[, (3 + j)] <- ""
            colnames(table.pwm)[3 + j] <- motifs[j]
          }
          ## Select for each motif +- 10 pb
          COUP <- character()
          DLX1 <- character()
          DLX2 <- character()
          ETS <- character()
          MADS <- character()
          MEIS <- character()
          PAX <- character()
          PBX <- character()
          for (j in 1:nrow(win.type)){ # each window
            for (k in 1:length(signatures[[j]])){ # each match
              # Select the position of the match
              start <- as.integer(positions_start[[j]][k])
              end <- as.integer(positions_end[[j]][k])
              # If we have already count that match, we skip the iteration 
              if (j > 1) { if (start %in% positions_start[[j - 1]]) { next } }
              
              # Extend the match 10 pb to each side of the sequence
              seq <- all.seqs[which(all.seqs$seqID == table.pwm$seqID[j]),
                              "original.seq"]
              motif <- signatures[[j]][k]
              len <- as.integer(positions_end[[j]][k]) - as.integer(positions_start[[j]][[k]]) + 1 + 20
              seq <- substr(seq, (start - 10), (end + 10)) 
              # Proof to possible cases in which we cannot extend the sequence
              if (nchar(seq) != len) { next } 
              
              # If there's more than one signature of the same motif 
              # in that seq, we concatenate them by "_" 
              if (table.pwm[j, motif] == "") {
                table.pwm[j, motif] <- seq
              } else {
                table.pwm[j, motif] <- paste(table.pwm[j, motif],
                                             seq, sep = "_")
              }
              # If the match is on the negtive strand, we reverse the sequence
              if (strands[[j]][k] == "-") {
                seq <- as.character(reverseComplement(DNAString(seq)))
              }
              # Create character vector for each motif
              object <- eval(parse(text = motif))
              object <- c(object, seq)
              assign(motif, object)
            }
          }
          
          ## Create and save logos
          setwd(newDir)
          if (!(dir.exists("logo_plots"))){
            dir.create("logo_plots")
          }
          if (!(dir.exists(sprintf("logo_plots/%d_mots", n.motifs[n])))){
            dir.create(sprintf("logo_plots/%d_mots", n.motifs[n]))
          }
          
          for (j in 1:length(motifs)){
            logo <- eval(parse(text = motifs[j]))
            if (length(logo) == 0) { next }
            logo <- ggseqlogo(logo) +
              ggtitle(sprintf("%s logo from %d sequences", 
                              motifs[j], length(logo)))
            filename <- sprintf("./logo_plots/%d_mots/%s_logo_%s_%s_%d_motifs.png", 
                                n.motifs[n], motifs[j], 
                                types[i], species[h], n.motifs[n])
            ggsave(filename, plot = logo,
                   width = 300, height = 100, units = "mm")
          }
        } else { next } # if species is not elegans, don't make the logos
      }
    }
    ## Save type genes for future analysis
    filename <- paste0(newDir, sprintf("/%s_genes.Rda", 
                                       gsub("[.]", "_", types[i])))
    save(type, file = filename)
  }
  
  ## Write final table.types:
  setwd(newDir)
  filename <- "./percentage_win_all_types.tsv"
  write.table(table.types, file = filename, sep = "\t", 
              row.names = T, col.names = F, quote = F)
  
  ## Write final table.motifs:
  filename <- "./percentage_motifs_all_types.tsv"
  write.table(table.motifs, file = filename, sep = "\t", 
              row.names = F, quote = F)
  
} else {
  load("DA_genes.Rda")
}





######################################################
##         6. Plots of results analysis win         ##
######################################################
#######################################################
##   Results of the worm19 + random genes analysis   ##
#######################################################
## Clean environment and set working directory
rm(list = ls())

parentDir <- "/media/erick/OS/Bioinformatics/TFM/DA/worm19/random_genes"
setwd(parentDir)

## Load needed packages
library(ggplot2)
library(RColorBrewer)

## Load complete gene summary file
f <- list.files()
file <- f[grep("190522_worm19_complete_gene_summary_windows.tsv", f)]
gene.summary <- read.table(file, sep = "\t",
                           stringsAsFactors = F, 
                           header = T) 

# How many types are we getting?
types <- sort(unique(gene.summary$type[gene.summary$type != "other"
                                       & gene.summary$type != "neuronal"])) 
conservation.status <- 4
n.motifs <- c(6, 7, 8)

## For each type and number of motifs, do a plot.
for (n in n.motifs) {
  for (i in 1:length(types)){
    ## The background random which is selected in each moment is defined by file
    # All random genes
    file <- paste0("./all_types_windows/all_random/", # select background
                   types[i], "_all_random_groups_WithoutOverlaps.tsv") 
    all.random <- read.table(file, sep = "\t",
                             stringsAsFactors = F, header = T)
    colnames(all.random) <- c("gene", 
                              paste(n.motifs[length(n.motifs)]:n.motifs[1],
                                    "mots", sep="_"), "conservation")
    
    # Conserved random
    file <- paste0("./all_types_windows/conserved_random/", # select background
                   types[i], "_all_random_groups_WithoutOverlaps.tsv")
    conserved.random <- read.table(file, sep = "\t",
                                   stringsAsFactors = F, header = T)
    colnames(conserved.random) <- c("gene", 
                                    paste(n.motifs[length(n.motifs)]:n.motifs[1],
                                          "mots", sep="_"), "conservation")
    
    ## Select DA genes (actually type, not only DA)
    type  <- gene.summary[gene.summary$type == types[i],]
    assign(types[i], type)
    
    ## Classificate type genes according to the % of orthologous genes
    ## that have windows with n motifs
    conservation <-  table(type$num_orthologs[which(type$num_orthologs != 0)])
    type.genes <- data.frame()
    mots <- paste(n.motifs[1]:n.motifs[length(n.motifs)], "mots", sep="_")
    ## For the different status of conservation:
    rows.conserved <-  which(type$num_orthologs != 0)
    rows.conserved.win <- which(type[rows.conserved, 
                                     sprintf("ortho_with_win_%d_mots", n)]
                                == type$num_orthologs[rows.conserved])
    rows <- rows.conserved[which(rows.conserved %in% rows.conserved.win)]
    genes <- type[rows, "elegans"]
    genes <- paste(genes, collapse = "_")
    type.genes[nrow(type.genes) + 1,"genes"] <- genes
    for(m in mots){
      type.genes[nrow(type.genes), m] <- length(rows.conserved.win) / length(rows.conserved)
    }
    type.genes[nrow(type.genes), "conservation"] <- "conserved"
    
    ## Add info of all genes
    type.genes[nrow(type.genes) + 1, "conservation"] <- "any"
    rows <- c(1:nrow(type))
    genes <- type[rows, "elegans"]
    genes <- paste(genes, collapse = "_")
    type.genes[nrow(type.genes), 1] <- genes
    for(m in mots){
      col <- which(colnames(type) == sprintf("any_%s_elegans", 
                                             strsplit(m, "_")[[1]][1]))
      var <- type[rows, col] > 0
      type.genes[nrow(type.genes), m] <- sum(var)/length(rows)
    }
    
    
    #######################################################
    ## Plot distributions: all genes and conserved genes ##
    #######################################################
    ## Create dir if it not exists
    if (!(dir.exists("plots_windows"))) { 
      system("mkdir plots_windows")
      system("mkdir plots_windows/violin_plots")
    }
    
    ## How are our genes of each type distributed?
    summary(type.genes)
    all.score <- type.genes[which(type.genes$conservation == "any"), 
                            sprintf("%d_mots", n)]
    type.score <- type.genes[which(type.genes$conservation == "conserved"),
                             sprintf("%d_mots", n)]
    
    ## How are the random genes distributed?
    summary(all.random)
    summary(conserved.random)
    conserved.genes <- conserved.random[which(conserved.random$conservation == conservation.status),]
    summary(conserved.genes)
    conserved.genes <- sample(conserved.genes[, sprintf("%d_mots", n)], 
                              10000) 
    
    all.genes <- sample(all.random[, sprintf("%d_mots", n)], 
                        length(conserved.genes)) 
    
    ## Make data for plot
    scores <- c(all.genes, conserved.genes)
    category <- c(rep("all.genes", length(all.genes)), 
                  rep("conserved.genes", length(conserved.genes)))
    datos <- data.frame(scores = scores,
                        category = as.factor(category),
                        type = c(rep(all.score, 
                                     sum(category == "all.genes")),
                                 rep(type.score, 
                                     sum(category == "conserved.genes"))))
    
    ## Violin plot of type
    ## And distribution of all.genes and conserved.genes
    ## We can select palettes with:
    ## RColorBrewer::display.brewer.all()
    palette <- RColorBrewer::brewer.pal(name = "Paired", n=3)
    
    title <- sprintf("%s genes over random background with %d motifs",
                     gsub("[.]", " ", types[i]), n)
    p <- ggplot(datos, aes(x = category, y = scores,
                           fill = category))
    p <- p + geom_violin(trim = T, bw = 0.01, adjust = 3) + 
      geom_boxplot(width = 0.1, fill = "white") +
      geom_point(aes(x = category, y = type, color = "type"),
                 shape = 16, size = 4.5) +
      scale_color_manual(values = "#FD626F",
                         labels = c(sprintf("%s genes", 
                                            gsub("[.]", " ", types[i])))) +
      scale_fill_manual(values = palette[1:2], guide = FALSE) +
      scale_x_discrete(labels = c("All genes",
                                  "Conserved genes")) +
      scale_y_continuous(name = "% of total genes",
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      labs(color = "Gene expression") +
      ggtitle(title) +
      theme(axis.text.x = element_text(color = "#000000",
                                       size = 10),
            axis.title.x = element_text(size = 13),
            axis.text.y = element_text(color = "#000000",
                                       size = 10),
            axis.title.y = element_text(size = 13),
            plot.title = element_text(hjust = 0.5,
                                      size = 17),
            legend.position = "top",
            legend.key = element_blank(),
            panel.background = element_rect(fill = "white")) 
    
    p
    filename <- sprintf("./plots_windows/violin_plots/%s_expression_over_random_background_%d_mots.png",
                        gsub("[.]", "_", types[i]), n)
    ggsave(filename, plot = p,
           width = 200, height = 250, units = "mm")
    
  }
}


#######################################################
##     Make violin plots comparison for all types    ##
#######################################################
## Clean environment and set working directory
rm(list = ls())

setwd("/media/erick/OS/Bioinformatics/TFM/DA/worm19/random_genes")

## Load needed packages
library(ggplot2)
library(RColorBrewer)

## Load complete gene summary file
f <- list.files()
file <- f[grep("190522_worm19_complete_gene_summary_windows.tsv", f)]
gene.summary <- read.table(file, sep="\t", stringsAsFactors=F, header=T) 

types <- c("ASE", "DA", "GABAergic", "RIA", "SDQ.ALN.PLN", "Touch.receptor")
conservation.status <- 4
n.motifs <- c(6, 7, 8)
categories <- c("all", "conserved")

## For each type and number of motifs, do a plot.
for (h in categories) {
  for (n in n.motifs) {
    all.data <- data.frame()
    for (i in 1:length(types)){
      ## The background random which is selected in each moment is defined by file
      # All random genes
      file <- paste0("./all_types_windows/all_random/", # select background
                     types[i], "_all_random_groups_WithoutOverlaps.tsv") 
      all.random <- read.table(file, sep = "\t",
                               stringsAsFactors = F, header = T)
      colnames(all.random) <- c("gene", 
                                paste(n.motifs[length(n.motifs)]:n.motifs[1], "mots", sep="_"),
                                "conservation")
      
      # Conserved random
      file <- paste0("./all_types_windows/conserved_random/", # select background
                     types[i], "_all_random_groups_WithoutOverlaps.tsv")
      conserved.random <- read.table(file, sep = "\t",
                                     stringsAsFactors = F, header = T)
      colnames(conserved.random) <- c("gene", 
                                      paste(n.motifs[length(n.motifs)]:n.motifs[1], "mots", sep="_"),
                                      "conservation")
      
      ## Select DA genes (actually type, not only DA)
      type  <- gene.summary[gene.summary$type == types[i],]
      assign(types[i], type)
      
      ## Classificate type genes according to the % of orthologous genes
      ## that have windows with n motifs
      conservation <-  table(type$num_orthologs[which(type$num_orthologs != 0)])
      type.genes <- data.frame()
      mots <- paste(n.motifs[1]:n.motifs[length(n.motifs)], "mots", sep="_")
      ## For the different status of conservation:
      rows.conserved <-  which(type$num_orthologs != 0)
      rows.conserved.win <- which(type[rows.conserved, sprintf("ortho_with_win_%d_mots", n)] == type$num_orthologs[rows.conserved])
      rows <- rows.conserved[which(rows.conserved %in% rows.conserved.win)]
      genes <- type[rows, "elegans"]
      genes <- paste(genes, collapse = "_")
      type.genes[nrow(type.genes) + 1,"genes"] <- genes
      for(m in mots){
        type.genes[nrow(type.genes), m] <- length(rows.conserved.win) / length(rows.conserved)
      }
      type.genes[nrow(type.genes), "conservation"] <- "conserved"
      
      ## Add info of all genes
      type.genes[nrow(type.genes) + 1, "conservation"] <- "any"
      rows <- c(1:nrow(type))
      genes <- type[rows, "elegans"]
      genes <- paste(genes, collapse = "_")
      type.genes[nrow(type.genes), 1] <- genes
      for(m in mots){
        col <- which(colnames(type) == sprintf("any_%s_elegans", strsplit(m, "_")[[1]][1]))
        var <- type[rows, col] > 0
        type.genes[nrow(type.genes), m] <- sum(var)/length(rows)
      }
      
      
      #######################################################
      ## Plot distributions: all genes and conserved genes ##
      #######################################################
      ## Create dir if it not exists
      if (!(dir.exists("plots_windows"))) { 
        system("mkdir plots_windows")
        system("mkdir plots_windows/violin_plots")
      }
      
      ## How are our genes of each type distributed?
      summary(type.genes)
      all.score <- type.genes[which(type.genes$conservation == "any"), sprintf("%d_mots", n)]
      type.score <- type.genes[which(type.genes$conservation == "conserved"), sprintf("%d_mots", n)]
      
      ## How are the random genes distributed?
      conserved.genes <- conserved.random[which(conserved.random$conservation == conservation.status),]
      conserved.genes <- sample(conserved.genes[, sprintf("%d_mots", n)],  nrow(conserved.genes))
      all.genes <- sample(all.random[, sprintf("%d_mots", n)], length(conserved.genes))
      
      ## Make data for plot
      scores <- c(all.genes, conserved.genes)
      category <- c(rep("all.genes", length(all.genes)), rep("conserved.genes", length(conserved.genes)))
      datos <- data.frame(scores = scores,
                          category = as.factor(category),
                          type = c(rep(all.score, sum(category == "all.genes")),
                                   rep(type.score, sum(category == "conserved.genes"))))
      
      ## Make new df for all violin plots comparison
      new.df <- subset(datos, category == paste(h, "genes", sep = "."))
      new.df[,"category"] <- types[i]
      if (nrow(all.data) != 0){
        all.data <- rbind(all.data, new.df)
      } else { all.data <- new.df }
    }
    ## Retouch names
    all.data$category[which(all.data$category == "Touch.receptor")] <- gsub("[.]", " ", 
                                                                            all.data$category[which(all.data$category == "Touch.receptor")])
    all.data$category[which(all.data$category == "SDQ.ALN.PLN" )] <- gsub("[.]", "/", 
                                                                          all.data$category[which(all.data$category == "SDQ.ALN.PLN" )])
    
    ## Reorder data
    df <- data.frame()
    df <- all.data[which(all.data$category == "DA"), ]
    df <- rbind(df, all.data[which(all.data$category == "RIA"),] )
    df <- rbind(df, all.data[which(all.data$category == "ASE"),] )
    df <- rbind(df, all.data[which(all.data$category == "Touch receptor"),] )
    df <- rbind(df, all.data[which(all.data$category == "GABAergic"),] )
    df <- rbind(df, all.data[which(all.data$category == "SDQ/ALN/PLN"),] )
    rownames(df) <- 1:nrow(df)
    
    ## Violin plot of all neuronal types
    # We can select palettes with:
    # RColorBrewer::display.brewer.all()
    palette <- RColorBrewer::brewer.pal(name="Set3", n=10)
    
    title <- sprintf("Motif signature in %s genes over random background with %d motifs", h, n)
    p <- ggplot(df, aes(x = category, y = scores,
                        fill = category))
    p <- p + geom_violin(trim = T, bw = 0.01, adjust = 3) + 
      geom_boxplot(width = 0.1, fill = "white") +
      geom_point(aes(x = category, 
                     y = type, 
                     color = "type"),
                 shape = 16, size = 4.5) +
      scale_color_manual(values = "#FD626F",   #"#E22121",
                         labels = "% of genes with windows") +
      scale_fill_manual(values = c(palette[9], 
                                   palette[1], 
                                   rep(palette[9], 4)), 
                        guide = FALSE) +
      scale_x_discrete(limits = unique(df$category),
                       labels = paste(unique(df$category), 
                                      "neurons", sep = "\n")) +
      scale_y_continuous(name = "% of total genes",
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      labs(color = "") +
      ggtitle(title) +
      theme(axis.text.x = element_text(color = "#000000",
                                       size = 10),
            axis.title.x = element_text(size = 13),
            axis.text.y = element_text(color = "#000000",
                                       size = 10),
            axis.title.y = element_text(size = 13),
            plot.title = element_text(hjust = 0.5,
                                      size = 18),
            legend.position = "top",
            legend.key = element_blank(),
            panel.background = element_rect(fill = "white"))  
    p
    filename <- sprintf("./plots_windows/violin_plots/%s_types_over_random_background_%d_mots.png", 
                        R.utils::capitalize(h), n)
    ggsave(filename, plot = p,
           width = 33, height = 20, units = "cm")
    
    ## Print percentile of the value for each type:
    if (n == 8){
      print(sprintf("%s genes comparison with %d motifs:", R.utils::capitalize(h), n))
      classes <- unique(df$category)
      for (class in classes){
        print(sprintf(">> %s percentile:", class))
        class.scores <- df[which(df$category == class),]
        per <- ecdf(class.scores$scores)(unique(class.scores$type))
        print(per)
        cat("\n")
      }
    }
  }
}


#######################################################
##   Results of the worm19 + random genes analysis   ##
#######################################################
## Clean environment and set working directory
rm(list = ls())

# parentDir <- "/home/esousa/TFM/DA/worm19/random_genes"
parentDir <- "/media/erick/OS/Bioinformatics/TFM/DA/worm19/random_genes"
setwd(parentDir)

## Load needed packages
library(ggplot2)
library(RColorBrewer)

## Load complete gene summary file
f <- list.files()
file <- f[grep("190522_worm19_complete_gene_summary_windows.tsv", f)]
gene.summary <- read.table(file, sep = "\t",
                           stringsAsFactors = F, 
                           header = T) 
# How many types are we getting?
types <- "DA"; i = 1 # only DA
conservation.status <- 4
n.motifs <- c(6, 7, 8)
categories <- c("all", "conserved")


## For each type and number of motifs, do a plot.
for (h in categories) {
  all.data <- data.frame()
  for (n in n.motifs) {
    ## The background random which is selected in each moment is defined by file
    # All random genes
    file <- paste0("./all_types_windows/all_random/", # select background
                   types[i], "_all_random_groups_WithoutOverlaps.tsv") 
    all.random <- read.table(file, sep = "\t",
                             stringsAsFactors = F, header = T)
    colnames(all.random) <- c("gene", 
                              paste(n.motifs[length(n.motifs)]:n.motifs[1], "mots", sep="_"),
                              "conservation")
    
    # Conserved random
    file <- paste0("./all_types_windows/conserved_random/", # select background
                   types[i], "_all_random_groups_WithoutOverlaps.tsv")
    conserved.random <- read.table(file, sep = "\t",
                                   stringsAsFactors = F, header = T)
    colnames(conserved.random) <- c("gene", 
                                    paste(n.motifs[length(n.motifs)]:n.motifs[1], "mots", sep="_"),
                                    "conservation")
    
    ## Select DA genes (actually type, not only DA)
    type  <- gene.summary[gene.summary$type == types[i],]
    assign(types[i], type)
    
    ## Classificate type genes according to the % of orthologous genes
    ## that have windows with n motifs
    conservation <-  table(type$num_orthologs[which(type$num_orthologs != 0)])
    type.genes <- data.frame()
    mots <- paste(n.motifs[1]:n.motifs[length(n.motifs)], "mots", sep="_")
    ## For the different status of conservation:
    rows.conserved <-  which(type$num_orthologs != 0)
    rows.conserved.win <- which(type[rows.conserved, sprintf("ortho_with_win_%d_mots", n)] == type$num_orthologs[rows.conserved])
    rows <- rows.conserved[which(rows.conserved %in% rows.conserved.win)]
    genes <- type[rows, "elegans"]
    genes <- paste(genes, collapse = "_")
    type.genes[nrow(type.genes) + 1,"genes"] <- genes
    for(m in mots){
      type.genes[nrow(type.genes), m] <- length(rows.conserved.win) / length(rows.conserved)
    }
    type.genes[nrow(type.genes), "conservation"] <- "conserved"
    
    ## Add info of all genes
    type.genes[nrow(type.genes) + 1, "conservation"] <- "any"
    rows <- c(1:nrow(type))
    genes <- type[rows, "elegans"]
    genes <- paste(genes, collapse = "_")
    type.genes[nrow(type.genes), 1] <- genes
    for(m in mots){
      col <- which(colnames(type) == sprintf("any_%s_elegans", strsplit(m, "_")[[1]][1]))
      var <- type[rows, col] > 0
      type.genes[nrow(type.genes), m] <- sum(var)/length(rows)
    }
    
    
    #######################################################
    ## Plot distributions: all genes and conserved genes ##
    #######################################################
    ## Create dir if it not exists
    if (!(dir.exists("plots_windows"))) { 
      system("mkdir plots_windows")
      system("mkdir plots_windows/violin_plots")
    }
    
    ## How are our genes of each type distributed?
    summary(type.genes)
    all.score <- type.genes[which(type.genes$conservation == "any"), sprintf("%d_mots", n)]
    type.score <- type.genes[which(type.genes$conservation == "conserved"), sprintf("%d_mots", n)]
    
    ## How are the random genes distributed?
    conserved.genes <- conserved.random[which(conserved.random$conservation == conservation.status),]
    conserved.genes <- sample(conserved.genes[, sprintf("%d_mots", n)],  nrow(conserved.genes))
    all.genes <- sample(all.random[, sprintf("%d_mots", n)], length(conserved.genes))
    
    ## Make data for plot
    scores <- c(all.genes, conserved.genes)
    category <- c(rep("all.genes", length(all.genes)), rep("conserved.genes", length(conserved.genes)))
    datos <- data.frame(scores = scores,
                        category = as.factor(category),
                        type = c(rep(all.score, sum(category == "all.genes")),
                                 rep(type.score, sum(category == "conserved.genes"))))
    
    
    ## Make new df for all violin plots comparison
    new.df <- subset(datos, category == paste(h, "genes", sep = "."))
    new.df[,"category"] <- sprintf("%d motif signature", n)
    if (nrow(all.data) != 0){
      all.data <- rbind(all.data, new.df)
    } else { all.data <- new.df }
  }
  
  ## Reorder data
  df <- data.frame()
  df <- all.data[which(all.data$category == "8 motif signature"), ]
  df <- rbind(df, all.data[which(all.data$category == "7 motif signature"),] )
  df <- rbind(df, all.data[which(all.data$category == "6 motif signature"),] )
  rownames(df) <- 1:nrow(df)
  
  ## Violin plot of all neuronal types
  # We can select palettes with:
  # RColorBrewer::display.brewer.all()
  palette <- RColorBrewer::brewer.pal(name="Set3", n=10)
  
  if (h == "all") { title <- "DA signature" 
  } else { title <- sprintf("%s DA signature", R.utils::capitalize(h)) }
  
  p <- ggplot(df, aes(x = category, y = scores,
                      fill = category))
  p <- p + geom_violin(trim = T, bw = 0.01, adjust = 3) +  
    geom_boxplot(width = 0.1, fill = "white") +
    geom_point(aes(x = category, 
                   y = type, 
                   color = "type"),
               shape = 16, size = 4.5) +
    scale_color_manual(values = "#FD626F",   #"#E22121",
                       labels = "% of genes with windows") +
    scale_fill_manual(values = rep(palette[1], length(n.motifs)), 
                      guide = FALSE) +
    scale_x_discrete(limits = unique(df$category),
                     labels = unique(df$category)) +
    scale_y_continuous(name = "% of total genes",
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(color = "",
         x = "") +
    ggtitle(title) +
    theme(axis.text.x = element_text(color = "#000000",
                                     size = 10),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(color = "#000000",
                                     size = 10),
          axis.title.y = element_text(size = 13),
          plot.title = element_text(hjust = 0.5,
                                    size = 18),
          legend.position = "top",
          legend.key = element_blank(),
          panel.background = element_rect(fill = "white"))  
  p
  filename <- sprintf("./plots_windows/violin_plots/%s_DA_signature.png", 
                      R.utils::capitalize(h))
  ggsave(filename, plot = p,
         width = 30, height = 20, units = "cm")
  
  print(sprintf("%s DA genes comparison with %d motifs:", R.utils::capitalize(h), n))
  classes <- unique(df$category)
  for (class in classes){
    print(sprintf(">> %s percentile:", class))
    class.scores <- df[which(df$category == class),]
    per <- ecdf(class.scores$scores)(unique(class.scores$type))
    print(per)
    cat("\n")
  }
}

