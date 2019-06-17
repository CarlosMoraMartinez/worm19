# worm19
This repository has been created for my Master's dissertation in Bioinformatics. It contains all the tools that have been used for the analysis of the regulatory dopaminergic signature of *C. elegans*.

Some of the scripts are intended to be used individually, because they expect already downloaded data (the ChIP-seq scripts) or they download their own data (the Single Cell script). The main script is **w19_complete_pipeline.R**, and it contains the whole pipeline of scripts that have been used for the search of regulatory windows and the posterior analysis of the data. This script relies on data, such as the PWM or the reference sequences, that can be found into the **data** folder.

A brief description of the scripts is provided here:
- The folder **w18** contains the R package for the search of the selected transcription-factor binding motifs. Inside you can found a script named **make_w18_package.R** with instructions on how to install it.
- First, the file **Analysis_Single_Cell_Cao.R** uses the data from Cao et al. (2017) to get a list of DEG in dopaminergic neurons. It also gets lists of DEG in other neuronal types: ASE, GABAergic, RIA, SDQ/ALN/PLN and Touch receptor neurons. These gene lists are intended to be use as a control when looking for the dopaminergic signature.
- The file **ChIP-seq_comparison.R** works with BED files downloaded from ENCODE, and loads them from the computer directly before analysing them.
- The file **Split_multiFASTA_into_chromosomes.py** is used to split the .fna (multiFASTA) file of the WBcel235 genome version into individual FASTA files, each of them corresponding to a chromosome.
- The file **Get_peak_sequences_from_BED.py** uses a given BED file and the genome FASTA files (split into chromosomes from the previous script) as an input; and writes a FASTA file with the sequences corresponding to the peaks of the BED file.
- The file **w19_complete_pipeline.R**...
