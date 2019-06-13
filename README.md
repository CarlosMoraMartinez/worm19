# worm19
This repository has been created for my Master's dissertation. It contains all the tools that have been used for the analysis of the regulatory dopaminergic signature of *C. elegans*.

Some of the scripts are intended to be used individually, because they expect already downloaded data (the ChIP-seq scripts) or they download their own data (the Single Cell one). Most of them are intended to work as a pipeline, in this precise order:
1. Item 1
2. Item 2
3. Item 3

A brief description of the scripts is provided here:
- First, the file "Analysis_Single_Cell_Cao.R" uses the data from Cao et al. (2017) to get a list of DEG in dopaminergic neurons. It also gets lists of DEG in other neuronal types: ASE, GABAergic, RIA, SDQ/ALN/PLN and Touch receptor neurons. These gene lists are intended to be use as a control when looking for the dopaminergic signature.
- The file "ChIP-seq_comparison.R" works with BED files downloaded from ENCODE, and loads them from the computer directly before analysing them.
- The file "Split_multiFASTA_into_chromosomes.py" is used to split the .fna file of the WBcel235 genome version into individual .fa files, each of them corresponding to a chromosome.
- The file "Get_peak_sequences_from_BED.py" uses a given BED file and the genome files split into chromosomes (from the previous script) as an input; and writes a FASTA file with the sequences corresponding to the peaks of the BED file.
- ...
