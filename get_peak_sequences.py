##### get_peak_sequences.py
"""
This script gets the sequences from the peaks on a bed file.

Autor: Erick Sousa
Created: 04/03/2019
"""
## Load modules
import os
path = os.getcwd()

## Ask for the BED file
print ("#" * 70)
print ("This script gets the sequences from each peak on a BED file".center(70))
print ("and writes them on a .txt file".center(70))
print ("#" * 70)

while True:
    try:
        gene = str(input("\n> Enter ChIP-seq target gene name: "))
        print("\nChoose ChIP-seq stage from the following: "
        "embryo, L1, L2, L3, L4, YA.")
        stage = str(input("> Enter stage: "))

        ## Load bed file
        print ("\nLoading BED file...")
        chip_path = path + "/" + gene
        file = gene + "_" + stage + ".bed"
        f_in = open((chip_path + "/" + file), "r")
        break
    except FileNotFoundError:
        print ("\nFile not found! Please check file and try again.")

bed_file = f_in.read()
f_in.close()
bed_file = bed_file.splitlines(False)
bed = []
for i in bed_file:
    line = i.split("\t")
    bed.append(line)

## Load chromosomes
print ("Loading genome files...")
genome_path = path + "/" + "WBcel235" + "/" + "FASTA_chr"
files = os.listdir(genome_path)
genome = []
for i in files:
    f_in = open(genome_path + "/" + i, "r")
    file = f_in.readlines()
    f_in.close()
    lineas = []
    for j in range(len(file)):
        lineas.append(file[j].strip())
    header = lineas[0]
    seq = "".join(lineas[1:])
    genome.append(seq)
# Now genome is a list which contains all chromosomes:
# chrI, chrII, chrIII, chrIV, chrV, chrX

## Iterate through bed peaks
print ("Getting bed peaks...")
headers = []
seqs = []
for i in range(len(bed)):
    chromosome = bed[i][0]
    start = int(bed[i][1]) - 1
    end = int(bed[i][2])
    header = ">%s:%d-%d " %(chromosome, start, end)
    if header not in headers:
        headers.append(header)

        if chromosome == 'chrI': seq = genome[0][start:end]
        if chromosome == 'chrII': seq = genome[1][start:end]
        if chromosome == 'chrIII': seq = genome[2][start:end]
        if chromosome == 'chrIV': seq = genome[3][start:end]
        if chromosome == 'chrV': seq = genome[4][start:end]
        if chromosome == 'chrX': seq = genome[5][start:end]
        seqs.append(seq)

## Write output
name = gene + "_" + stage + "_sequences.txt"
f_out = open(chip_path + "/" + name, "w")
for i in range(len(seqs)):
    f_out.write(headers[i] + "\n")
    f_out.write(seqs[i] + "\n")
f_out.close()

print ("Done!")
print ("\nThe output file is", name, "\n")
