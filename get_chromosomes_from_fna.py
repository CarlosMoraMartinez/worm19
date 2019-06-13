##### get_chromosomes_from_fna.py

"""
This script is used for split multiFASTA (.fna) files into multiple
FASTA (.fa) files, each of them with one of the chromosomes.

Autor: Erick Sousa
Created: 04/03/2019
"""

## Load modules
import os
path = os.getcwd()

## Create output dir:
genome = 'WBcel235'
out_path = path + "/" + genome + "/" + "FASTA_chr"
if not (os.path.exists(out_path)):
    try:
        os.mkdir(out_path)
    except OSError:
        print ("Creation of the directory %s failed" % out_path)
    else:
        print ("Successfully created the directory %s " % out_path)

## Get sequences from a .bed file:
ruta = path + "/" + genome + "/" + "GCA_000002985.3_WBcel235_genomic.fna"
f_in = open(ruta, "r")
fichero = f_in.read()
f_in.close()
lineas = fichero.splitlines(False)
pos_chr = []
## Get starting position from each chromosome
for i in range(len(lineas)):
    if lineas[i][0] == '>':
        pos_chr.append(i)

## For each chromosome:
for i in range(len(pos_chr)):
    # Get sequence
    if (i != (len(pos_chr) - 1)):
        sequence = lineas[pos_chr[i]:pos_chr[i+1]]
    else:
        sequence = lineas[pos_chr[i]:len(lineas)]
    # Get chromosome name
    header = lineas[pos_chr[i]].split()
    name = "chr_" + header[len(header) - 1] + ".fa"
    f_out = open((out_path + "/" + name), "w")
    # Write output file
    for j in sequence:
        f_out.write(j + "\n")
    print ("%s correctly saved" % name)
    f_out.close()
