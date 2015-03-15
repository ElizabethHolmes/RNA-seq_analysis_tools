#!/usr/bin/env python

# dealing with command line input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-genes", type=str)
parser.add_argument("-gtf", type=str)
args = parser.parse_args()
genes_filename = args.genes
gtf_filename = args.gtf

# making list of genes from genes file
file = open(genes_filename)
genes = []
for line in file:
    gene = line.rstrip()
    genes.append(gene)
file.close()

# extracting co-ordinates from gtf file
delimiter = '\t'
start_co_ordinates = []
stop_co_ordinates = []
scaffolds = []
for gene in genes:
    file = open(gtf_filename)
    co_ordinates = []
    for line in file:
        split = line.split(delimiter)
        info = split[8]
        info = info.split(';')
        gene_id = info[0]
        gene_id = gene_id.split('"')
        gene_id = gene_id[1]
        if gene_id == gene:
            start_co_ordinate = int(split[3])
            stop_co_ordinate = int(split[4])
            co_ordinates.append(start_co_ordinate)
            co_ordinates.append(stop_co_ordinate)
            scaffold = split[0]
    start_co_ordinates.append(min(co_ordinates)-1)
    stop_co_ordinates.append(max(co_ordinates))
    scaffolds.append(scaffold)
    file.close()

# writing output
output = open('gene_co-ordinates.bed', 'w')
i = 0
for gene in genes:
    output.write(scaffolds[i] + '\t' + str(start_co_ordinates[i]) + '\t' + str(stop_co_ordinates[i]) + '\t' + gene + '\n')
    i = i + 1
output.close()
