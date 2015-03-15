#!/usr/bin/env python

# dealing with command line input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-gtf", type=str, required=True)
parser.add_argument("-genome", type=str, required=True)
parser.add_argument("-genes", type=str, required=True)
args = parser.parse_args()
gtf = args.gtf
genome = args.genome
genes_file = args.genes

# creating list of genes
file = open(genes_file)
genes=[]
for line in file:
    gene=line.rstrip()
    genes.append(gene)
file.close()

# getting gene co_ordinates
gene_start_co_ordinates = dict.fromkeys(genes)
gene_stop_co_ordinates = dict.fromkeys(genes)
filename = gtf
delimiter = '\t'
for gene in genes:
    file = open(filename)
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
    start_co_ordinate = min(co_ordinates)-1
    stop_co_ordinate = max(co_ordinates)
    gene_start_co_ordinates[gene] = start_co_ordinate
    gene_stop_co_ordinates[gene] = stop_co_ordinate
    file.close()

# getting isoform co_ordinates for bed files
genes_and_transcripts = dict.fromkeys(genes)
gene_scaffolds_for_bed = dict.fromkeys(genes)

all_intron_sequences_list = []
exon_co_ordinates_for_bed = []
transcripts_for_each_gene = []
exon_co_ordinates_for_each_transcript = []
all_transcripts = []

transcript_id_previous = 0
gene_id_previous = 0

filename = gtf
file = open(filename)

for line in file:
    split = line.split(delimiter)
    info = split[8]
    info = info.split(';')
    gene_id = info[0]
    gene_id = gene_id.split('"')
    gene_id = gene_id[1]
    if gene_id in genes:
        scaffold = split[0]
        transcript_id = info[1]
        transcript_id = transcript_id.split('"')
        transcript_id = transcript_id[1]
        exon_start_co_ordinate = int(split[3])
        gene_scaffolds_for_bed[gene_id] = scaffold
        if transcript_id != transcript_id_previous: # i.e. first time transcript is encountered
            if transcript_id_previous != 0:
                exon_co_ordinates_for_bed.append(exon_co_ordinates_for_each_transcript)
                exon_co_ordinates_for_each_transcript = []
                transcripts_for_each_gene.append(transcript_id_previous)
            transcript_id_previous = transcript_id
            exon_stop_co_ordinate = int(split[4])
            exon_co_ordinates_for_each_transcript.append(exon_start_co_ordinate)
            exon_co_ordinates_for_each_transcript.append(exon_stop_co_ordinate)
            all_transcripts.append(transcript_id)
        else: # i.e. transcript has been encountered before
            exon_stop_co_ordinate = int(split[4])
            exon_co_ordinates_for_each_transcript.append(exon_start_co_ordinate)
            exon_co_ordinates_for_each_transcript.append(exon_stop_co_ordinate)
        if gene_id != gene_id_previous:
            if gene_id_previous !=0:
                genes_and_transcripts[gene_id_previous] = transcripts_for_each_gene
                transcripts_for_each_gene = []
            gene_id_previous = gene_id

exon_co_ordinates_for_bed.append(exon_co_ordinates_for_each_transcript)
transcripts_for_each_gene.append(transcript_id)
genes_and_transcripts[gene_id_previous] = transcripts_for_each_gene

# creating bed files
for gene in genes:
    output = open(gene + '.bed', 'w')
    scaffold = gene_scaffolds_for_bed[gene]
    start_co_ordinate = gene_start_co_ordinates[gene]
    stop_co_ordinate = gene_stop_co_ordinates[gene]
    output.write(scaffold + '\t' + str(start_co_ordinate) + '\t' + str(stop_co_ordinate) + '\t' + gene + ' genomic sequence' + '\n')
    for transcript in genes_and_transcripts[gene]:
            for i in range(0, len(exon_co_ordinates_for_bed[all_transcripts.index(transcript)])/2):
                output.write(scaffold + '\t')
                output.write(str(exon_co_ordinates_for_bed[all_transcripts.index(transcript)][2*i]-1) + '\t')
                output.write(str(exon_co_ordinates_for_bed[all_transcripts.index(transcript)][2*i + 1]) + '\t')
                output.write(transcript + '\n')
    output.close()

# obtaining fasta sequences using bedtools
import subprocess
for gene in genes:
    filename = gene + '.bed'
    subprocess.call(["bedtools", "getfasta", "-fi", genome, "-bed", filename, "-fo", gene + "_exons.fa", "-name"])

# writing fasta sequences with gaps to files
for gene in genes:
    input = open(gene + '_exons.fa')
    output = open(gene + '.fa', 'w')
    header = input.readline()
    output.write(header)
    sequence = ''
    transcript_previous = ''
    for line in input:
        if line.startswith('>'):
            for i in range(0, len(sequence)):
                 if i%60 == 0 and i != 0:
                    output.write('\n')
                 output.write(sequence[i])
            output.write('\n' + line)
            transcript = line.rstrip()
            transcript = transcript.split('>')
            transcript = transcript[1]
            transcript = transcript.split('-')
            transcript = transcript[0]
            transcript_previous = transcript
            sequence = ''
            break
        else:
            line_sequence = line.rstrip()
            sequence = sequence + line_sequence
        
    transcript_co_ordinates = exon_co_ordinates_for_bed[all_transcripts.index(transcript)]
    exon_count = 0
    
    for line in input:
        if line.startswith('>'):
            transcript = line.rstrip()
            transcript = transcript.split('>')
            transcript = transcript[1]
            transcript = transcript.split('-')
            transcript = transcript[0]
            if transcript == transcript_previous:
                 exon_count = exon_count + 1
                 start_co_ordinate = transcript_co_ordinates[2*exon_count-1] + 1
                 stop_co_ordinate = transcript_co_ordinates[2*exon_count]
                 for i in range(start_co_ordinate, stop_co_ordinate):
                    sequence = sequence + '-'
            else:
                 if transcript_co_ordinates[2*exon_count+1] != gene_stop_co_ordinates[gene]:
                    for i in range (transcript_co_ordinates[2*exon_count+1], gene_stop_co_ordinates[gene]):
                        sequence = sequence + '-'
                 for i in range(0, len(sequence)):
                    if i%60 == 0 and i != 0:
                        output.write('\n')
                    output.write(sequence[i])
                 output.write('\n' + line)
                 transcript_previous = transcript
                 transcript_co_ordinates = exon_co_ordinates_for_bed[all_transcripts.index(transcript)]
                 exon_count = 0
                 sequence= ''
        else:
            if exon_count == 0:
                if transcript_co_ordinates[0]-1 != gene_start_co_ordinates[gene]:
                    for i in range (gene_start_co_ordinates[gene], transcript_co_ordinates[0]-1):
                        sequence = sequence + '-'
            line_sequence = line.rstrip()
            sequence = sequence + line_sequence

    if transcript_co_ordinates[2*exon_count+1] != gene_stop_co_ordinates[gene]:
        for i in range (transcript_co_ordinates[2*exon_count+1], gene_stop_co_ordinates[gene]):
            sequence = sequence + '-'
    for i in range(0, len(sequence)):
        if i%60 == 0 and i != 0:
            output.write('\n')
        output.write(sequence[i])
output.close()

# removing files
import os
for gene in genes:
    filename = gene + '.bed'
    os.remove(filename)
    filename = gene + '_exons.fa'
    os.remove(filename)
