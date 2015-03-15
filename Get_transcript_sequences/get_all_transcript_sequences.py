#!/usr/bin/env python

# dealing with command line input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-gtf", type=str, required=True)
parser.add_argument("-genome", type=str, required=True)
args = parser.parse_args()
gtf = args.gtf
genome = args.genome

# creating bed file
input = open(gtf)
output = open('all_exons.bed', 'w')
input.readline()
for line in input:
    split = line.split('\t')
    scaffold = split[0]
    start_co_ordinate = int(split[3]) - 1
    stop_co_ordinate = split[4]
    info = split[8]
    output.write(scaffold + '\t' + str(start_co_ordinate) + '\t' + stop_co_ordinate + '\t' + info + '\n')
input.close()
output.close()

# obtaining fasta sequences using bedtools
import subprocess
subprocess.call(["bedtools", "getfasta", "-fi", genome, "-bed", "all_exons.bed", "-fo", "all_exons.fa", "-name"])

# stitching together exon sequences
output = open('all_transcripts.fa', 'w')
input = open('all_exons.fa')
first_line = input.readline()
info = first_line.split('"')
transcript = info[3]
output.write(info[0] + info[1] + info[2] + info[3])
if '; gene name ' in info:
    output.write(info[6] + info [7])
output.write('\n')
transcript_previous = transcript
sequence = ''
for line in input:
    if line.startswith('>'):
        info = line.split('"')
        transcript = info[3]
        if transcript != transcript_previous:
            sequence = list(sequence)
            for i in range(0, len(sequence)):
                if i%60 == 0 and i != 0:
                    output.write('\n')
                output.write(sequence[i])
            output.write('\n' + info[0] + info[1] + info[2] + info[3])
            if '; gene_name ' in info:
                output.write(info[6] + info [7])
            output.write('\n')
            transcript_previous = transcript
            sequence = ''
    else:
        line_sequence = line.rstrip()
        sequence = sequence + line_sequence

sequence = list(sequence)
for i in range(0, len(sequence)):
    if i%60 == 0 and i != 0:
        output.write('\n')
    output.write(sequence[i])
