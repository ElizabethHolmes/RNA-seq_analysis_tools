#!/usr/bin/env python

# dealing with command line input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t",type=float, default=0)
parser.add_argument("-above_threshold", nargs="+", type=str, required=True)
parser.add_argument("-absent", nargs="+", type=str, required=True)
parser.add_argument("-other", nargs="+", type=str)
parser.add_argument("-gtf", type=str, required=True)
parser.add_argument("-genome", type=str, required=True)
parser.add_argument("-min_specific_intron", type=int, default=10)
parser.add_argument("-t_absent", type=float, default=0)
parser.add_argument("-junctions", nargs="+", type=str)
parser.add_argument("-expression", type=str, required=True)
args = parser.parse_args()
threshold = args.t
above_threshold = args.above_threshold
absent = args.absent
other = args.other
gtf = args.gtf
genome = args.genome
min_specific_intron = args.min_specific_intron
absent_threshold = args.t_absent
junctions = args.junctions
expression = args.expression

# making directories
import os
import shutil
if os.path.exists('Sequences'):
	shutil.rmtree('Sequences')
os.mkdir('Sequences')

# making lists of isoforms with expression above threshold in above_threshold samples
filename = expression
initial_isoforms = []
initial_genes = []
initial_annotation_IDs = []
initial_loci = []
delimiter = '\t'
for sample in above_threshold:
    isoform_list = []
    gene_list = []
    locus_list = []
    annotation_ID_list = []
    file = open(filename)
    file.readline()
    for line in file:
        split = line.split(delimiter)
        isoform = split[0]
        gene = split[1]
        annotation_ID = split[2]
        locus = split[3]
        sample1 = split[4]
        sample2 = split[5]
        sample1_expression = split[7]
        sample1_expression = float(sample1_expression)
        sample2_expression = split[8]
        sample2_expression = float(sample2_expression)
        if (sample1 == sample and sample1_expression > threshold) or (sample2 == sample and sample2_expression > threshold):
            if isoform not in isoform_list:
                isoform_list.append(isoform)
                gene_list.append(gene)
                locus_list.append(locus)
                annotation_ID_list.append(annotation_ID)
    file.close()
    initial_isoforms.append(isoform_list)
    initial_genes.append(gene_list)
    initial_annotation_IDs.append(annotation_ID_list)
    initial_loci.append(locus_list)

# making lists of isoforms with no expression in absent samples
for sample in absent:
    isoform_list = []
    gene_list = []
    locus_list = []
    annotation_ID_list = []
    file = open(filename)
    file.readline()
    for line in file:
        split = line.split(delimiter)
        isoform = split[0]
        gene = split[1]
        annotation_ID = split[2]
        locus = split[3]
        sample1 = split[4]
        sample2 = split[5]
        sample1_expression = split[7]
        sample1_expression = float(sample1_expression)
        sample2_expression = split[8]
        sample2_expression = float(sample2_expression)
        if (sample1 == sample and sample1_expression <= absent_threshold) or (sample2 == sample and sample2_expression <= absent_threshold):
            if isoform not in isoform_list:
                isoform_list.append(isoform)
                gene_list.append(gene)
                locus_list.append(locus)
                annotation_ID_list.append(annotation_ID)
    file.close()
    initial_isoforms.append(isoform_list)

# extracting isoforms in all lists and for which the whole gene is not specific
final_isoforms = []
final_genes = []
final_loci = []
final_annotation_IDs = []
sample_number = len(initial_isoforms)
initial_isoform_list = initial_isoforms[0]
initial_gene_list = initial_genes[0]
initial_loci_list = initial_loci[0]
initial_annotation_IDs_list = initial_annotation_IDs[0]
for isoform in initial_isoform_list:
    initial_gene = initial_gene_list[initial_isoform_list.index(isoform)]
    count = 0
    for list in initial_isoforms:
        if isoform in list:
            count = count + 1
    if count == sample_number:
        file = open(filename)
        file.readline()
        for line in file:
            split = line.split(delimiter)
            gene = split[1]
            sample1 = split[4]
            sample2 = split[5]
            sample1_expression = split[7]
            sample1_expression = float(sample1_expression)
            sample2_expression = split[8]
            sample2_expression = float(sample2_expression)
            if gene == initial_gene and ((sample1 in absent and sample1_expression > absent_threshold) or (sample2 in absent and sample2_expression > absent_threshold)):
                final_isoforms.append(isoform)
                final_genes.append(initial_gene_list[initial_isoform_list.index(isoform)])
                final_annotation_IDs.append(initial_annotation_IDs_list[initial_isoform_list.index(isoform)])
                final_loci.append(initial_loci_list[initial_isoform_list.index(isoform)])
                break
        file.close()

# writing file headers
output = open('isoform_all_specific_expression_profiles.txt', 'w')
output.write('Locus\tAnnotation gene ID\tCuffdiff gene ID\tCuffdiff transcript ID\t')
for sample in above_threshold:
    output.write(sample + '\t')
if other:
	for sample in other:
		output.write(sample + '\t')
output.write('\n')

# writing expression values
for final_isoform in final_isoforms:
    gene = final_genes[final_isoforms.index(final_isoform)]
    annotation_ID = final_annotation_IDs[final_isoforms.index(final_isoform)]
    locus = final_loci[final_isoforms.index(final_isoform)]
    output.write(locus + '\t' + annotation_ID + '\t' + gene + '\t' + final_isoform + '\t')
    for sample in above_threshold:
        file = open(filename)
        file.readline()
        for line in file:
            split = line.split(delimiter)
            isoform = split[0]
            sample1 = split[4]
            sample2 = split[5]
            sample1_expression = split[7]
            sample2_expression = split[8]
            if sample1 == sample and isoform == final_isoform:
                output.write(sample1_expression + '\t')
                break
            elif sample2 == sample and isoform == final_isoform:
                output.write(sample2_expression + '\t')
                break
        file.close()
    if other:
    	for sample in other:
        	file = open(filename)
        	file.readline()
        	for line in file:
            		split = line.split(delimiter)
            		isoform = split[0]
            		sample1 = split[4]
            		sample2 = split[5]
            		sample1_expression = split[7]
            		sample2_expression = split[8]
            	if sample1 == sample and isoform == final_isoform:
                	output.write(sample1_expression + '\t')
                	break
            	elif sample2 == sample and isoform == final_isoform:
                	output.write(sample2_expression + '\t')
                	break
        	file.close()
    output.write('\n')
output.close()

# writing list of genes with no duplicates
final_genes2 = []
for gene in final_genes:
    if gene not in final_genes2:
        final_genes2.append(gene)

# getting gene co_ordinates
gene_start_co_ordinates = dict.fromkeys(final_genes2)
gene_stop_co_ordinates = dict.fromkeys(final_genes2)
filename = gtf
for gene in final_genes2:
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

# creating isoform number sequences for comparison to determine exon segments missing in above_threshold samples
# getting isoform co_ordinates for bed files
genes_and_transcripts = dict.fromkeys(final_genes2)
gene_scaffolds_for_bed = dict.fromkeys(final_genes2)

all_intron_sequences_list = []
exon_co_ordinates_for_bed = []
transcripts_for_each_gene = []
intron_sequences_for_each_transcript = []
intron_sequence = ''
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
    if gene_id in final_genes2:
        scaffold = split[0]
        transcript_id = info[1]
        transcript_id = transcript_id.split('"')
        transcript_id = transcript_id[1]
        exon_start_co_ordinate = int(split[3])
        gene_scaffolds_for_bed[gene_id] = scaffold
        if transcript_id != transcript_id_previous: # i.e. first time transcript is encountered
            if transcript_id_previous != 0:
                all_intron_sequences_list.append(intron_sequences_for_each_transcript)
                exon_co_ordinates_for_bed.append(exon_co_ordinates_for_each_transcript)
                intron_sequences_for_each_transcript = []
                exon_co_ordinates_for_each_transcript = []
                transcripts_for_each_gene.append(transcript_id_previous)
            transcript_id_previous = transcript_id
            exon_stop_co_ordinate = int(split[4])
            exon_co_ordinates_for_each_transcript.append(exon_start_co_ordinate)
            exon_co_ordinates_for_each_transcript.append(exon_stop_co_ordinate)
            all_transcripts.append(transcript_id)
        else: # i.e. transcript has been encountered before
            for i in range(exon_stop_co_ordinate+1, exon_start_co_ordinate):
                intron_sequence = intron_sequence + str(i) + '.'
            intron_sequences_for_each_transcript.append(intron_sequence)
            intron_sequence = ''
            exon_stop_co_ordinate = int(split[4])
            exon_co_ordinates_for_each_transcript.append(exon_start_co_ordinate)
            exon_co_ordinates_for_each_transcript.append(exon_stop_co_ordinate)
        if gene_id != gene_id_previous:
            if gene_id_previous !=0:
                genes_and_transcripts[gene_id_previous] = transcripts_for_each_gene
                transcripts_for_each_gene = []
            gene_id_previous = gene_id

if transcripts_for_each_gene:
    all_intron_sequences_list.append(intron_sequences_for_each_transcript)
    exon_co_ordinates_for_bed.append(exon_co_ordinates_for_each_transcript)
    if transcript_id not in transcripts_for_each_gene:
        transcripts_for_each_gene.append(transcript_id)
    genes_and_transcripts[gene_id_previous] = transcripts_for_each_gene

# comparing isoform number sequences to identify transcripts with exon segments missing in above_threshold samples
# making a list of all isoforms for genes with specific intron isoforms
intron_co_ordinates_for_bed = []
genes_and_isoforms_for_bed = []
isoform_scaffolds_for_bed = []
specific_intron_isoforms = []
specific_intron_isoforms_with_replicates = [] # for isoforms that have more than one specific intron
specific_intron_genes_with_replicates = [] # for isoforms that have more than one specific intron
specific_intron_genes = []
all_isoforms_list = []
for isoform in final_isoforms:
    gene = final_genes[final_isoforms.index(isoform)]
    intron_sequences = all_intron_sequences_list[all_transcripts.index(isoform)]
    all_isoforms = genes_and_transcripts[gene]
    other_intron_sequences = []
    for other_isoform in all_isoforms:
        if other_isoform not in final_isoforms:
            final_list = []
            for list in all_intron_sequences_list[all_transcripts.index(other_isoform)]:
                final_list.append(list)
            other_intron_sequences.append(final_list)
    for sequence in intron_sequences:
        length = sequence.count('.')
        if length >= min_specific_intron:
            count = 0
            for other_intron_sequence in other_intron_sequences:
                if sequence not in other_intron_sequence:
                    count = count + 1
            if count == len(other_intron_sequences):
                sequence_split = sequence.split('.')
                intron_co_ordinates_for_bed.append(sequence_split[0])
                intron_co_ordinates_for_bed.append(sequence_split[len(sequence_split)-2])
                genes_and_isoforms_for_bed.append(gene + ', ' + isoform)
                isoform_scaffolds_for_bed.append(gene_scaffolds_for_bed[gene])
                specific_intron_isoforms_with_replicates.append(isoform)
                specific_intron_genes_with_replicates.append(gene)
                if isoform not in specific_intron_isoforms:
                    specific_intron_isoforms.append(isoform)
                if gene not in specific_intron_genes:
                    specific_intron_genes.append(gene)
                    for entry in all_isoforms:
                        all_isoforms_list.append(entry)

# writing file headers
output = open('gene_expression_profiles.txt', 'w')
output.write('Locus\tAnnotation gene ID\tCuffdiff gene ID\tCuffdiff transcript ID\tSpecific?\t')
for sample in above_threshold:
    output.write(sample + '\t')
if other:
	for sample in other:
    		output.write(sample + '\t')
for sample in absent:
    output.write(sample + '\t')
output.write('\n')

# writing expression values
filename = expression
for isoform in all_isoforms_list:
    info_written = 'false'
    for sample in above_threshold:
        file = open(filename)
        file.readline()
        for line in file:
             split = line.split(delimiter)
             isoform_ID = split[0]
             gene = split[1]
             annotation_ID = split[2]
             locus = split[3]
             sample1 = split[4]
             sample2 = split[5]
             sample1_expression = split[7]
             sample2_expression = split[8]
             if isoform == isoform_ID:
                if info_written == 'false':
                    output.write(locus + '\t' + annotation_ID + '\t' + gene + '\t' + isoform + '\t')
                    info_written = 'true'
                    if isoform in specific_intron_isoforms:
                        output.write('Yes\t')
                    else:
                        output.write(' \t')
                if sample1 == sample:
                    output.write(sample1_expression + '\t')
                    break
                elif sample2 == sample:
                    output.write(sample2_expression + '\t')
                    break
        file.close()

    if other:
	for sample in other:
        	file = open(filename)
        	file.readline()
        	for line in file:
             		split = line.split(delimiter)
             		isoform_ID = split[0]
		        gene = split[1]
             		sample1 = split[4]
             		sample2 = split[5]
             		sample1_expression = split[7]
             		sample2_expression = split[8]
             		if isoform == isoform_ID:
                		if sample1 == sample:
                    			output.write(sample1_expression + '\t')
                    			break
                		elif sample2 == sample:
                    			output.write(sample2_expression + '\t')
                    			break
        	file.close()

    for sample in absent:
        file = open(filename)
        file.readline()
        for line in file:
             split = line.split(delimiter)
             isoform_ID = split[0]
             gene = split[1]
             sample1 = split[4]
             sample2 = split[5]
             sample1_expression = split[7]
             sample2_expression = split[8]
             if isoform == isoform_ID:
                if sample1 == sample:
                    output.write(sample1_expression + '\t')
                    break
                elif sample2 == sample:
                    output.write(sample2_expression + '\t')
                    break
        file.close()
    output.write('\n')
output.close()

# creating bed files
for gene in specific_intron_genes:
    output = open('Sequences/' + gene + '.bed', 'w')
    scaffold = gene_scaffolds_for_bed[gene]
    start_co_ordinate = gene_start_co_ordinates[gene]
    stop_co_ordinate = gene_stop_co_ordinates[gene]
    output.write(scaffold + '\t' + str(start_co_ordinate) + '\t' + str(stop_co_ordinate) + '\t' + gene + ' genomic sequence' + '\n')
    for transcript in genes_and_transcripts[gene]:
        if transcript in specific_intron_isoforms:
            for i in range(0, len(exon_co_ordinates_for_bed[all_transcripts.index(transcript)])/2):
                output.write(scaffold + '\t')
                output.write(str(exon_co_ordinates_for_bed[all_transcripts.index(transcript)][2*i]-1) + '\t')
                output.write(str(exon_co_ordinates_for_bed[all_transcripts.index(transcript)][2*i + 1]) + '\t')
                output.write(transcript + '-specific_intron_isoform' + '\n')
    for transcript in genes_and_transcripts[gene]:
        if transcript not in specific_intron_isoforms:
            for i in range(0, len(exon_co_ordinates_for_bed[all_transcripts.index(transcript)])/2):
                output.write(scaffold + '\t')
                output.write(str(exon_co_ordinates_for_bed[all_transcripts.index(transcript)][2*i]-1) + '\t')
                output.write(str(exon_co_ordinates_for_bed[all_transcripts.index(transcript)][2*i + 1]) + '\t')
                output.write(transcript + '\n')
    output.close()

# obtaining fasta sequences using bedtools
import subprocess
for gene in specific_intron_genes:
    filename = 'Sequences/' + gene + '.bed'
    subprocess.call(["bedtools", "getfasta", "-fi", genome, "-bed", filename, "-fo", "Sequences/" + gene + "_exons.fa", "-name"])

# writing fasta sequences with gaps to files
output2 = open('Sequences/all_genes_genomic_sequences.fa', 'w')
all_sequences = dict.fromkeys(specific_intron_isoforms)
for gene in specific_intron_genes:
    input = open('Sequences/' + gene + '_exons.fa')
    output = open('Sequences/' + gene + '.fa', 'w')
    header = input.readline()
    output.write(header)
    output2.write(header)
    sequence = ''
    transcript_previous = ''
    for line in input:
        if line.startswith('>'):
            for i in range(0, len(sequence)):
                 if i%60 == 0 and i != 0:
                    output.write('\n')
                    output2.write('\n')
                 output.write(sequence[i])
                 output2.write(sequence[i])
            output.write('\n' + line)
            output2.write('\n')
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
                 if transcript_previous in specific_intron_isoforms:
                     all_sequences[transcript_previous] = sequence
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
    if transcript in specific_intron_isoforms:
        all_sequences[transcript] = sequence
output.close()
output2.close()

# writing file headers
output = open('specific_intron_information.txt', 'w')
output.write('Locus\tAnnotation gene ID\tCuffdiff gene ID\tCuffdiff transcript ID\t')
for sample in above_threshold:
    output.write(sample + ' expression\t')
if other:
	for sample in other:
    		output.write(sample + ' expression\t')
for sample in above_threshold:
    output.write(sample + ' junction coverage\t')
if other:
	for sample in other:
    		output.write(sample + ' junction coverage\t')
output.write('Intron start in gene\tIntron stop in gene\tExon junction in flanking sequence\tOther exon junctions in flanking sequence\tFlanking sequence')
output.write('\n')

# writing expression values, junction coverage values, intron co-ordinates, exon junction co-ordinates, and flanking sequence
i = 0
for isoform in specific_intron_isoforms_with_replicates:
    # expression values
    input = open('isoform_all_specific_expression_profiles.txt')
    input.readline()
    for line in input:
        if isoform in line:
            expression_values = line.rstrip()
            output.write(expression_values + '\t')
    input.close()
    # junction coverage values
    for file in junctions:
        input = open(file)
        input.readline()
        for line in input:
            split = line.split(delimiter)
            scaffold = split[0]
            block_sizes = split[10]
            block_sizes = block_sizes.split(',')
            block1_size = block_sizes[0]
            block2_size = block_sizes[1]
            junction_start_co_ordinate = int(split[1]) + int(block1_size)
            junction_stop_co_ordinate = int(split[2]) - int(block2_size)
            if scaffold == isoform_scaffolds_for_bed[i] and junction_start_co_ordinate == int(intron_co_ordinates_for_bed[2*i])-1 and junction_stop_co_ordinate == int(intron_co_ordinates_for_bed[2*i+1]):
                output.write(split[4])
                break
        output.write('\t')
        input.close()
    # intron co-ordinates
    gene = specific_intron_genes_with_replicates[i]
    intron_start_co_ordinate = int(intron_co_ordinates_for_bed[2*i]) - gene_start_co_ordinates[gene]
    intron_stop_co_ordinate = int(intron_co_ordinates_for_bed[2*i+1]) - gene_start_co_ordinates[gene]
    output.write(str(intron_start_co_ordinate)+'\t'+str(intron_stop_co_ordinate)+'\t')
    # exon junction co-ordinates and flanking sequence
    sequence = all_sequences[isoform]
    count_ungapped = 0
    count_gapped = 0
    sequence_ungapped = ''
    base_previous = ''
    other_exon_junction_co_ordinates = []
    junction_co_ordinate = 0
    for base in sequence:
        count_gapped = count_gapped + 1
        if base != "-":
            count_ungapped = count_ungapped + 1
            sequence_ungapped = sequence_ungapped + base
            if count_gapped == intron_stop_co_ordinate+1:
                if count_ungapped > 500:
                    start_co_ordinate = count_ungapped - 500
                    junction_co_ordinate = 501
                else:
                    start_co_ordinate = 1
                    junction_co_ordinate = count_ungapped
            elif base_previous == "-" and count_ungapped !=1:
                other_exon_junction_co_ordinates.append(count_ungapped)
            if count_ungapped == start_co_ordinate + junction_co_ordinate + 499 and count_gapped > intron_stop_co_ordinate:
                break
        base_previous=base
    output.write(str(junction_co_ordinate)+'\t')
    if other_exon_junction_co_ordinates:
        if len(other_exon_junction_co_ordinates)>1:
            for j in range (0,len(other_exon_junction_co_ordinates)-1):
                if other_exon_junction_co_ordinates[j] > start_co_ordinate-1:
                    output.write(str(other_exon_junction_co_ordinates[j]-start_co_ordinate)+', ')
        if other_exon_junction_co_ordinates[len(other_exon_junction_co_ordinates)-1] > start_co_ordinate-1:
            output.write(str(other_exon_junction_co_ordinates[len(other_exon_junction_co_ordinates)-1]-start_co_ordinate)+'\t')
        else:
            output.write('\t')
    else:
        output.write('\t')
    for k in range(start_co_ordinate-1, len(sequence_ungapped)):
        output.write(sequence_ungapped[k])
    output.write('\n')
    i = i + 1
output.close()

# writing intron co-ordinates bed file
output = open('intron_co-ordinates.bed', 'w')
for i in range(0, len(genes_and_isoforms_for_bed)):
    output.write(isoform_scaffolds_for_bed[i] + '\t' + intron_co_ordinates_for_bed[2*i] + '\t' + intron_co_ordinates_for_bed[2*i+1] + '\t' + genes_and_isoforms_for_bed[i] + '\n')
output.close()

# removing files
for gene in specific_intron_genes:
    filename = 'Sequences/' + gene + '.bed'
    os.remove(filename)
    filename = 'Sequences/' + gene + '_exons.fa'
    os.remove(filename)
os.remove('isoform_all_specific_expression_profiles.txt')
