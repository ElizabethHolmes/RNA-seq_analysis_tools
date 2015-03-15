#!/usr/bin/env python

# dealing with command line input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t",type=int)
parser.add_argument("-above_threshold", nargs="+", type=str)
parser.add_argument("-absent", nargs="+", type=str)
parser.add_argument("-other", nargs="+", type=str)
parser.add_argument("filename", type=str)
args = parser.parse_args()
threshold = args.t
above_threshold = args.above_threshold
absent = args.absent
other = args.other
filename = args.filename

# making lists of isoforms with expression above threshold in above_threshold samples
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
        if (sample1 == sample and sample1_expression == 0) or (sample2 == sample and sample2_expression == 0):
            if isoform not in isoform_list:
                isoform_list.append(isoform)
                gene_list.append(gene)
                locus_list.append(locus)
                annotation_ID_list.append(annotation_ID)
    file.close()
    initial_isoforms.append(isoform_list)

# extracting isoforms in all lists
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
    count = 0
    for list in initial_isoforms:
        if isoform in list:
            count = count + 1
    if count == sample_number:
            final_isoforms.append(isoform)
            final_genes.append(initial_gene_list[initial_isoform_list.index(isoform)])
            final_annotation_IDs.append(initial_annotation_IDs_list[initial_isoform_list.index(isoform)])
            final_loci.append(initial_loci_list[initial_isoform_list.index(isoform)])

# writing file headers
output = open(‘transcript_expression_profiles.txt', 'w')
output.write('Locus\tAnnotation gene ID\tCuffdiff gene ID\tCuffdiff transcript ID\t')
for sample in above_threshold:
    output.write(sample + '\t')
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
