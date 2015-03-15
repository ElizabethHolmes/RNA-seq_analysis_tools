#!/usr/bin/env python

# dealing with command line input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t",type=int, default=0)
parser.add_argument("-above_threshold", nargs="+", type=str)
parser.add_argument("-absent", nargs="+", type=str)
parser.add_argument("-other", nargs="+", type=str)
parser.add_argument("-t_absent", type=int, default=0)
parser.add_argument("-expression", type=str)
args = parser.parse_args()
threshold = args.t
above_threshold = args.above_threshold
absent = args.absent
other = args.other
absent_threshold = args.t_absent
expression = args.expression

# making lists of genes with expression above threshold in above_threshold samples
initial_genes = []
initial_annotation_IDs = []
initial_loci = []
delimiter = '\t'
for sample in above_threshold:
    gene_list = []
    locus_list = []
    annotation_ID_list = []
    file = open(expression)
    file.readline()
    for line in file:
        split = line.split(delimiter)
        gene = split[0]
        annotation_ID = split[2]
        locus = split[3]
        sample1 = split[4]
        sample2 = split[5]
        sample1_expression = split[7]
        sample1_expression = float(sample1_expression)
        sample2_expression = split[8]
        sample2_expression = float(sample2_expression)
        if (sample1 == sample and sample1_expression > threshold) or (sample2 == sample and sample2_expression > threshold):
            if gene not in gene_list:
                gene_list.append(gene)
                locus_list.append(locus)
                annotation_ID_list.append(annotation_ID)
    file.close()
    initial_genes.append(gene_list)
    initial_annotation_IDs.append(annotation_ID_list)
    initial_loci.append(locus_list)

# making lists of genes with no expression in absent samples
for sample in absent:
    gene_list = []
    file = open(expression)
    file.readline()
    for line in file:
        split = line.rstrip().split(delimiter)
        gene = split[0]
        annotation_ID = split[2]
        locus = split[3]
        sample1 = split[4]
        sample2 = split[5]
        sample1_expression = split[7]
        sample1_expression = float(sample1_expression)
        sample2_expression = split[8]
        sample2_expression = float(sample2_expression)
        if (sample1 == sample and sample1_expression <= absent_threshold) or (sample2 == sample and sample2_expression <= absent_threshold):
            if gene not in gene_list:
                gene_list.append(gene)
                locus_list.append(locus)
                annotation_ID_list.append(annotation_ID)
    file.close()
    initial_genes.append(gene_list)

# extracting genes in all lists
final_genes = []
final_loci = []
final_annotation_IDs = []
sample_number = len(initial_genes)
initial_gene_list = initial_genes[0]
initial_loci_list = initial_loci[0]
initial_annotation_IDs_list = initial_annotation_IDs[0]
for gene in initial_gene_list:
    count = 0
    for list in initial_genes:
        if gene in list:
            count = count + 1
    if count == sample_number:
            final_genes.append(gene)
            final_annotation_IDs.append(initial_annotation_IDs_list[initial_gene_list.index(gene)])
            final_loci.append(initial_loci_list[initial_gene_list.index(gene)])

# writing file headers
output = open('gene_expression_profiles.txt', 'w')
output.write('Locus\tAnnotation gene ID\tCuffdiff gene ID\t')
for sample in above_threshold:
    output.write(sample + '\t')
if other:
	for sample in other:
    		output.write(sample + '\t')
for sample in absent:
    output.write(sample + '\t')
output.write('\n')

# writing expression values
for final_gene in final_genes:
    annotation_ID = final_annotation_IDs[final_genes.index(final_gene)]
    locus = final_loci[final_genes.index(final_gene)]
    output.write(locus + '\t' + annotation_ID + '\t' + final_gene + '\t')
    for sample in above_threshold:
        file = open(expression)
        file.readline()
        for line in file:
            split = line.split(delimiter)
            gene = split[0]
            sample1 = split[4]
            sample2 = split[5]
            sample1_expression = split[7]
            sample2_expression = split[8]
            if sample1 == sample and gene == final_gene:
                    output.write(sample1_expression + '\t')
                    break
            elif sample2 == sample and gene == final_gene:
                    output.write(sample2_expression + '\t')
                    break
        file.close()

    if other:
		for sample in other:
        		file = open(expression)
        		file.readline()
        		for line in file:
            			split = line.split(delimiter)
            			gene = split[0]
            			sample1 = split[4]
            			sample2 = split[5]
            			sample1_expression = split[7]
            			sample2_expression = split[8]
            			if sample1 == sample and gene == final_gene:
                			output.write(sample1_expression + '\t')
                			break
            			elif sample2 == sample and gene == final_gene:
                			output.write(sample2_expression + '\t')
                			break
        		file.close()

    for sample in absent:
	file = open(expression)
	file.readline()
	for line in file:
	    split = line.split(delimiter)
	    gene = split[0]
	    sample1 = split[4]
	    sample2 = split[5]
            sample1_expression = split[7]
  	    sample2_expression = split[8]
	    if sample1 == sample and gene == final_gene:
	 	output.write(sample1_expression + '\t')
		break
	    elif sample2 == sample and gene == final_gene:
		output.write(sample2_expression + '\t')
		break
    output.write('\n')
output.close()
