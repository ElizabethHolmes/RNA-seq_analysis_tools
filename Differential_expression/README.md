# Differential_expression
## About
Differential_expression is a combination of two Python scripts for analysing differential expression between samples using RNA-seq data. These scripts work on the output of the Tuxedo suite of tools for RNA-seq analysis [(Trapnell *et al.*, 2012)](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html) to identify genes or transcripts that are expressed above a user-specified threshold in one or more user-specified samples and at or below a user-specified threshold in one or more other user-specified samples. Genes are identified with `differential_expression_genes.py` and transcripts are identified with `differential_expression_transcripts.py`.

For example, Differential_expression could be used to identify putative male-specifically expressed genes, which are expressed in males, but not in females. It may be advisable to use a threshold higher than 0 FPKM to detect genuine expression of genes, to reduce false positives, and to use a threshold higher than 0 FPKM to confirm genuine absence of expression of genes, to reduce false negatives. False positives or false negatives could be caused by e.g. errors in transcript assembly or sample cross-contamination. So in this case the user might specify a threshold of 10 FPKM for male samples (genes must have an expression level above 10 FPKM in these samples), and 0.5 FPKM for female samples (genes must not an expression level above 0.5 FPKM in these samples).

The information output by Differential_expression is already available in the output of Cuffdiff, but not in a convenient format; Differential_expression extracts only the genes of interest and lists information on these genes in a convenient tabular format.

** PLEASE NOTE: The Differential_expression scripts are some of my first scripts written as a beginner programmer and so I apologise if the code is inelegant, unconventional or otherwise sub-optimal; they work for the intended purpose and I provide them in case they might be useful to others, but with no guarantees. ** 

## Citation
Differential_expression is not yet associated with a paper; to cite it please use:

    Sutton, ER. (2015). Differential_expression [Software]. 
    Available at https://github.com/ElizabethSutton/RNA-seq_differential_analysis_tools/Differential_expression.

## Requirements
Differential_expression requires Python. It has been tested only with Python 2.7.3.

Differential_expression works on the output of the Tuxedo suite of tools for RNA-seq analysis. A file output by Cuffdiff is required. It has been tested only with Cuffdiff v2.1.1.  

## Usage
Both Differential_expression scripts have a number of required and optional arguments, detailed below.

### *Required arguments*
* `above_threshold` - space-separated list of samples for which genes/transcripts output must have an FPKM value above `t` (see 'Optional arguments' section); the sample names must be written exactly as in the isoform_exp.diff file (see below)
* `absent` - space-separated list of samples for which genes/transcripts output must not have an FPKM value above `t_absent` (see 'Optional arguments' section); the sample names must be written exactly as in the isoform_exp.diff file (see below)
* `expression` - path to isoform_exp.diff file output by Cuffdiff

### *Optional arguments*
* `t` - threshold FPKM value - genes/transcripts output must have an FPKM value above this in samples listed in the `above_threshold` argument; default value is 0
* `t_absent` - threshold absent FPKM value; genes/transcripts output must not have an FPKM value above this in samples listed in the `absent` argument; default value is 0 
* `other` - space-separated list of samples for which expression information is desired but for which expression values do not influence the genes/transcript output; the sample names must be written exactly as in the isoform_exp.diff file

### *Example*
To identify putative male-specifically expressed genes as per the example in the 'About' section, the command might be:

    ./differential_expression_genes.py -t 10 -t_absent 0.5 -above_threshold male1 male2 -absent female1 female2 -expression isoform_exp.diff 

## Output
Both Differential_expression scripts produce one output file, called `gene_expression_profiles.txt` for `differential_expression_genes.py`, and `transcript_expression_profiles.txt` for `differential_expression_transcripts.py`.
This file lists the genes/transcripts identified as fitting the user-specified criteria (e.g. male-specifically expressed) and their expression values in the user-specified samples in tabular format.
