# Differential_splicing
## About
Differential_splicing is a Python script designed to analyse differential splicing between samples using RNA-seq data [(Trapnell *et al.*, 2015)](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html). It works on the output of the Tuxedo suite of tools for RNA-seq analysis to identify introns in spliceforms that are expressed above a user-specified threshold in one or more user-specified samples and at or below a user-specified threshold in one or more other user-specified samples.

For example, Differential_splicing could be used to identify putative male-specific introns, which are present in spliceforms expressed in males, but absent from spliceforms of the same genes expressed in females. It may be advisable to use a threshold higher than 0 FPKM to detect genuine presence of spliceforms, to reduce false positives, and to use a threshold higher than 0 FPKM to confirm genuine absence of spliceforms, to reduce false negatives. False positives or false negatives could be caused by e.g. errors in transcript assembly or sample cross-contamination. So in this case the user might specify a threshold of 10 FPKM for male samples (introns must be in transcripts with an expression level above 10 FPKM in these samples), and 0.5 FPKM for female samples (introns must not be in transcripts with an expression level above 0.5 FPKM in these samples).

** PLEASE NOTE: Differential_splicing is one of my first scripts written as a beginner programmer and so I apologise if the code is inelegant, unconventional or otherwise sub-optimal; it works for the intended purpose and I provide it in case it might be useful to others, but with no guarantees. ** 

## Citation
Differential_splicing is not yet associated with a paper; to cite it please use:

    Sutton, ER. (2015). Differential_splicing [Software]. 
    Available at https://github.com/ElizabethSutton/Differential-splicing.

## Requirements
Differential_splicing requires Python. It has been tested only with Python 2.7.3.

Differential_splicing takes a few hours to run on a standard-sized dataset and so is probably best run on a high performance server with multiple cores.

Differential_splicing works on the output of the Tuxedo suite of tools for RNA-seq analysis. Files output by Cuffdiff, Cuffmerge and TopHat are required. It has been tested only with Cuffdiff v2.1.1, Cuffmerge v1.0.0 and TopHat v2.0.9.  

## Usage
Differential_splicing has a number of required and optional arguments, detailed below.

### *Required arguments*
* `above_threshold` - space-separated list of samples for which introns output must be in spliceforms with an FPKM value above `t` (see 'Optional arguments' section); the sample names must be written exactly as in the isoform_exp.diff file (see below)
* `absent` - space-separated list of samples for which introns output must not be in spliceforms with an FPKM value above `t_absent` (see 'Optional arguments' section); the sample names must be written exactly as in the isoform_exp.diff file (see below)
* `genome` - path to genome FASTA file
* `expression` - path to isoform_exp.diff file output by Cuffdiff
* `gtf` - path to merged.gtf file output by Cuffmerge
* `junctions` - space-separated list of paths to junctions.bed files output by Tophat for `above_threshold` and `other` samples, with the paths for the `above_threshold` samples listed first and in the same order as for the `above_threshold` argument, and the paths for the `other` samples listed second and in the same order as for the `other` argument

### *Optional arguments*
* `t` - threshold FPKM value - introns output must be in spliceforms with an FPKM value above this value in samples listed in the `above_threshold` argument; default value is 0
* `t_absent` - threshold absent FPKM value; introns output must not be in transcripts with an FPKM value above this in samples listed in the `absent` argument, and there must be an alternative transcript for the gene in question with an FPKM value above this in at least one of the samples listed in the `absent` argument (otherwise the program would identify introns in genes with sample-specific expression as well as genes with sample-specific splicing); default value is 0 
* `other` - space-separated list of samples for which expression information is desired but for which expression values do not influence the introns output; the sample names must be written exactly as in the isoform_exp.diff file
* `min_specific_intron` - introns output must have a length in nucleotides above this value; default value is 10

### *Example*
To identify putative male-specific introns as per the example in the 'About' section, the command might be:

    ./differential_splicing.py -t 10 -t_absent 0.5 -above_threshold male1 male2 -absent female1 female2 -genome genome.fa -expression ..isoform_exp.diff -gtf merged.gtf -junctions TopHat_files/Male1/junctions.bed TopHat_files/Male2/junctions.bed TopHatfiles/Female1/junctions.bed TopHatfiles/Female2/junctions.bed 

## Output

