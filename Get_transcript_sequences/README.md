# Get_transcript_sequences
## About
Get_transcript_sequences is a pair of Python scripts for obtaining transcript sequences from the output of the Tuxedo suite of tools for RNA-seq analysis [(Trapnell *et al.*, 2012)](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html), though it could also ge used on gtf files from other sources. The `get_all_transcript_sequences.py` script obtains all predicted transcript sequences, which might be useful for e.g. creating a BLAST database to BLAST other sequences against. The `get_transcript_sequences_aligned.py` script obtains predicted transcript sequences in an aligned format for genes of interest, which might be useful for e.g. primer design to test genes of interest.

## Citation
To cite RNA-seq_analysis tools, please cite the following paper: [Sutton ER, Yu Y, Shimeld SM, White-Cooper H, Alphey L: Identification of genes for engineering the male germline of *Aedes aegypti* and *Ceratitis capitata*. BMC Genomics 2016 17:948.](https://www.ncbi.nlm.nih.gov/pubmed/27871244)

## Requirements
Get_transcript_sequences requires Python. It has been tested only with Python 2.7.3. It also requires [bedtools](https://github.com/arq5x/bedtools2). It has been tested only with bedtools version 2.16.2.

## Usage
Both Get_transcript_sequences scripts have two required arguments, and `get_transcript_sequences_aligned.py` has one additional argument, detailed below.

### *Required arguments*
* `genome` - path to genome FASTA file
* `gtf` - path to gtf file

Additional argument for `get_transcript_sequences_aligned.py`:
* `genes` - path to text file containing list of genes (in a single column) for which transcript sequences are desired 

### *Example*
To obtain sequences for all transcripts predicted in the gtf file, the command might be:

    ./get_all_transcript_sequences.py -gtf merged.gtf -genome genome.fa

To obtain aligned transcript sequence for selected genes, the command might be:

    ./get_transcript_sequences_aligned.py -gtf merged.gtf -genome genome.fa -genes genes.txt

## Output
The `get_all_transcript_sequences.py` script produces one output file, called `all_transcripts.fa`. This is a multi-FASTA file containing sequences for all the transcripts predicted in the gtf file.

The `get_transcript_sequences_aligned.py` script produces multiple output files, one for each gene listed in the input text file. Each is a multi-FASTA files containing sequences for the genomic sequence followed by all the transcripts predicted in the gtf file for that gene. The transcript sequences have dashes for nucleotides in introns, so that all sequences line up with each other when imported into programs for visualisation (e.g. [Geneious](http://www.geneious.com/)). 
