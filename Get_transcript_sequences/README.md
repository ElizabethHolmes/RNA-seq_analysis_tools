# Get_transcript_sequences
## About
Get_transcript_sequences is a pair of Python scripts for obtaining transcript sequences from the output of the Tuxedo suite of tools for RNA-seq analysis [(Trapnell *et al.*, 2012)](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html), though it could also ge used on gtf files from other sources. The `get_all_transcript_sequences.py` script obtains all predicted transcript sequences, which might be useful for e.g. creating a BLAST database to BLAST other sequences against. The `get_transcript_sequences_aligned.py` script obtains predicted transcript sequences in an aligned format for genes of interest, which might be useful for e.g. primer design to test genes of interest.

** PLEASE NOTE: The Get_trancript_sequences scripts are some of my first scripts written as a beginner programmer and so I apologise if the code is inelegant, unconventional or otherwise sub-optimal; they work for the intended purpose and I provide them in case they might be useful to others, but with no guarantees. ** 

## Citation
Get_transcript_sequences is not yet associated with a paper; to cite it please use:

    Sutton, ER. (2015). Get_transcript_sequences [Software]. 
    Available at https://github.com/ElizabethSutton/RNA-seq_analysis_tools/Get_transcript_sequences.

## Requirements
Get_transcript_sequences requires Python. It has been tested only with Python 2.7.3. It also requires [bedtools](https://github.com/arq5x/bedtools2). It has been tested only with bedtools version 2.16.2.

## Usage
Both Get_transcript_sequences scripts have two required arguments, and `get_transcript_sequences_aligned.py` has one additional argument detailed below.

### *Required arguments*
* `genome` - path to genome FASTA file
* `gtf` - path to merged.gtf file output by Cuffmerge

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
