# **"Goliath" Scripts**
## This repository contains various small programming efforts that I have used or currently use in my genomics and bioinformatics research. Feel free to check them out and use as needed if you think they might be useful to you.

## **Current Scripts:**

### **sequenceEdits**
<p>Python script developed to aid Monterey County Public Health Lab in manual curation of COVID-19 genome sequences generated on Clear Labs ONT sequencing platform, but were rejected from GISAID. The program is meant to provide some brief summary data on ambiguous bases, provide a means of automating some of the manual assessment that typically goes into variant/sequence curation, and supports manual editing of the query sequence based using positions from a mulitiple sequence alignment. It use MAFFT alignment of reference (typically Wuhan strain) and query sequence (assembled genome) and base coverage frequencies at variable sites (calculated from bam file of sequencing reads mapped back to reference genome) to generate a predicted base call for all variable sites. Dependencies include a local installation of MAFFT, Biopython, Pysam, and Numpy.</p>

### **fastaExtract.py**
<p>Python script that you can use to extract sequences from a larger file by feeding the program a list of sequences that you either was to include or exclude from your extraction. Meant to be run on command line and output needs to be redirected (with a '>') to a new file to be saved. 
  Standard usage is aa follows: python fastaExtract.py -i/e -list <list of sequences you either want to extract or avoid> -fasta <path to fasta file></p>
  
  ### **GCSkew.py***
  <p>Python script that you can feed a genome and window size and get percent GC and GC skew calculations for non-overlapping, window-sized intervals from a genome. Useful for creating Circos plots, especially in Circa.</>
