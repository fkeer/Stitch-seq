# Stitch-seq
Code base for Stitch-seq  
Written by Frances Keer, some initial code written by Victoria Popic.  
Advisors: Aziz Al'Khafaji, Paul Blainey  


Run the Snakefile on single-end fastqs, then load into Jupyter notebooks for analysis. 
CD81 example files included for demo. 


Install instructions: Install all dependencies. Create directory containing single-end fastq files. Create misc directory with necessary supplementary files (grnas.txt for whitelists, fasta files for all targets). In misc folder, use bowtie2 to generate a custom reference for the targets based on the fasta files. Run Snakefile from upper directory. Jupyter notebook: transfer counts_results files to a folder containing a labeled gRNA CSV and run. 

Install time: 10min


Software Dependencies:
Command line:
Python v3.10
cutadapt v4.4
bowtie2 v2.5.2
samtools v1.17
Snakemake v7.28.2
tqdm v4.65.0
pigz v2.8
polyleven v0.8
pysam v0.21.0
biopython v1.81
NGmerge v0.3

Jupyter notebook v7.3.3:
Python v3.13.3
anndata v0.12.2
gseapy v1.1.9
matplotlib v3.10.1
numpy v2.2.4
pandas v2.2.3
scanpy v1.11.4
scipy v1.15.2
seaborn v0.13.2
sklearn v1.7.1
statsmodels v0.14.5


Demo:
Follow installation instructions. Copy CD81.fastq, grnas.txt, and cd81-mrna.fasta to misc folder. Use bowtie2 to make custom reference (bowtie2-build cd81.fasta goi). Run Snakefile on CD81.fastq (snakemake -p --cores all results/cd81_counts_report-mrna.txt).  
