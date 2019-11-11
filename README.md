# PRAP

## 1.Introduction
PRAP is a platform independent Python3 tool used to analyze pan-resistome characteristics for multiple genomes. It accepts various format of sequence files as input files. The identification of antibiotic resistance genes is mainly based on the Comprehensive Antibiotic Resistance Database (CARD) and ResFinder database.

## 2.What can PRAP do
1)	Antibiotic resistance genes (ARGs) identification
2)	Pan-resistome feature analysis
3)	Classifying and analyzing for identified ARGs
4)	Analysis of ARGs associated with given antibiotics

## 3.Installation
### [1] Install Python3 v3.6+ for Win/Mac/Linux (https://www.python.org/)

### [2] Install Python3 accessory packages:
**For PRAP, several Python3 packages are required**
#### Biopython v1.7+
#### NumPy v1.15+
#### Pandas v0.23+
#### SciPy v1.1+
#### Matplotlib v3.0+
#### Seaborn v0.9
#### Scikit-learn v0.19+

### [3] Install blast v2.7.1+ for Win/Mac/Linux 
Blast+ is available at https://blast.ncbi.nlm.nih.gov/Blast.cgi 
**you need to modify the directory of where you install blast** 
For example, if the blastn programs are in C:/blast+/bin 
please change the directory of “blast+” in “settings.txt” like: 
[blast+=C:/blast+/bin/] #please use “/” to separate the directory 

## 4.Files Preparation
### [1]Input files
#### Raw reads sequence: FASTQ files (example: A.fastq)
#### Protein sequence: FASTA amino acid files (example: A.faa)
#### Nucleotide sequence: FASTA nucleic acid files (example: A.fna)
#### GenBank annotation files (example: A.gb)

### [2]Setting files
**The setting file in PRAP package (/PRAP/settings.txt) is used to set parameters for individual modules.**

### [3]Phenotype files
#### “ar_phenotype.csv” or "res_phenotype.csv"

## 5.Usage
### [1]Run PRAP

#### A simple example:
**python install_dir/PRAP.py -m module_name -indir input_directory -outdir output_directory**

#### Raw reads as input files (files with a ".fastq" extension):
[R1] analysis with CARD nucleotide database

[R2] analysis with ResFinder nucleotide database  

#### Genbank files as input files (files with a ".gb" extension):
[G1-1] analysis with CARD nucleotide database  

[G1-2] analysis with CARD protein database  

[G2-1] analysis with ResFinder nucleotide database  

[G2-2] analysis with ResFinder protein database  


#### Nucleotide seq as input files (files with a ".fna" extension):
[N1] analysis with CARD nucleotide database  

[N2] analysis with ResFinder nucleotide database  


#### Protein seq as input files (files with a ".faa" extension):
[P1] analysis with CARD protein database  

[P2] analysis with ResFinder protein database  


#### Preprocessing module:
[a] "CDSex.py": extract coding sequence from genbank files;
	form both protein and nucleotide fasta files;

#### Gene identification modules:
[b-1] "ArKmer.py": find resistance genes from raw reads;  

[b-2] "ArBlastn.py": find resistance genes from ".fna" files;  

[b-3] "ArBlastp.py": find resistance genes from ".faa" files;  

[c-1] "ResKmer.py": find resistance genes from raw reads;  

[c-2] "ResBlastn.py": find resistance genes from ".fna" files;  

[c-3] "ResBlastn.py": find resistance genes from ".fna" files;  


#### Analysis modules:
[d] "Pangenome.py": pan-resistome analysis (mainly analyze for pan-genome features)  

[e] "PanAccess.py": pan & accessory resistome analysis (mainly classify and statistical analysis for ARGs)  

[f-1] "ArMatrix.py": analysis associated genes for each kind of antibiotics in CARD database  

[f-2] "ResMatrix.py": analysis associated ARGs for each kind of antibiotics in ResFinder database
