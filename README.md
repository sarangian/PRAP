# PRAP

##1.Introduction
PRAP is a platform independent Python3 tool used to analyze pan-resistome characteristics for multiple genomes. It accepts various format of sequence files as input files. The identification of antibiotic resistance genes is mainly based on the Comprehensive Antibiotic Resistance Database (CARD) and ResFinder database.

##2.What can PRAP do
1)	Antibiotic resistance genes (ARGs) identification
2)	Pan-resistome feature analysis
3)	Classifying and analyzing for identified ARGs
4)	Analysis of ARGs associated with given antibiotics

##3.	Installation
###[1] Install Python3 v3.6+ for Win/Mac/Linux (https://www.python.org/)

###[2] Install Python3 accessory packages:
For PRAP, several Python3 packages are required. 
####Biopython v1.7+
####NumPy v1.15+
####Pandas v0.23+
#####SciPy v1.1+
#####Matplotlib v3.0+
#####Seaborn v0.9
#####Scikit-learn v0.19+

###[3] Install blast v2.7.1+ for Win/Mac/Linux 
Blast+ is available at https://blast.ncbi.nlm.nih.gov/Blast.cgi
**you need to modify the directory of where you install blast**
For example, if the blastn programs are in C:/blast+/bin
please change the directory of “blast+” in “settings.txt” like:
[blast+=C:/blast+/bin/] #please use “/” to separate the directory
