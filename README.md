# Analyse mutations 
This pipeline has been set up in order to compute the mutations and their frequencies caused by CRISPR from a collection of different sample of maize protoplast. This pipeline is in *python* and can be used with *python 3.0*

## Program version
This pipeline use different other tools, here is the version of the different software used. 
* PEAR : version 0.9.10
* FastQC : version 0.11.7
* Fastq-mcf : version 1.04.676
* Needle (EMBOSS) : version 6.6.0.0

## Argument of the pipeline
This pipeline take some arguments as entry some are needed and other are optional. 
* -h or --help : show help message and exit
* -f or --forward (__required__): forward file of the reads for the assembly 
* -r or --reverse (__required__): reverse file of the reads for the assembly
* -a or --adapters (__required__): fasta sequence of the adapters
* -t or --target (__required__): fasta sequence of the target sequence for the alignment
* -p or --primer (__required__): fasta sequence of the primer (as to be the same number of file as the sequence target)
* -pam (__required__): position of the first base of the pam sequence (used to know the crispr range)
* -o or --directory (optional): path to the output results (default : /output_program/)
* -cutoff (optional)(0-1): cutoff for the frequency logo (default : 0.0001)
* -score (optional)(0-350): minimal alignment score required for the analyse (default : 200), can't be superior to 350
* -d or --direction (__required__) (sens,antisens,forward,reverse): direction of the pam sequence (as to be the same number of file as the number of pam position)

## Code line 
`./Projet_CRISPR_MAIS/analyse_mutations.py -f NP06_S6_R1_001.fastq -r NP06_S6_R2_001.fastq -t NP06_ZmKAK1_KAK1_233_L1750.fasta -a data/adaptors.fa -p primers_NP06.txt -pam 93 -d sens`
