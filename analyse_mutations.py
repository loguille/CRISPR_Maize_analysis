#!/usr/bin/env python3

import os
import argparse
import sys
import os.path
import re

#Change this path to import the different program used from an other directory
sys.path.insert(1, '/projects/SeedDev/CRISPR_NGS_2019/loic/CRISPMais/Projet_CRISPR_MAIS')

import countFrequenceFromSRS
import tag_handling
import transform_tag
import creation_logo
import countFrequenceFromSRS_with2PAM
import compare_srs_file


parser = argparse.ArgumentParser()
parser.add_argument('-s','--score', help = 'choose a minimal score of alignment between 0 and 350',default=200,choices=range(0,351),metavar="[0-350]" ,type = int)
parser.add_argument('-f','--forward', help = 'path to the forward sequence file (required)', type = str,required=True)
parser.add_argument('-r','--reverse', help = 'path to the reverse sequence file (required)',type = str, required=True)
parser.add_argument('-t','--target',help= 'path to the target sequence (required)',type =str, required=True,nargs='+')
parser.add_argument('-p','--primer',help ='path to the primers sequence in fasta format (required) as to be the same number of file as the sequence target',type = str , required =True,nargs='+')
parser.add_argument('-pam',help ='position of the first base of the pam sequence as to be calculate from the sequence without the primers (required)',type = int,required = True,nargs='+')
parser.add_argument('-a','--adapters',help ='path to the adapters (required)',type = str,required = True)
parser.add_argument('-o','--directory', help ='name of results directory',type =str, default = 'output_program')
parser.add_argument('-cutoff',help = 'Cutoff of frequency for the logo', type = float, default =0.0001)
parser.add_argument('-d','--direction',help='Give the direction of the pam sequence (use in order to calculate the CRISPR range', choices=['sens','antisens','forward','reverse'],required = True,nargs='+',type = str)


args = parser.parse_args()

#raise an exception if pam position as more than two number 
if len(args.pam) > 2 :
    print("Pam argument take one or two position but not more")
    sys.exit(1)

score = args.score
read1 = args.forward
read2 = args.reverse
sequence_target = args.target
primers = args.primer
pam_position = args.pam
adapters = args.adapters
directory = args.directory
cutoff = args.cutoff
direction = args.direction

#raise an error if the number of primer and sequence target are not equal
if len(primers) != len(sequence_target):
    print("Number of primers and sequence target file as to be the same")
    sys.exit(1)

# raise an error if the cutoff is superior to 1
if cutoff > 1 :
    print("Cutoff cannot be superior to 1")
    sys.exit(1)

# raise an error if the number of pam position and the number of direction are not equal
if len(pam_position) != len(direction):
    print("Number of pam position and direction have to be the same")
    sys.exit(1)

#get working directory

path=os.path.abspath(sys.argv[0])
scriptDir=os.path.dirname(path)


#Path to the different program used 
PathToPEAR='/usr/local/bin/pear'
PathToFASTQ_MCF='/usr/bin/fastq-mcf'
PathToFASTQC='/usr/bin/fastqc'
PathToBowtieBuild='/usr/bin/bowtie2-build'
PathToBowtie2='/usr/bin/bowtie2'
PathToMatrix='/home/jjust01/src/napoly/repeats/bin/EDNAFULL83'

#function to count the number of reads containing a N in the sequence 
def count_number_sequence_withN(fileIN):
    regex=re.compile(r"^[ATCGN]*N[ATCGN]*$")
    count=0
    for i in open(fileIN,'r'):
        if regex.match(i):
            count+=1
    return(count)

#function to count the number of reads in a fasta format
def count_number_sequence(fileIN):
    regex=re.compile(r"^[ATCGN]*$")    
    count=0
    for i in open(fileIN,'r'):
        if regex.match(i):
            count+=1
    return(count)

#function to count the number of line in a file 
def count_nb_lines(fileIN):
    count = -1
    for i in open(fileIN,'r'):
        count += 1
    return(count)

#create log directory
path_log=scriptDir+'/'+directory

if os.path.exists(path_log) :
    if len(os.listdir(path_log)) != 0 :
        print('\n#############################\nWARNING : DIRECTORY NOT EMPTY\n#############################\n')
else :
    os.makedirs(path_log,exist_ok=True)

#get filename W/O extension
filename=os.path.basename(read1)
filename_wo_extension=filename.split('.')[0]

print('\n######################\nSTARTING PEAR\n######################\n')
os.makedirs(path_log+"/PEAR", exist_ok =True)

Pear_output_file=path_log+'/PEAR/'+filename_wo_extension

#assembled paired end sequences
os.system(PathToPEAR+" -f "+read1+" -r "+read2+" -o "+Pear_output_file+" -y 2G -p 0.001 -j10")

#count the percent of reads unassembled
count_of_lines_Pear_output = sum(1 for l in open(Pear_output_file+'.assembled.fastq'))
extension=os.path.splitext(read1)
if extension[1]=='.gz':
    os.system("gunzip "+read1)
    read1=extension[0]
    count_of_lines_Input_file=sum(1 for l in open(read1))
    os.system("gzip "+read1)
    read1=sys.argv[1]
else:
    count_of_lines_Input_file=sum(1 for l in open(read1))

file_count_percent_of_reads_kept=path_log+"/stderr.out"

if int(count_of_lines_Pear_output)/int(count_of_lines_Input_file)<0.87:
    print("Percent of assembly read is inferior to the cutoff (0.87), check the manual of PEAR to see if a solution exist")
    sys.exit(1)
else :
    with open(file_count_percent_of_reads_kept,'w')as read_kept:
        read_kept.write("########################\nAssembly step using PEAR\n########################\n")
        read_kept.write("Number of lines from input file : "+str(count_of_lines_Input_file/4)+'\n')
        read_kept.write("Number of reads from PEAR output file : "+str(count_of_lines_Pear_output/4)+'\n')
        read_kept.write("Percent reads assembly: "+str(count_of_lines_Pear_output/count_of_lines_Input_file)+"\n")



print('\n######################\nSTARTING FASTQC AND CUT ADAPTERS\n######################\n')
os.makedirs(path_log+"/QC", exist_ok=True)
Input_file_QC=path_log+"/PEAR/"+filename_wo_extension+".assembled.fastq"
Output_file_html=filename_wo_extension+".assembled.trimmed.fastq.html"
Output_file_trimmed=filename_wo_extension+".assembled.trimmed.fastq"

print("#### Trimming phase ####")
os.system(PathToFASTQ_MCF+" -l 100 -o "+path_log+"/QC/"+Output_file_trimmed+" "+adapters+" "+Input_file_QC+" -q 0")
print("#### FASTQC phase ####")
os.system(PathToFASTQC+" "+path_log+"/QC/"+Output_file_trimmed)

number_of_reads_wN_kept=count_number_sequence(path_log+"/QC/"+Output_file_trimmed)
number_of_reads_from_PEAR=count_number_sequence(Pear_output_file+".assembled.fastq")
with open(file_count_percent_of_reads_kept,'a') as read_kept:
    read_kept.write("#######################\nTrimming with FASTQ-MCF\n#######################\n")
    read_kept.write("Number of reads kept after trimming : "+str(number_of_reads_wN_kept)+'\n')
    read_kept.write("Percent of reads kept after trimming :"+str(number_of_reads_wN_kept/number_of_reads_from_PEAR)+"\n")
    read_kept.write("Report FASTQC can be found in QC/ repository\n")



print('\n###############\nTAG HANDLING\n#################\n')
os.makedirs(path_log+"/TAG", exist_ok=True)
Input_file_tag=path_log+"/QC/"+filename_wo_extension+".assembled.trimmed.fastq"
Output_file_with_N=path_log+"/QC/"+filename_wo_extension+"withN.assembled.trimmed.fastq"
Output_file_tag=path_log+"/TAG/"+filename_wo_extension+".assembled.trimmed.untagged.fastq"


number_of_sequencewithN=count_number_sequence_withN(Input_file_tag)
number_of_sequence_total=count_number_sequence(Input_file_tag)
if number_of_sequencewithN > 0 :
    os.system(" grep -B1 -A2 ^[ATCGN]*N[ATCGN]*$ "+Input_file_tag+" > temporary_file")
    os.system(" sed '/^--$/d' temporary_file > "+Output_file_with_N)
    os.system(" grep -B1 -A2 ^[ATCG]*$ "+Input_file_tag+"> temporary_file")
    os.system( "sed '/^--$/d' temporary_file > "+Input_file_tag)
    os.system("rm temporary_file")
    with open(file_count_percent_of_reads_kept,'a') as read_kept:
        read_kept.write("############\nTAG handling\n############\n")
        read_kept.write("Number of sequence found with a N in it (these sequences will be eliminated for the analyse but kept in a file in the QC folder ):"+str(number_of_sequencewithN)+"\n")
        if number_of_sequencewithN/number_of_sequence_total > 0.01 :
            read_kept.write("Warning number of sequence with a N in it is superior to 1% :"+str(number_of_sequencewithN/number_of_sequence_total)+"\n")
        else:
            read_kept.write("Percent of sequence with a N : "+str(number_of_sequencewithN/number_of_sequence_total)+"\n")
else :
    with open(file_count_percent_of_reads_kept,'a') as read_kept:
        read_kept.write("No sequence with N found\n")

temporary_file = path_log+"/TAG/temp.txt"
output_file_count_tag = path_log+"/TAG/"+filename_wo_extension+".count_untagged.fasta"
os.system(r'''grep  ^[ATCGN]*$ '''+Input_file_tag+r''' |sort|uniq -c|awk -F' ' 'BEGIN {OFS=""} {print $2,":",$1}'|sed 's/^\ *//g' > '''+temporary_file)
tag_handling.handle_tag(temporary_file,output_file_count_tag)
os.system('rm '+temporary_file)

nb_unique_seq = count_nb_lines(output_file_count_tag)
with open(file_count_percent_of_reads_kept,'a') as read_kept:
    read_kept.write("Number of unique sequence : "+str(nb_unique_seq)+"\nThis file can be found in TAG/ repository\n")


print('handling tags done')
print('\n###############\nCOLLAPSE AND COUNTING TAGS\n#################\n')
os.makedirs(path_log+'/collapser', exist_ok =True)


CollapseFile=path_log+'/collapser/'+filename_wo_extension+".assembled.trimmed.untagged.collapser.fasta"
outputFile=path_log+'/count/'+filename_wo_extension+".collapseTagcount.txt"

transform_tag.write_file(output_file_count_tag, CollapseFile)

print('collapse and counting tags done')

print('\n###############\nALIGNMENT WITH NEEDLEMAN & WUNSCH\n#################\n')
#repository creation
os.makedirs(path_log+'/Aligner',exist_ok = True)
fileINA=path_log+'/collapser/'+filename_wo_extension+'.assembled.trimmed.untagged.collapser.fasta'


#We realigne the sequence with 2 or more insertions/deletions with a needleman & wunsch algorithm
for seq_tar in sequence_target :
    print(seq_tar)
    seq_name = os.path.basename(seq_tar).split('.')[0]
    sequence_realigned=path_log+'/Aligner/'+seq_name+'.aligned.needle.srs'
    print(sequence_realigned)
    os.system("needle -snucleotide1 Y -sformat fasta -asequence "+seq_tar+" -snucleotide2 Y -sformat2 fasta -bsequence "+fileINA+" -datafile "+PathToMatrix+" -gapopen 10.0 -gapextend 0.5 -awidth 99999 -outfile "+sequence_realigned)



with open(file_count_percent_of_reads_kept,'a') as read_kept:
    read_kept.write("#################\nAlignment with NW\n#################\nAlignment with Needlemman & Wunsch done, the alignment file can be found in Aligner/ repository\n")

if len(sequence_target) > 1:
    aln = {}
    for srs_file in os.listdir(path_log+'/Aligner/') :
        file_srs = path_log+'/Aligner/'+srs_file
        new_file_srs = path_log+'/Aligner/'+srs_file.split('.')[0]+'.best_ali_score.srs'
        aln = compare_srs_file.readfiles(aln,file_srs,new_file_srs)
        os.system('rm '+file_srs)
    compare_srs_file.compare_score(aln)


print('\n###################\nFREQUENCY CALCULATION\n######################')
#repository creation
os.makedirs(path_log+'/frequence', exist_ok =True)
frequence_mutations_file = path_log+'/frequence/'+filename_wo_extension+'_frequence_mut.txt'

os.makedirs(path_log+'/logo', exist_ok = True)
with open(file_count_percent_of_reads_kept,'a') as read_kept :
    read_kept.write("#####################\nFrequency calculation\n#####################\n")

for i in range(len(os.listdir(path_log+'/Aligner/'))) :
    frequence_mutations_file = path_log+'/frequence/'+os.path.basename(sequence_target[i]).split('.')[0]+'_frequence_mut.txt'
    print(frequence_mutations_file)
    fileIn = path_log+'/Aligner/'+os.listdir(path_log+'/Aligner/')[i]
    if len(pam_position) == 2:
        countFrequenceFromSRS_with2PAM.write_mutations(frequence_mutations_file,sequence_target[i],fileIn,pam_position,primers[i],file_count_percent_of_reads_kept,score,direction)
    else :
        countFrequenceFromSRS.write_mutations(frequence_mutations_file,sequence_target[i],fileIn,pam_position,primers[i],file_count_percent_of_reads_kept,score,direction) 
    long_file_mut = path_log+'/frequence/'+os.path.basename(sequence_target[i]).split('.')[0]+'_frequence_mut.long.txt'
    nb_mut = count_nb_lines(long_file_mut)
    with open(file_count_percent_of_reads_kept,'a') as read_kept:
        read_kept.write("Number of unique mutations : "+os.path.basename(sequence_target[i]).split('.')[0]+' '+str(nb_mut)+'\n')
    output_logo = path_log+'/logo/'+os.path.basename(sequence_target[i]).split('.')[0]+'logo.svg'
    creation_logo.construct_logo(long_file_mut,primers[i],sequence_target[i],cutoff,output_logo,pam_position,direction)





print('\n###################\nCREATION LOGO\n#####################')

