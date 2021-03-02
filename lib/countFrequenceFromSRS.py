#!/usr/bin/python

import sys
import re
import calcul_new_score
from collections import OrderedDict


class Mutation():
    def __init__(self,type,position,sequence):
        self.type = type
        self.position = position
        self.sequence = sequence
        self.frequence = 0

class Alignment():
    def __init__(self,alignmentseqcible,alignment,alignmentread):
        self.alignmentseqcible = alignmentseqcible
        self.alignment = alignment
        self.alignmentread = alignmentread

class Read():
    def __init__(self,seq):
        self.seq =  seq 
        self.mutation = []
        self.freq = 0
        self.correct_freq = 0

    def addMutation(self,type,position,sequence):
        self.mutation.append(Mutation(type,position,sequence))


def compareAlign(fileSRS,score_thr):
    '''
    This function compare the two alignments and can print the different mutation existing between 
    the two sequences
    '''
    ali, nbseq, corrected_nb_seq = calcul_new_score.parseFilesrs(fileSRS)
    read_mut={}
    seq_w_low_score = 0
    nb_alignment =0
    for read in ali:
        insertion = ''
        deletion = ''
        nb_alignment += int(ali[read].number_of_seq)
        if ali[read].new_score >= score_thr :
            target_seq = ali[read].align_ref
            aln = ali[read].align
            read_seq = ali[read].align_seq
            read_mut[ali[read].align_seq]=Read(ali[read].align_seq)
            read_mut[ali[read].align_seq].freq += int(ali[read].number_of_seq)
            read_mut[ali[read].align_seq].correct_freq += int(ali[read].correct_nb_seq)
            #here we add the number of sequence to the class Read 
            length_insertion = 0
            '''
            Here we create two different string that contain the insertion or deletion from the target
            sequence, we print this two strings if they are not empty, for the position we add +2 the first base 
            of insertion/deletion correspond to a number in the sequence
            maximal indice for a mismatch/deletion is the last base of the target sequence, and the last base +1 
            for an insertion 
            '''
            for i in range(len(aln)):
                if aln[i] == '.' :
                    '''
                    We add insertion and deletion to the mutation only if the position doesn't correspond to a blank
                    and if insertion or deletion is not empty. 
                    if we don't do that we will have many insertion or deletion one character we have to wait for the
                    end of the mutation
                    '''
                    if insertion :
                        length_insertion += len(insertion)
                        '''
                        length_insertion is used to know the correct position in the reference sequence
                        '''
                        read_mut[ali[read].align_seq].addMutation("insertion",position-length_insertion+2,insertion)
                        '''
                        here we add the mutation to the class Read with the function addMutation which create 
                        class Mutation composed of different information showed above
                        '''
                        insertion = ''
                    elif deletion :
                        read_mut[ali[read].align_seq].addMutation("deletion",position-length_insertion-len(deletion)+2,deletion)
                        deletion = ''
                    mutation = read_seq[i]
                    if length_insertion == 0 :
                        read_mut[ali[read].align_seq].addMutation('mismatch',i+1-length_insertion,mutation)
                    else :
                        read_mut[ali[read].align_seq].addMutation('mismatch',i+1-length_insertion,mutation)
                elif aln[i] == ' ' :
                    if target_seq[i] == '-' :
                        insertion += read_seq[i]
                        position = i
                    else :
                        deletion += target_seq[i]
                        position = i
                else :
                    if insertion :
                        length_insertion += len(insertion)
                        read_mut[ali[read].align_seq].addMutation("insertion",position-length_insertion+2,insertion)
                        insertion = ''                
                    elif deletion :
                        read_mut[ali[read].align_seq].addMutation("deletion",position-length_insertion-len(deletion)+2,deletion)
                        deletion = ''
                    else :
                        continue
            if insertion :
                length_insertion += len(insertion)
                read_mut[ali[read].align_seq].addMutation("insertion",position-length_insertion+2,insertion)
                insertion = ''
            if deletion :
                read_mut[ali[read].align_seq].addMutation("deletion",position-length_insertion-len(deletion)+2,deletion)
                deletion = ''
        else :
            seq_w_low_score += int(ali[read].number_of_seq)
    print('Number of sequence with too low score : '+str(seq_w_low_score))
    return(read_mut, nbseq, corrected_nb_seq, seq_w_low_score, nb_alignment)
                


def concat_mutation_2(dict_read):
    '''
    This function create a dictionnary of list that contains each different mutations for each position
    and the number of time a mutation exists in our data 
    '''
    dico_mutation={}
    for element in dict_read :
        for i in dict_read[element].mutation:
            if i.position not in dico_mutation:
                dico_mutation[i.position] = [[i.type,i.sequence,dict_read[element].freq,dict_read[element].correct_freq]]
            else :
                value = False
                for j in range(len(dico_mutation[i.position])) :
                    if dico_mutation[i.position][j][1] == i.sequence :
                        if dico_mutation[i.position][j][0] == i.type :
                            '''
                            case where this is the same type of mutations at the same position, we add the 
                            number of reads
                            '''
                            dico_mutation[i.position][j][2] += dict_read[element].freq
                            dico_mutation[i.position][j][3] += dict_read[element].correct_freq
                            value = True
                            break
                            '''
                            In this case we don't have to continue because we know that the mutation already exist
                            and we have added the number of reads that have this mutations, so we can break the loop
                            '''
                        else :
                            '''
                            we continue to see if the case before exists  
                            '''
                            continue
                if value == False :
                    '''
                    if the case doesn't exists we add a new mutation 
                    '''
                    dico_mutation[i.position].append([i.type,i.sequence,dict_read[element].freq,dict_read[element].correct_freq])                    
    return(dico_mutation)


def count_frequence(mutations,nbSeq):
    '''
    This function can count the frequence of the different mutation based on the number of reads that contains 
    the mutation on the total number of reads
    '''
    for keys in mutations.keys():
        for i in range(len(mutations[keys])):
            mutations[keys][i][2] = float(mutations[keys][i][2])/float(nbSeq)
    return(mutations)

def open_target_file(sequence_target):
    '''
    This function can read in the file of sequence and return the complete reference sequence
    '''
    with open(sequence_target,'r') as infile:
        for line in infile: 
            line=line.rstrip()
            if not line.startswith('>'):
                sequence = line
    return(sequence)

def open_primers_file(primers_file):
    '''
    This function read the primers file and return a dictionnary that contain the sequence
     of the primer forward and reverse
    '''
    primers_dict={}
    with open(primers_file,'r') as primer:
        for line in primer :
            line = line.rstrip()
            if line.startswith('>'):
                typePrimers=line.split('_')[1]
                if typePrimers not in primers_dict :
                    primers_dict[typePrimers] = []
                continue
            primers_dict[typePrimers].append(line)
    return(primers_dict)


def return_starting_ending_position(pam_position,direction):
    if direction == 'sens' or direction == 'forward':
        start = pam_position-20
        end = pam_position
    else :
        start = pam_position+4
        end = pam_position+24
    return(start,end)


def write_mutations(fileOut, sequence,srs_file,pam_position,primers_file, file_error,score,direction):
    '''
    This function write the different mutations in one compress output file and one extended output file
    sequence_target correspond to the sequence of reference used to align the read.
    '''
    Pam_position1 = pam_position[0]
    direction1 = direction[0]
    dict_read, nbsequence, correctedNbseq, seq_score_to_low, nb_ali = compareAlign(srs_file,score)
    mut = concat_mutation_2(dict_read) 
    mut_dict = count_frequence(mut,nbsequence)
    sequence_target = open_target_file(sequence)
    with open(file_error,'a') as error : 
        error.write('Number of sequences with too low score : '+str(seq_score_to_low)+'\n')
        error.write('Frequence of sequences with too low score : '+str(float(seq_score_to_low/nb_ali))+'\n')
    '''
    we use the function open_primers_file to get a dictionnary that have as key forward and reverse and as 
    value the sequence of the primers
    '''
    primers = open_primers_file(primers_file)
    for key in primers :
        if key == 'forward':
            len_primer_forward = len(primers[key][0])
        else :
            len_primer_reverse = len(primers[key][0])
    sort_dico_mut = {}
    sort_dico_mut = OrderedDict(sorted(mut_dict.items(), key = lambda t: t[0]))
    start_pam1,end_pam1 = return_starting_ending_position(Pam_position1,direction1)

    '''
    First we have to sort the dictionnary by the key to have the different mutation by ascending position
    '''
    insertion = ''#get the different insertion
    insertion_tot_freq = 0#get the total frequency
    insertion_freq = ''#get all the frequencies separated by a coma
    corrected_ins_freq = ''#get all the corrected frequencies separated by a coma 
    corrected_ins_tot_freq = 0#get the total corrected frequency
    corrected_ins_count = ''#get the corrected count separated by a coma
    corrected_ins_tot_count = 0#get the total corrected count 
    list_length_ins = ''

    '''
    The seven variables below are used to construct the two files containing concatened versions of yes or no for 
    the crispr range
    We fill these variables with the same information as the concatened file. 
    '''
    insertion_in_crispr_range = ''  #Here we use the same type of variable as before but with in and out the crispr range
    insertion_tot_freq_in_crispr_range = 0
    insertion_freq_in_crispr_range = ''
    len_ins_in_crispr_range = 0
    list_length_ins_in_crispr_range = ''
    nb_ins_in_crispr_range = 0
    count_ins_in_crispr_range = ''
    total_count_ins_in_crispr_range = 0
    corrected_freq_ins_in_crispr_range = ''
    corrected_tot_freq_ins_in_crispr_range = 0
    corrected_count_ins_in_crispr_range = ''
    corrected_tot_count_ins_in_crispr_range = 0

    insertion_out_crispr_range = ''
    insertion_tot_freq_out_crispr_range = 0
    insertion_freq_out_crispr_range = '' 
    len_ins_out_crispr_range = 0
    list_length_ins_out_crispr_range = ''
    nb_ins_out_crispr_range =0
    count_ins_out_crispr_range = ''
    total_count_ins_out_crispr_range = 0 
    corrected_freq_ins_out_crispr_range = ''
    corrected_tot_freq_ins_out_crispr_range = 0
    corrected_count_ins_out_crispr_range = ''
    corrected_tot_count_ins_out_crispr_range = 0

    #Var used to fill the lines containing mismatch for the different concatened file
    mismatch = ''
    mismatch_tot_freq = 0
    mismatch_freq = ''
    corrected_mis_freq = ''
    corrected_mis_tot_freq = 0
    corrected_mis_count = ''
    corrected_mis_tot_count = 0
    list_length_mis = ''

    mismatch_in_crispr_range = ''
    mismatch_tot_freq_in_crispr_range = 0
    mismatch_freq_in_crispr_range = ''
    count_mismatch_in_crispr_range = ''
    total_count_mis_in_crispr_range =0
    corrected_freq_mis_in_crispr_range = ''
    corrected_tot_freq_mis_in_crispr_range = 0
    corrected_count_mis_in_crispr_range = ''
    corrected_tot_count_mis_in_crispr_range = 0
    list_length_mis_in_crispr_range = ''

    mismatch_out_crispr_range = ''
    mismatch_tot_freq_out_crispr_range = 0
    mismatch_freq_out_crispr_range = ''
    count_mismatch_out_crispr_range = ''
    total_count_mis_out_crispr_range = 0
    corrected_freq_mis_out_crispr_range = ''
    corrected_tot_freq_mis_out_crispr_range = 0
    corrected_count_mis_out_crispr_range = ''
    corrected_tot_count_mis_out_crispr_range = 0
    list_length_mis_out_crispr_range = ''

    #var used to fill the lines containing deletion for the different concatened file 
    deletion = ''
    deletion_tot_freq = 0
    deletion_freq = ''
    corrected_del_freq = ''
    corrected_del_tot_freq = 0
    corrected_del_count = ''
    corrected_del_tot_count = 0
    list_length_del = ''

    deletion_in_crispr_range = ''
    deletion_tot_freq_in_crispr_range = 0
    deletion_freq_in_crispr_range = ''
    len_del_in_crispr_range = 0
    list_length_del_in_crispr_range = ''
    nb_del_in_crispr_range = 0
    count_del_in_crispr_range = ''
    total_count_del_in_crispr_range = 0
    corrected_freq_del_in_crispr_range = ''
    corrected_tot_freq_del_in_crispr_range = 0
    corrected_count_del_in_crispr_range = ''
    corrected_tot_count_del_in_crispr_range = 0

    deletion_out_crispr_range = ''
    deletion_tot_freq_out_crispr_range = 0
    deletion_freq_out_crispr_range = ''
    len_del_out_crispr_range = 0
    list_length_del_out_crispr_range = ''
    nb_del_out_crispr_range = 0
    count_del_out_crispr_range = ''
    total_count_del_out_crispr_range = 0
    corrected_freq_del_out_crispr_range = ''
    corrected_tot_freq_del_out_crispr_range = 0
    corrected_count_del_out_crispr_range = ''
    corrected_tot_count_del_out_crispr_range = 0

    len_ins = 0 #get the length of the different insertion 
    len_del = 0 #get the length of the different deletion
    nb_ins = 0 #count the number of insertion
    nb_del = 0 #count the number of deletion 
    count_ins = ''
    count_del = ''
    count_mismatch = ''
    crispr_range_mism = '' #contain yes or no or both of us for the concatenate file 
    crispr_range_ins = ''
    crispr_range_del = ''
    bool_crispr_range = False #bool in order to know if the data range for crispr contain only No or Yes or both of us 
    sum_crisp = ''
    total_count_ins = 0 #count the total number of indel or mismatch 
    total_count_del = 0
    total_count_mis = 0
    ''' 
    Here we construct the different name file from the name passed in as entry
    '''
    unconcatened_file = ''
    file_in_crispr_range = ''
    file_out_crispr_range = ''
    list_element = fileOut.split('.')[:-1]
    for l in list_element :
        unconcatened_file += l+'.'
        file_in_crispr_range += l+'.'
        file_out_crispr_range += l+'.'
    unconcatened_file += 'long.txt'
    file_in_crispr_range += 'in_crispr_range.txt'
    file_out_crispr_range += 'out_crispr_range.txt'
    with open(fileOut,'w') as writing_file, open(unconcatened_file,'w') as long_file, open(file_in_crispr_range,'w') as yes_file ,open(file_out_crispr_range,'w') as no_file:
        writing_file.write('Start position\tEnd position\tMutations Size\tMaximal mutation size\tType of mutation\tReference\tSequence read\tFrequency\tTotal frequency\tCount\tTotal Count\tCorrected frequency\tCorrected total frequency\tCorrected count\tCorrected total count\tMutation induced by CRISPR\tSummary mutation induced by CRISPR\n')
        long_file.write('Start position\tEnd position\tType of mutation\tReference\tSequence read\tFrequency\tMutations length\tCount\tCorrected frequency\tCorrected count\tMutation induced by CRISPR\n')
        yes_file.write('Start position\tEnd position\tMutations Size\tMaximal mutation size\tType of mutation\tReference\tSequence read\tFrequency\tTotal frequency\tCount\tTotal Count\tCorrected frequency\tCorrected total frequency\tCorrected count\tCorrected total count\n')
        no_file.write('Start position\tEnd position\tMutations Size\tType of mutation\tMaximal mutation size\tReference\tSequence read\tFrequency\tTotal frequency\tCount\tTotal Count\tCorrected frequency\tCorrected total frequency\tCorrected count\tCorrected total count\n')        
        for keys in sort_dico_mut.keys():
            if keys-len_primer_forward < 1 or keys > len(sequence_target)-len_primer_reverse :
                continue
            else:
                for i in range(len(sort_dico_mut[keys])):
                    if sort_dico_mut[keys][i][0] == 'insertion' :
                        '''
                        20 correspond to the range for the valid position of CRISPR so this range goes from 
                        pam position -20 to pam position
                        '''
                        if keys-len_primer_forward >= start_pam1  and keys-len_primer_forward <= end_pam1 :
                            long_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len(sort_dico_mut[keys][i][1]))+'\t'+sort_dico_mut[keys][i][0]+'\t/\t'+sort_dico_mut[keys][i][1]+'\t'+str(sort_dico_mut[keys][i][2])+'\t'+str(len(sort_dico_mut[keys][i][1]))+'\t'+str(int(sort_dico_mut[keys][i][2]*nbsequence))+'\t'+str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+'\t'+str(sort_dico_mut[keys][i][3])+'\tYes\n')
                            crispr_range_ins += 'Yes,'
                            insertion_in_crispr_range += sort_dico_mut[keys][i][1]+','#fill the variable with 'in' only if the mutation is in the crispr range
                            insertion_tot_freq_in_crispr_range += sort_dico_mut[keys][i][2]
                            insertion_freq_in_crispr_range += str(sort_dico_mut[keys][i][2])+','
                            count_ins_in_crispr_range += str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                            total_count_ins_in_crispr_range += int(sort_dico_mut[keys][i][2]*nbsequence)
                            corrected_freq_ins_in_crispr_range += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                            corrected_tot_freq_ins_in_crispr_range += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                            corrected_count_ins_in_crispr_range += str(sort_dico_mut[keys][i][3])+','
                            corrected_tot_count_ins_in_crispr_range += sort_dico_mut[keys][i][3]
                            list_length_ins_in_crispr_range += str(len(sort_dico_mut[keys][i][1]))+','
                            
                            if len(sort_dico_mut[keys][i][1]) > len_ins_in_crispr_range :
                                len_ins_in_crispr_range = len(sort_dico_mut[keys][i][1])
                                nb_ins_in_crispr_range += 1
                        else :
                            #here we want to concatenate the information about all the insertion that start at the same position but aren't in the crispr range
                            long_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len(sort_dico_mut[keys][i][1]))+'\t'+sort_dico_mut[keys][i][0]+'\t/\t'+sort_dico_mut[keys][i][1]+'\t'+str(sort_dico_mut[keys][i][2])+'\t'+str(len(sort_dico_mut[keys][i][1]))+'\t'+str(int(sort_dico_mut[keys][i][2]*nbsequence))+'\t'+str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+'\t'+str(sort_dico_mut[keys][i][3])+'\tNo\n')
                            crispr_range_ins += 'No,'
                            insertion_out_crispr_range += sort_dico_mut[keys][i][1]+',' # here we fill the object only if the mutation isn't in the crispr range 
                            insertion_tot_freq_out_crispr_range += sort_dico_mut[keys][i][2] #we add the frequency to have the total frequency at the end for all the insertion that start at the same position 
                            insertion_freq_out_crispr_range += str(sort_dico_mut[keys][i][2])
                            count_ins_out_crispr_range +=str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                            total_count_ins_out_crispr_range += int(sort_dico_mut[keys][i][2]*nbsequence)
                            corrected_freq_ins_out_crispr_range += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                            corrected_tot_freq_ins_out_crispr_range += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                            corrected_count_ins_out_crispr_range += str(sort_dico_mut[keys][i][3])+','
                            corrected_tot_count_ins_out_crispr_range += sort_dico_mut[keys][i][3]
                            list_length_ins_out_crispr_range += str(len(sort_dico_mut[keys][i][1]))+','

                            if len(sort_dico_mut[keys][i][1]) > len_ins_out_crispr_range :
                                #here we want to have the longest insertion and the number of insertion that start at one position 
                                len_ins_out_crispr_range = len(sort_dico_mut[keys][i][1])
                                nb_ins_out_crispr_range += 1
                        insertion += sort_dico_mut[keys][i][1]+',' #we concatenate all the insertion to a string 
                        insertion_freq += str(sort_dico_mut[keys][i][2])+',' #we concatenate all the frequence to a string
                        insertion_tot_freq += sort_dico_mut[keys][i][2]
                        count_ins += str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                        total_count_ins += int(sort_dico_mut[keys][i][2]*nbsequence)
                        corrected_ins_count += str(sort_dico_mut[keys][i][3])+','
                        corrected_ins_tot_count += sort_dico_mut[keys][i][3]
                        corrected_ins_freq += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                        corrected_ins_tot_freq += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                        list_length_ins += str(len(sort_dico_mut[keys][i][1]))+','
                        
                        if len(sort_dico_mut[keys][i][1]) > len_ins : #keep the longest insertion and the number of insertion for one line
                            len_ins = len(sort_dico_mut[keys][i][1])
                            nb_ins += 1

                    elif sort_dico_mut[keys][i][0] == 'mismatch' :
                        if keys-len_primer_forward >= start_pam1  and keys-len_primer_forward <= end_pam1 :
                            long_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward)+'\t'+sort_dico_mut[keys][i][0]+'\t'+sequence_target[keys-1]+'\t'+sort_dico_mut[keys][i][1]+'\t'+str(sort_dico_mut[keys][i][2])+'\t0\t'+str(int(sort_dico_mut[keys][i][2]*nbsequence))+'\t'+str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+'\t'+str(sort_dico_mut[keys][i][3])+'\tYes\n')
                            crispr_range_mism += 'Yes,'
                            mismatch_in_crispr_range += sort_dico_mut[keys][i][1]+','
                            mismatch_tot_freq_in_crispr_range += sort_dico_mut[keys][i][2]
                            mismatch_freq_in_crispr_range += str(sort_dico_mut[keys][i][2])+','
                            count_mismatch_in_crispr_range += str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                            total_count_mis_in_crispr_range += int(sort_dico_mut[keys][i][2]*nbsequence)
                            corrected_freq_mis_in_crispr_range += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                            corrected_tot_freq_mis_in_crispr_range += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                            corrected_count_mis_in_crispr_range += str(sort_dico_mut[keys][i][3])+','
                            corrected_tot_count_mis_in_crispr_range += sort_dico_mut[keys][i][3]
                            list_length_mis_in_crispr_range += '0,'

                        else:
                            long_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward)+'\t'+sort_dico_mut[keys][i][0]+'\t'+sequence_target[keys-1]+'\t'+sort_dico_mut[keys][i][1]+'\t'+str(sort_dico_mut[keys][i][2])+'\t0\t'+str(int(sort_dico_mut[keys][i][2]*nbsequence))+'\t'+str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+'\t'+str(sort_dico_mut[keys][i][3])+'\tNo\n')
                            crispr_range_mism += 'No,'
                            mismatch_out_crispr_range += sort_dico_mut[keys][i][1]+','
                            mismatch_tot_freq_out_crispr_range += sort_dico_mut[keys][i][2]
                            mismatch_freq_out_crispr_range += str(sort_dico_mut[keys][i][2])+','
                            count_mismatch_out_crispr_range += str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                            total_count_mis_out_crispr_range += int(sort_dico_mut[keys][i][2]*nbsequence)
                            corrected_freq_mis_out_crispr_range += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                            corrected_tot_freq_mis_out_crispr_range += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                            corrected_count_mis_out_crispr_range += str(sort_dico_mut[keys][i][3])+','
                            corrected_tot_count_mis_out_crispr_range += sort_dico_mut[keys][i][3]
                            list_length_mis_out_crispr_range += '0,'

                        mismatch += sort_dico_mut[keys][i][1]+','
                        mismatch_tot_freq += sort_dico_mut[keys][i][2]
                        mismatch_freq += str(sort_dico_mut[keys][i][2])+','
                        count_mismatch += str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                        total_count_mis += int(sort_dico_mut[keys][i][2]*nbsequence)
                        corrected_mis_count += str(sort_dico_mut[keys][i][3])+','
                        corrected_mis_tot_count += sort_dico_mut[keys][i][3]
                        corrected_mis_freq += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                        corrected_mis_tot_freq += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                        list_length_mis += '0,'

                    elif sort_dico_mut[keys][i][0] == 'deletion' :
                        if keys-len_primer_forward >= start_pam1  and keys-len_primer_forward <= end_pam1 or keys-len_primer_forward+len(sort_dico_mut[keys][i][1]) >= start_pam1  and keys-len_primer_forward+len(sort_dico_mut[keys][i][1]) <= end_pam1 :
                            long_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len(sort_dico_mut[keys][i][1]))+'\t'+sort_dico_mut[keys][i][0]+'\t'+sort_dico_mut[keys][i][1]+'\t/\t'+str(sort_dico_mut[keys][i][2])+'\t'+str(-len(sort_dico_mut[keys][i][1]))+'\t'+str(int(sort_dico_mut[keys][i][2]*nbsequence))+'\t'+str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+'\t'+str(sort_dico_mut[keys][i][3])+'\tYes\n')
                            crispr_range_del += 'Yes,'
                            deletion_in_crispr_range += sort_dico_mut[keys][i][1]+','
                            deletion_tot_freq_in_crispr_range += sort_dico_mut[keys][i][2]
                            deletion_freq_in_crispr_range += str(sort_dico_mut[keys][i][2])+','
                            count_del_in_crispr_range += str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                            total_count_del_in_crispr_range += int(sort_dico_mut[keys][i][2]*nbsequence)
                            corrected_freq_del_in_crispr_range += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                            corrected_tot_freq_del_in_crispr_range += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                            corrected_count_del_in_crispr_range += str(sort_dico_mut[keys][i][3])+','
                            corrected_tot_count_del_in_crispr_range += sort_dico_mut[keys][i][3]
                            list_length_del_in_crispr_range += str(len(sort_dico_mut[keys][i][1]))+','

                            if len(sort_dico_mut[keys][i][1]) > len_del_in_crispr_range :
                                len_del_in_crispr_range = len(sort_dico_mut[keys][i][1])
                                nb_del_in_crispr_range += 1
                        else :
                            long_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len(sort_dico_mut[keys][i][1]))+'\t'+sort_dico_mut[keys][i][0]+'\t'+sort_dico_mut[keys][i][1]+'\t/\t'+str(sort_dico_mut[keys][i][2])+'\t'+str(-len(sort_dico_mut[keys][i][1]))+'\t'+str(int(sort_dico_mut[keys][i][2]*nbsequence))+'\t'+str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+'\t'+str(sort_dico_mut[keys][i][3])+'\tNo\n')
                            crispr_range_del += 'No,'
                            deletion_out_crispr_range += sort_dico_mut[keys][i][1]+','
                            deletion_tot_freq_out_crispr_range += sort_dico_mut[keys][i][2]
                            deletion_freq_out_crispr_range += str(sort_dico_mut[keys][i][2])+','
                            count_del_out_crispr_range += str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                            total_count_del_out_crispr_range += int(sort_dico_mut[keys][i][2]*nbsequence)
                            corrected_freq_del_out_crispr_range += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                            corrected_tot_freq_del_out_crispr_range += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                            corrected_count_del_out_crispr_range += str(sort_dico_mut[keys][i][3])+','
                            corrected_tot_count_del_out_crispr_range += sort_dico_mut[keys][i][3]
                            list_length_del_out_crispr_range += str(len(sort_dico_mut[keys][i][1]))+','

                            if len(sort_dico_mut[keys][i][1]) > len_del_out_crispr_range:
                                len_del_out_crispr_range = len(sort_dico_mut[keys][i][1])
                                nb_del_out_crispr_range += 1
                        deletion += sort_dico_mut[keys][i][1]+','
                        deletion_tot_freq += sort_dico_mut[keys][i][2]
                        deletion_freq += str(sort_dico_mut[keys][i][2])+','
                        count_del += str(int(sort_dico_mut[keys][i][2]*nbsequence))+','
                        total_count_del += int(sort_dico_mut[keys][i][2]*nbsequence)
                        corrected_del_count += str(sort_dico_mut[keys][i][3])+','
                        corrected_del_tot_count += sort_dico_mut[keys][i][3]
                        corrected_del_freq += str(float(sort_dico_mut[keys][i][3])/float(correctedNbseq))+','
                        corrected_del_tot_freq += float(sort_dico_mut[keys][i][3])/float(correctedNbseq)
                        list_length_del += str(len(sort_dico_mut[keys][i][1]))+','

                        if len(sort_dico_mut[keys][i][1]) > len_del :
                            len_del = len(sort_dico_mut[keys][i][1])
                            nb_del += 1

                if mismatch :
                    '''
                    Here we want calculate the summary of the data range for crispr, we loop on each data range
                    crispr to see if the list contains only Yes or No or both 
                    '''
                    summ_mism = crispr_range_mism.split(',')
                    res = list(filter(None,summ_mism)) 
                    first_word = res[0]
                    if len(res)-1 > 0 :#if the list contain more than one element we have to look in the list
                        for i in range(1,len(res)):
                            if res[i] != first_word :
                                bool_crispr_range = True 
                                break
                            else : 
                                continue
                    else :#else we have only one element and the summary is easy 
                        bool_crispr_range = False
                    if bool_crispr_range == True :
                        sum_crisp = 'Yes, No' 
                    elif first_word == 'Yes' :
                        sum_crisp = 'Yes'
                    else :
                        sum_crisp = 'No'
                    writing_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward)+'\t'+list_length_mis+'\t'+str(0)+'\tmismatch\t'+sequence_target[keys-1]+'\t'+str(mismatch)+'\t'+str(mismatch_freq)+'\t'+str(mismatch_tot_freq)+'\t'+count_mismatch+'\t'+str(total_count_mis)+'\t'+corrected_mis_freq+'\t'+str(corrected_mis_tot_freq)+'\t'+corrected_mis_count+'\t'+str(corrected_mis_tot_count)+'\t'+crispr_range_mism+'\t'+sum_crisp+'\n')
                    '''
                    Here we just clear the different objects that we have construct before 
                    '''
                    bool_crispr_range = False
                    mismatch = ''
                    mismatch_tot_freq = 0
                    mismatch_freq = ''
                    count_mismatch = ''
                    crispr_range_mism = ''
                    total_count_mis = 0
                    corrected_mis_count = ''
                    corrected_mis_tot_count = 0
                    corrected_mis_freq = ''
                    corrected_mis_tot_freq = 0
                    list_length_mis = ''

                if insertion :
                    summ_ins = crispr_range_ins.split(',')
                    res = list(filter(None,summ_ins)) 
                    first_word = res[0]
                    if len(res)-1 > 0 :
                        for i in range(1,len(res)):
                            if res[i] != first_word :
                                bool_crispr_range = True 
                                break
                            else : 
                                continue
                    else :
                        bool_crispr_range = False
                    if bool_crispr_range == True :
                        sum_crisp = 'Yes, No' 
                    elif  first_word == 'Yes' :
                        sum_crisp = 'Yes'
                    else :
                        sum_crisp = 'No'
                    if nb_ins > 1: #we print the longest base for the ending on the indel if the number of indels is superior to one
                        writing_file.write(str(keys-len_primer_forward)+'\t<'+str(keys-len_primer_forward+len_ins-1)+'\t'+list_length_ins+'\t<'+str(len_ins)+'\tinsertion\t/\t'+str(insertion)+'\t'+str(insertion_freq)+'\t'+str(insertion_tot_freq)+'\t'+count_ins+'\t'+str(total_count_ins)+'\t'+corrected_ins_freq+'\t'+str(corrected_ins_tot_freq)+'\t'+corrected_ins_count+'\t'+str(corrected_ins_tot_count)+'\t'+crispr_range_ins+'\t'+sum_crisp+'\n')
                    else :
                        writing_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len_ins-1)+'\t'+list_length_ins+'\t'+str(len_ins)+'\tinsertion\t/\t'+str(insertion)+'\t'+str(insertion_freq)+'\t'+str(insertion_tot_freq)+'\t'+count_ins+'\t'+str(total_count_ins)+'\t'+corrected_ins_freq+'\t'+str(corrected_ins_tot_freq)+'\t'+corrected_ins_count+'\t'+str(corrected_ins_tot_count)+'\t'+crispr_range_ins+'\t'+sum_crisp+'\n')
                    bool_crispr_range = False
                    insertion = ''
                    insertion_tot_freq = 0
                    insertion_freq = ''
                    nb_ins = 0
                    len_ins = 0
                    count_ins = ''
                    crispr_range_ins = ''
                    total_count_ins = 0
                    corrected_ins_tot_count = 0
                    corrected_ins_tot_freq = 0
                    corrected_ins_freq = ''
                    corrected_ins_count = ''
                    list_length_ins = ''

                if deletion :
                    summ_del = crispr_range_del.split(',')
                    res = list(filter(None,summ_del)) 
                    first_word = res[0]
                    if len(res)-1 > 0 :
                        for i in range(1,len(res)):
                            if res[i] != first_word :
                                bool_crispr_range = True 
                                break
                            else : 
                                continue
                    else :
                        bool_crispr_range = False
                    if bool_crispr_range == True :
                        sum_crisp = 'Yes, No' 
                    elif  first_word == 'Yes' :
                        sum_crisp = 'Yes'
                    else :
                        sum_crisp = 'No'
                    if nb_del > 1:
                        writing_file.write(str(keys-len_primer_forward)+'\t<'+str(keys-len_primer_forward+len_del-1)+'\t'+list_length_del+'\t<'+str(len_del)+'\tdeletion\t'+str(deletion)+'\t/\t'+str(deletion_freq)+'\t'+str(deletion_tot_freq)+'\t'+count_del+'\t'+str(total_count_del)+'\t'+corrected_del_freq+'\t'+str(corrected_del_tot_freq)+'\t'+corrected_del_count+'\t'+str(corrected_del_tot_count)+'\t'+crispr_range_del+'\t'+sum_crisp+'\n')
                    else :
                        writing_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len_del-1)+'\t'+list_length_del+'\t'+str(len_del)+'\tdeletion\t'+str(deletion)+'\t/\t'+str(deletion_freq)+'\t'+str(deletion_tot_freq)+'\t'+count_del+'\t'+str(total_count_del)+'\t'+corrected_del_freq+'\t'+str(corrected_del_tot_freq)+'\t'+corrected_del_count+'\t'+str(corrected_del_tot_count)+'\t'+crispr_range_del+'\t'+sum_crisp+'\n')
                    bool_crispr_range = False
                    deletion_tot_freq = 0
                    deletion = ''
                    deletion_freq = ''
                    nb_del = 0
                    len_del = 0
                    count_del = ''
                    crispr_range_del = ''
                    total_count_del = 0
                    corrected_del_tot_count = 0
                    corrected_del_tot_freq = 0
                    corrected_del_freq = ''
                    corrected_del_count = ''
                    list_length_del = ''

                if mismatch_in_crispr_range  :
                    yes_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward)+'\t'+list_length_mis_in_crispr_range+'\t0'+'\tmismatch\t'+sequence_target[keys-1]+'\t'+str(mismatch_in_crispr_range)+'\t'+str(mismatch_freq_in_crispr_range)+'\t'+str(mismatch_tot_freq_in_crispr_range)+'\t'+count_mismatch_in_crispr_range+'\t'+str(total_count_mis_in_crispr_range)+'\t'+corrected_freq_mis_in_crispr_range+'\t'+str(corrected_tot_freq_mis_in_crispr_range)+'\t'+corrected_count_mis_in_crispr_range+'\t'+str(corrected_tot_count_mis_in_crispr_range)+'\n')
                    mismatch_in_crispr_range = ''
                    mismatch_tot_freq_in_crispr_range = 0
                    mismatch_freq_in_crispr_range = ''
                    count_mismatch_in_crispr_range = ''
                    total_count_mis_in_crispr_range = 0
                    corrected_freq_mis_in_crispr_range = ''
                    corrected_tot_freq_mis_in_crispr_range = 0
                    corrected_count_mis_in_crispr_range = ''
                    corrected_tot_count_mis_in_crispr_range = 0
                    list_length_mis_in_crispr_range = ''

                if mismatch_out_crispr_range != '':
                    no_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward)+'\t'+list_length_mis_out_crispr_range+'\t0'+'\tmismatch\t'+sequence_target[keys-1]+'\t'+str(mismatch_out_crispr_range)+'\t'+str(mismatch_freq_out_crispr_range)+'\t'+str(mismatch_tot_freq_out_crispr_range)+'\t'+count_mismatch_out_crispr_range+'\t'+str(total_count_mis_out_crispr_range)+'\t'+corrected_freq_mis_out_crispr_range+'\t'+str(corrected_tot_freq_mis_out_crispr_range)+'\t'+corrected_count_mis_out_crispr_range+'\t'+str(corrected_tot_count_mis_out_crispr_range)+'\n')
                    mismatch_out_crispr_range = ''
                    mismatch_tot_freq_out_crispr_range = 0
                    mismatch_freq_out_crispr_range = ''
                    count_mismatch_out_crispr_range = ''
                    total_count_mis_out_crispr_range = 0
                    corrected_freq_mis_out_crispr_range = ''
                    corrected_tot_freq_mis_out_crispr_range = 0
                    corrected_count_mis_out_crispr_range = ''
                    corrected_tot_count_mis_out_crispr_range = 0
                    list_length_mis_out_crispr_range = ''

                if insertion_in_crispr_range :
                    if nb_ins_in_crispr_range > 1 :
                        yes_file.write(str(keys-len_primer_forward)+'\t<'+str(keys-len_primer_forward+len_ins_in_crispr_range-1)+'\t'+list_length_ins_in_crispr_range+'\t<'+str(len_ins_in_crispr_range)+'\tinsertion\t/\t'+str(insertion_in_crispr_range)+'\t'+str(insertion_freq_in_crispr_range)+'\t'+str(insertion_tot_freq_in_crispr_range)+'\t'+count_ins_in_crispr_range+'\t'+str(total_count_ins_in_crispr_range)+'\t'+corrected_freq_ins_in_crispr_range+'\t'+str(corrected_tot_freq_ins_in_crispr_range)+'\t'+corrected_count_ins_in_crispr_range+'\t'+str(corrected_tot_count_ins_in_crispr_range)+'\n')
                    else :
                        yes_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len_ins_in_crispr_range-1)+'\t'+list_length_ins_in_crispr_range+'\t'+str(len_ins_in_crispr_range)+'\tinsertion\t/\t'+str(insertion_in_crispr_range)+'\t'+str(insertion_freq_in_crispr_range)+'\t'+str(insertion_tot_freq_in_crispr_range)+'\t'+count_ins_in_crispr_range+'\t'+str(total_count_ins_in_crispr_range)+'\t'+corrected_freq_ins_in_crispr_range+'\t'+str(corrected_tot_freq_ins_in_crispr_range)+'\t'+corrected_count_ins_in_crispr_range+'\t'+str(corrected_tot_count_ins_in_crispr_range)+'\n')
                    insertion_in_crispr_range = ''
                    insertion_tot_freq_in_crispr_range = 0
                    insertion_freq_in_crispr_range = ''
                    count_ins_in_crispr_range = ''
                    total_count_ins_in_crispr_range = 0
                    nb_ins_in_crispr_range = 0
                    len_ins_in_crispr_range = 0
                    corrected_freq_ins_in_crispr_range = ''
                    corrected_tot_freq_ins_in_crispr_range = 0
                    corrected_count_ins_in_crispr_range = ''
                    corrected_tot_count_ins_in_crispr_range = 0
                    list_length_ins_in_crispr_range = ''
                
                if insertion_out_crispr_range :
                    if nb_ins_out_crispr_range > 1 :
                        no_file.write(str(keys-len_primer_forward)+'\t<'+str(keys-len_primer_forward+len_ins_out_crispr_range-1)+'\t'+list_length_ins_out_crispr_range+'\t<'+str(len_ins_out_crispr_range)+'\tinsertion\t/\t'+str(insertion_out_crispr_range)+'\t'+str(insertion_freq_out_crispr_range)+'\t'+str(insertion_tot_freq_out_crispr_range)+'\t'+count_ins_out_crispr_range+'\t'+str(total_count_ins_out_crispr_range)+'\t'+corrected_freq_ins_out_crispr_range+'\t'+str(corrected_tot_freq_ins_out_crispr_range)+'\t'+corrected_count_ins_out_crispr_range+'\t'+str(corrected_tot_count_ins_out_crispr_range)+'\n')
                    else :
                        no_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len_ins_out_crispr_range-1)+'\t'+list_length_ins_out_crispr_range+'\t'+str(len_ins_out_crispr_range)+'\tinsertion\t/\t'+str(insertion_out_crispr_range)+'\t'+str(insertion_freq_out_crispr_range)+'\t'+str(insertion_tot_freq_out_crispr_range)+'\t'+count_ins_out_crispr_range+'\t'+str(total_count_ins_out_crispr_range)+'\t'+corrected_freq_ins_out_crispr_range+'\t'+str(corrected_tot_freq_ins_out_crispr_range)+'\t'+corrected_count_ins_out_crispr_range+'\t'+str(corrected_tot_count_ins_out_crispr_range)+'\n')
                    insertion_out_crispr_range = ''
                    insertion_tot_freq_out_crispr_range = 0
                    insertion_freq_out_crispr_range = ''
                    count_ins_out_crispr_range = ''
                    total_count_ins_out_crispr_range = 0
                    nb_ins_out_crispr_range = 0
                    len_ins_out_crispr_range = 0
                    corrected_freq_ins_out_crispr_range = ''
                    corrected_tot_freq_ins_out_crispr_range = 0
                    corrected_count_ins_out_crispr_range = ''
                    corrected_tot_count_ins_out_crispr_range = 0
                    list_length_ins_out_crispr_range = ''


                if deletion_in_crispr_range :
                    if nb_del_in_crispr_range > 1 :
                        yes_file.write(str(keys-len_primer_forward)+'\t<'+str(keys-len_primer_forward+len_del_in_crispr_range-1)+'\t'+list_length_del_in_crispr_range+'\t<'+str(len_del_in_crispr_range)+'\tdeletion\t'+str(deletion_in_crispr_range)+'\t/\t'+str(deletion_freq_in_crispr_range)+'\t'+str(deletion_tot_freq_in_crispr_range)+'\t'+count_del_in_crispr_range+'\t'+str(total_count_del_in_crispr_range)+'\t'+corrected_freq_del_in_crispr_range+'\t'+str(corrected_tot_freq_del_in_crispr_range)+'\t'+corrected_count_del_in_crispr_range+'\t'+str(corrected_tot_count_del_in_crispr_range)+'\n')
                    else :
                        yes_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len_del_in_crispr_range-1)+'\t'+list_length_del_in_crispr_range+'\t'+str(len_del_in_crispr_range)+'\tdeletion\t'+str(deletion_in_crispr_range)+'\t/\t'+str(deletion_freq_in_crispr_range)+'\t'+str(deletion_tot_freq_in_crispr_range)+'\t'+count_del_in_crispr_range+'\t'+str(total_count_del_in_crispr_range)+'\t'+corrected_freq_del_in_crispr_range+'\t'+str(corrected_tot_freq_del_in_crispr_range)+'\t'+corrected_count_del_in_crispr_range+'\t'+str(corrected_tot_count_del_in_crispr_range)+'\n')
                    deletion_in_crispr_range = ''
                    deletion_tot_freq_in_crispr_range = 0
                    deletion_freq_in_crispr_range = ''
                    count_del_in_crispr_range = ''
                    total_count_del_in_crispr_range = 0
                    nb_del_in_crispr_range = 0
                    len_del_in_crispr_range = 0
                    corrected_freq_del_in_crispr_range = ''
                    corrected_tot_freq_del_in_crispr_range = 0
                    corrected_count_del_in_crispr_range = ''
                    corrected_tot_count_del_in_crispr_range = 0
                    list_length_del_in_crispr_range = ''
                
                if deletion_out_crispr_range :
                    if nb_del_out_crispr_range > 1 :
                        no_file.write(str(keys-len_primer_forward)+'\t<'+str(keys-len_primer_forward+len_del_out_crispr_range-1)+'\t'+list_length_del_out_crispr_range+'\t<'+str(len_del_out_crispr_range)+'\tdeletion\t'+str(deletion_out_crispr_range)+'\t/\t'+str(deletion_freq_out_crispr_range)+'\t'+str(deletion_tot_freq_out_crispr_range)+'\t'+count_del_out_crispr_range+'\t'+str(total_count_del_out_crispr_range)+'\t'+corrected_freq_del_out_crispr_range+'\t'+str(corrected_tot_freq_del_out_crispr_range)+'\t'+corrected_count_del_out_crispr_range+'\t'+str(corrected_tot_count_del_out_crispr_range)+'\n')
                    else :
                        no_file.write(str(keys-len_primer_forward)+'\t'+str(keys-len_primer_forward+len_del_out_crispr_range-1)+'\t'+list_length_del_out_crispr_range+'\t'+str(len_del_out_crispr_range)+'\tdeletion\t'+str(deletion_out_crispr_range)+'\t/\t'+str(deletion_freq_out_crispr_range)+'\t'+str(deletion_tot_freq_out_crispr_range)+'\t'+count_del_out_crispr_range+'\t'+str(total_count_del_out_crispr_range)+'\t'+corrected_freq_del_out_crispr_range+'\t'+str(corrected_tot_freq_del_out_crispr_range)+'\t'+corrected_count_del_out_crispr_range+'\t'+str(corrected_tot_count_del_out_crispr_range)+'\n')
                    deletion_out_crispr_range = ''
                    deletion_tot_freq_out_crispr_range = 0
                    deletion_freq_out_crispr_range = ''
                    count_del_out_crispr_range = ''
                    total_count_del_out_crispr_range = 0
                    nb_del_out_crispr_range = 0
                    len_del_out_crispr_range = 0
                    corrected_freq_del_out_crispr_range = ''
                    corrected_tot_freq_del_out_crispr_range = 0
                    corrected_count_del_out_crispr_range = ''
                    corrected_tot_count_del_out_crispr_range = 0
                    list_length_del_out_crispr_range = ''
