#!/usr/bin/python

import sys
import re

srsFile=sys.argv[1]



class ParseNeedle():
    def __init__(self,name_seq, name_seq_target,length,number_of_seq,score,align_ref,align,align_seq,correct_nb_seq):
        self.name_seq = name_seq
        self.name_seq_target = name_seq_target
        self.length = length
        self.number_of_seq = number_of_seq
        self.score = score
        self.align_ref = align_ref
        self.align = align
        self.align_seq = align_seq
        self.correct_nb_seq = correct_nb_seq
        self.new_score = 0


def calcul_new_score(read): 
    for ali in read :
        aln = read[ali].align
        indel = ''
        score = 0
        for i in range(len(aln)) :
            if aln[i] == '.' :
                if indel :
                    score -= 5
                    indel = ''
                score -= 5
            elif aln[i] == '|':
                if indel :
                    score -= 5
                    indel = ''
                score += 4
            else :
                indel += 'i'
        score = int(score / float(read[ali].length) * 100)
        read[ali].new_score = score
    return(read)





def parseFilesrs(fileIN):
    read={}
    nb_seq=0
    total_nb_corrected = 0
    with open(fileIN ,'r') as srsfile:
        count = 1
        for line in srsfile:
            line = line.rstrip()
            if line.startswith('# 1:'):
                target_name = line.split(' ')[2]
            if line.startswith('# 2:'):
                sequence_name = line.split(' ')[2]
                number_of_sequence = line.split('_')[2]
                correct_numb_of_seq = line.split('_')[4]
                total_nb_corrected += int(correct_numb_of_seq)
                nb_seq += int(number_of_sequence)
            if line.startswith('# Length'):
                sequence_length = int(line.split(' ')[2])
            if line.startswith('# Score'):
                score_aln = float(line.split(' ')[2])
            if not line.startswith('#'):
                element = [line]
                if element != ['']:
                    targetSeq = element[0][21:21+sequence_length]
                    if count == 1:
                        aln_target = targetSeq
                        count += 1 
                    elif count == 2 :
                        aln = targetSeq
                        if len(aln) != len(aln_target):
                            i = len(aln_target) - len(aln)
                            while i != 0:
                                aln += ' '
                                i -= 1
                        count += 1
                    else :
                        aln_sequence = targetSeq
                        read[sequence_name] = ParseNeedle(sequence_name,target_name,sequence_length,number_of_sequence,score_aln,aln_target,aln,aln_sequence,correct_numb_of_seq)
                        count = 1
    calcul_new_score(read)
    return(read, nb_seq, total_nb_corrected)

#read,seq_nb = parseFilesrs(srsFile)
#for element in read :
#    print(read[element].new_score)                    