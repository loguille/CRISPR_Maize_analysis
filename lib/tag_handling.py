#!/usr/bin/python

import sys 
import re

#Program to handle tag which are the first and the last 5 nt of the sequence, create a dict that get the different sequence of tag and their numbers

def handle_tag(input_file, output_file) :
    with open(input_file,'r') as qcfile, open(output_file,'w') as tagfile:
        dico_seq_tag = {}
        for line in qcfile:
            line=line.rstrip()
            element = line.split(':')
            tag = element[0][:5]+'+'+element[0][-5:]
            cut_seq = element[0][5:-5]
            if cut_seq not in dico_seq_tag :
                dico_seq_tag[cut_seq] = [[tag,element[1]]]
            else :
                dico_seq_tag[cut_seq].append([tag,element[1]])
        tagfile.write('Sequence\tTag and count\tCorrect count by tag\tTotal count\n')            
        for keys,elem in dico_seq_tag.items() :
            construct_elem = ''
            number = 0
            tot_count = 0
            for i in range(len(elem)) :
                number += 1
                tot_count += int(elem[i][1])
                construct_elem += elem[i][0]+' '+str(elem[i][1])+','    
            tagfile.write(keys+'\t'+construct_elem+'\t'+str(number)+'\t'+str(tot_count)+'\n')
    
            
