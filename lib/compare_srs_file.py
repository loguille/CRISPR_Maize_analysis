#!/usr/bin/python

import sys
import collections

#this program will be run in order too compare alignment on different target sequence in srsfile, different srsfile will be create containing the best alignment based on the score

def readfiles(dico_ali,namefile,new_name_file):
    with open(namefile,'r') as first_ali :
        count = 0
        for line in first_ali :
            line = line.rstrip()
            if line.startswith('# 1:'):
                ref_seq = line
            if line.startswith('# 2:'):
                name_seq = line
            if line.startswith('# Score:'):
                score_line = line 
                score = float(line.split(' ')[2])
            if line.startswith('# Length:'):
                length = line
            if line.startswith('# Identity:'):
                identity= line
            if line.startswith('# Similarity:'):
                similarity = line
            if line.startswith('# Gaps:'):
                gaps =line
            if not line.startswith('#'):
                element = [line]
                if element != ['']:
                    if count == 0 :
                        aln_tar = element
                        count += 1
                    elif count == 1:
                        aln =element
                        count += 1
                    else :
                        aln_seq =element
                        count = 0
                        if name_seq not in dico_ali :
                            dico_ali[name_seq] = [[ref_seq,score_line,score,length,identity,similarity, gaps, aln_tar,aln,aln_seq,new_name_file]]
                        else :
                            dico_ali[name_seq].append([ref_seq,score_line,score,length,identity,similarity, gaps, aln_tar,aln,aln_seq,new_name_file])
    new_dict = collections.OrderedDict(sorted(dico_ali.items()))
    return(new_dict)

def compare_score(dico):
    for i in dico.keys():
        scoremax = 0
        for element in range(len(dico[i])):
            if dico[i][element][2] > scoremax :
                scoremax = dico[i][element][2]
                inte = element
        with open(dico[i][inte][10],'a') as fileout :
            fileout.write('#=======================================\n')
            fileout.write('#\n')
            fileout.write('# Aligned_sequences: 2\n')
            fileout.write(dico[i][inte][0]+'\n')
            fileout.write(i+'\n')
            fileout.write('# Matrix: EDNAFULL83\n')
            fileout.write('# Gap_penalty: 10.0\n')
            fileout.write('# Extend_penalty: 0.5\n')
            fileout.write('#\n')
            fileout.write(dico[i][inte][3]+'\n')
            fileout.write(dico[i][inte][4]+'\n')
            fileout.write(dico[i][inte][5]+'\n')
            fileout.write(dico[i][inte][6]+'\n')
            fileout.write(dico[i][inte][1]+'\n')
            fileout.write('#\n#\n')
            fileout.write('#=======================================\n')
            fileout.write('\n')
            fileout.write(dico[i][inte][7][0]+'\n')
            fileout.write(dico[i][inte][8][0]+'\n')
            fileout.write(dico[i][inte][9][0]+'\n')
            fileout.write('\n\n')