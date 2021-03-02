#!/usr/bin/python

import sys 
import re

#fileInput = sys.argv[1]
#fileOutput = sys.argv[2]

def write_file(input_file,output_file):
    with open(input_file ,'r') as inputF, open(output_file, 'w') as output :
        count = 0
        for line in inputF:
            line = line.rstrip()
            if not line.startswith('Sequence'):
                element_line = line.split('\t')
                count += 1
                output.write('>'+str(count)+'_count_'+str(element_line[3])+'_countCorrected_'+str(element_line[2])+'\n')
                output.write(element_line[0]+'\n')

#write_file(fileInput,fileOutput)
