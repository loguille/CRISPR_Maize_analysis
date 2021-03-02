#!/usr/bin/env python3

import sys
import math


class base():
    def __init__(self, type, text, freq):
        self.type = type
        self.freq = freq
        self.text = text
        self.color = ""

def read_and_cut_primers(filePrimers, target_sequence) :
    with open(filePrimers) as primer : 
        bool_forward = False
        bool_reverse = False
        for line in primer :
            line = line.rstrip()
            if line.startswith('>') :
                if line.split('_')[1] == 'forward' :
                    bool_forward = True
                    continue
                elif line.split('_')[1] == 'reverse' :
                    bool_reverse = True
                    continue
            else :
                if bool_forward == True :
                    forward_primer = len(line)
                    bool_forward = False 
                if bool_reverse == True :
                    reverse_primer = len(line)
                    bool_reverse = False
    with open(target_sequence) as sequence :
        for ligne in sequence :
            ligne = ligne.rstrip()
            if not ligne.startswith('>'):
                sequence_wo_primers = ligne[forward_primer : -reverse_primer]
    return(sequence_wo_primers)


def codeSVG(x, y, SVGid, ncl, color="#000000", opacity=1):
    """
    function to write a SVG letter, x correspond to the x coordinate, y to the y coordinate, SVGid correspond to
    a string that contain different denomination like seqRef seqProb or seqMut and a number corresponding to the 
    position, and  ncl contain the nucleotide A,T,C,G or D.
    """
    if str(SVGid).isdigit():
        SVGid = "text" + str(SVGid)
    return '''<text
    xml:space="preserve"
    style = "font-style:normal;font-weight:normal;font-size:12px;line-height:1.25;font-family:Courier;letter-spacing:0px;word-pacing:0px;fill:''' + str(color) + ''';fill-opacity:''' + str(opacity) + ''';stroke:none;stroke-width:0.26458332"
    x = "''' + str(x) + '''"
    y = "''' + str(y) + '''"
    id="''' + str(SVGid) + '''">''' + str(ncl) + '''</text>'''


def limit_crispr_zone(pam_position,direction):
    dico_limit_crispr = {}
    for i in range(len(pam_position)) :
        if direction[i] == 'sens' or direction[i] == 'forward' :
            dico_limit_crispr[pam_position[i]-20] = pam_position[i]
        else :
            dico_limit_crispr[pam_position[i]+4] = pam_position[i]+24
    print(dico_limit_crispr)
    return(dico_limit_crispr)

def construct_logo(freq_file,primers_file,reference_sequence,cutOff,fileOut,base_crispr,direction):
    sequence_wo_primers = read_and_cut_primers(primers_file,reference_sequence)
    start_end_dico = limit_crispr_zone(base_crispr,direction)
    reads = {}  # list of object base
    maxPos = 0
    with open(freq_file,"r") as freqFile:
        for line in freqFile:
            line = line.rstrip()
            if line.startswith("Start"):#skip header
                continue
            column = line.split('\t')
            type = column[2]
            if type == 'mismatch' or type == 'insertion' :
                seq = column[4]
            else :
                seq = column[3]
            pos = int(column[0])-1
            freq = float(column[5])
            cutOff = float(cutOff)
            if freq > cutOff:
                if pos not in reads.keys():
                    reads[pos] = []
                    maxPos = pos
                reads[pos].append(base(type, seq, freq))


    if maxPos < len(sequence_wo_primers):
        #to print seqRef if count table is poor
        maxPos = len(sequence_wo_primers)

    ########################
    # Calculs matrice sequences #
    ########################
    matrix = [{}, {}]

    # reference sequence
    pos = 0
    for ncl in sequence_wo_primers:
        matrix[0][pos] = (ncl)
        pos += 1
    pos = 0
    # other reads
    while pos < maxPos + 1:
        if pos not in matrix[1].keys():
            # create if new position
            matrix[1][pos] = {"A": 0, "T": 0, "C": 0, "G": 0, "D": 0}
        if pos not in reads.keys():
            # perfect match : no mutation
            pos += 1
            continue
        for objBase in reads[pos]:
            # reads[pos] is a list of objet base
            type = objBase.type
            # treatement mutation by mutation
            if type == 'mismatch':
                matrix[1][pos][objBase.text] += objBase.freq
            elif type == 'deletion':
                # objBase.text = distance = gap number
                long = len(objBase.text)
                for i in range(long):
                    if pos + i not in matrix[1].keys():
                        matrix[1][pos + i] = {"A": 0, "T": 0, "C": 0, "G": 0, "D": 0}
                    matrix[1][pos + i]["D"] += objBase.freq
            elif type == 'insertion':
                long = len(objBase.text)
                for i in range(long):
                    if pos + i not in matrix[1].keys():
                        matrix[1][pos + i] = {"A": 0, "T": 0, "C": 0, "G": 0, "D": 0}
                    matrix[1][pos + i][objBase.text[i]] += objBase.freq #for insertion they add the frequency to the corresponding base
            else:
                pass
        pos += 1

    ########################
    # Calculs document properties#
    ########################
    width = maxPos * 8 + 12  #choose arbitrary
    height = 100


    ########################
    # Header and footer #
    ########################
    header = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>


    <svg
    xmlns:dc="http://purl.org/dc/elements/1.1/"
    xmlns:cc="http://creativecommons.org/ns#"
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    xmlns:svg="http://www.w3.org/2000/svg"
    xmlns="http://www.w3.org/2000/svg"
    xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
    height = "''' + str(height) + '''"
    width = "''' + str(width) + '''"
    viewBox = "0 0 ''' + str(width) + ''' ''' + str(height) + '''"
    version="1.1"
    id="svg1"
    script:version="0.1. (2020-04-01)"
    sodipodi:docname="logo.svg">
    <metadata
    id="metadata5">
    <rdf:RDF>
    <cc:Work
        rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
        rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
    </cc:Work>
    </rdf:RDF>
    </metadata>
    <g
    inkscape:label="Calque 1"
    inkscape:groupmode="layer"
    id="layer1">
    '''

    footer = '''
    </g>
    </svg>'''


    ########################
    # Writting #
    ########################
    x = 0
    y = [75, 65, 55, 45, 35]
    id = 10

    # max to min
    colors = ["#7800ff", "#7800ff", "#7800ff", "#7800ff", "#7800ff"]


    f = open(fileOut, "w")

    f.write(header)

    for pos in range(0, maxPos + 1):
        # seqRef
        if pos in matrix[0].keys():#here they construct the reference sequence
            ncl = matrix[0][pos]
            SVGid = "seqRef" + str(id)
            for start,end in start_end_dico.items():
                if pos >= start and pos < end :
                    f.write(codeSVG(x, 100, SVGid, ncl, "#ff0000", opacity=1))
                else :
                    f.write(codeSVG(x, 100, SVGid, ncl, "#4d4d4d", opacity=0.7))
        # other bases
        id += 1
        order = []  # order list for ncl
        sum = 0 # sum frequencies
        for ncl in matrix[1][pos]:
            if matrix[1][pos][ncl] != 0:
                # pass frequencie equal 0
                sum += matrix[1][pos][ncl]
                order.append((matrix[1][pos][ncl], ncl))
        orderSorted = sorted(order, key=lambda tup: tup[0], reverse=True)
        # seq most probably
        # frequencie of seqRef in pos is equal to 1 - sum
        Sid = "seqProb" + str(id)
        if len(orderSorted) > 0 and orderSorted[0][0] > (1 - sum):
            f.write(codeSVG(x, 86, SVGid=Sid, ncl=orderSorted[0][1], color="#37ABC8"))
        elif pos in matrix[0].keys():
            f.write(codeSVG(x, 89, SVGid=Sid, ncl=matrix[0][pos], color="#4d4d4d"))
        else:
            pass
        id += 1
        # seq muted
        for c, elementOrdered in enumerate(orderSorted):
            Sid = "seqMut" + str(id)
            freq = elementOrdered[0]
            ncleotid = elementOrdered[1]
            if ncleotid != 'D':
                f.write(codeSVG(x, y=y[c], SVGid=Sid,ncl=ncleotid, color=colors[c], opacity=math.log(freq*1000,10)))
            else :
                f.write(codeSVG(x, y=y[c], SVGid=Sid,ncl=ncleotid, color="#036b28", opacity=math.log(freq*1000,10)))
            id += 1
        x += 8
        id = int(id / 10 + 1) * 10  # for start next id to the next tenth
    f.write(footer)
    f.close()
