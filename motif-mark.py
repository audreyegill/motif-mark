#!/usr/bin/env python

######################
###                ###
###    Modules     ###
###                ###
######################

import regex as re
import cairo
import argparse
import os
from math import ceil


########################
###                  ###
###    Arguments     ###
###                  ###
########################

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="Path to input fasta file.", type=str, required=True)
    parser.add_argument("-o","--out", help="Name of output .svg file. Defaults to <input file name>.svg", type=str, required=False)
    parser.add_argument("-m","--motifs", help="Path to a text file consisting only of motifs which are separated by a new line.", type=str, required=False)
    return parser.parse_args()

args = arguments()
fasta_file = args.file
motif_file = args.motifs

if args.out is None:
    fn = file.split('/')[-1].split('.')[0]+'.svg'
else:
    fn = args.out
    
#############################
###                       ###
###   Parse everything    ###
###                       ###
#############################    

def parse_files(fasta_file, motif_file):
    ''' 
    Input: a fasta file path and a motif file path
        Note: the motif file must have a single motif per line
    Returns: seq_max = maximum sequence length
             gene = {fasta header: single-line sequence} dictionary
             exon = {fasta header: exon position} dictionary
             regex_dict = {motif: IUPAC regex match} dictionary'''
    
    ### Initialize first loop
    with open(fasta_file,"r") as fasta:
        read = ''
        first = True
        gene = {}
        head = False
        seq_max = 0
        # This loop results in {header:sequence} dictionary
        # and provides a maximum sequence length
        for line in fasta:
            if line.startswith('>'):
                if head:
                    gene[head] = read
                head = line.strip('>').strip()
                read = ""
            else:
                read += line.strip()
                if len(read) > seq_max:
                    seq_max = len(read)
            gene[head] = read
    
    ### Initialize second loop 
    exons = {}
    # This loop creates a {head:[exon positions]} dictionary
    for head in gene.keys():
        position = []
        for exon in re.finditer("[A-Z]+", gene[head]):
            position.append(exon.span())
        exons[head] = position
    
    ### Initialize third loop 
    regex_dict = {}
    # This loop results in a {motif: IUPAC regex match} dictionary 
    with open(motif_file) as motifs:
        iupac = {"A":"[Aa]", "T":"[TtUu]", "C":"[Cc]", "G":"[Gg]", 
            "U":"[TtUu]", "M":"[AaCc]", "R":"[AaGg]", "W":"[AaTtUu]",
            "S":"[CcGg]", "Y":"[CcTtUu]", "K":"[GgTtUu]", 
            "V":"[AaCcGg]", "H":"[AaCcTtUu]", "D":"[AaGgTtUu]", 
            "B":"[CcGgTtUu]", "N":"[AaTtCcGgUu]"}
        for line in motifs:
            regex_motif = ''
            line = line.strip().upper()
            for char in line:
                regex_motif += iupac[char]
            regex_dict[line] = regex_motif

            
    return seq_max, gene, exons, regex_dict

seq_max, gene, exons, regex_dict = parse_files(fasta_file, motif_file)



###############################
###                         ###
###    Set up cairo bits    ###
###                         ###
###############################

width, height = ceil((seq_max + 00)/100)*100, 300*len(regex_dict)
count = 1 # gene count

surface = cairo.SVGSurface(fn, width, height)
context = cairo.Context(surface)
context.set_font_size(25)

colors = [[0.008, 0.93, 0.56], # sea foam?
          [0.93, 0.81, 0.008], # gold
          [0.24, 0.62, 1.0], # sky blue
          [0.92, 0.18, 0.51], # fushia?
          [0.13, 0.65, 0.13], # greenish
          [0.93,0.10,0.10], # red
          [0.17, 0.31, 1.0], # blue
          [0.18,0.92,0.91], # teal?
          [0.95, 0.58, 0.06], # orange
          [0.64, 0.93, 0.008], # lime
          [0.48, 0.07, 0.68], # grape
          [0.6, 0.125, 0.125], # brick?
          [1, 1, 0] # just straight up yellow 
         ]

############################
###                      ###
###    Error messages    ###
###                      ###
############################

#if condition:
#    sys.exit("Message")

if len(regex_dict) > len(colors):
	exit('Sorry, you have too many motifs :( only 13 max motifs supported')

if '.fa' not in fasta_file:
    exit('Sorry I can only read fastas :(')
	
if ".txt" not in motif_file:
    exit('Sorry I only accept motifs in .txt format :(')

#############################
###                       ###
###    Get to plotting    ###
###                       ###
#############################


for head, seq in gene.items():
    # update y coords for each gene
    y = 250 + 250*(count - 1)
    
    # Gene line
    context.set_source_rgb(0,0,0)
    context.set_line_width(10)
    context.move_to(100, y)
    context.line_to(100 + len(seq), y)
    context.stroke()
    
    # Add exon
    start = exons[head][0][0]
    end = exons[head][0][1]
    context.set_line_width(50)
    context.move_to(100 + start, y)
    context.line_to(100 + end, y)
    context.stroke()
    
    # Add motif matches
    motif_pos = {}
    for motif in regex_dict:
        motif_pos[motif] = []
        for match in re.finditer(regex_dict[motif], seq):
            motif_pos[motif].append(match.span())
    
    # initialize motif count (number motifs - nummo)
    nummo = 0
    
    # Plot motif matches
    for motif in motif_pos.keys():
        mocols = colors[nummo]
        
        # plot motifs on gene
        context.set_source_rgb(mocols[0],mocols[1], mocols[2])
        pos = motif_pos[motif]
        for tuple in pos:
            start = tuple[0]
            end = tuple[1]
            context.set_line_width(50)
            context.move_to(100 + start, y)
            context.line_to(100 + end, y)
            context.stroke()
        nummo += 1
    
    # Plot gene name
    context.set_source_rgb(0,0,0)
    context.move_to(100, y - 50)
    context.show_text(head)
    context.select_font_face("")
    
    count += 1


# Print the motif names in their colors at the top of the image
redo_nummo = 0
name_len = 0
# The above for loop did not like having this as part of its functionality
# it's fine here
for motif in motif_pos.keys():
    mocols = colors[redo_nummo]
    context.set_source_rgb(mocols[0],mocols[1], mocols[2])
    context.move_to(100*(1+redo_nummo)+name_len*5, 100)
    context.show_text(motif)
    redo_nummo += 1
    name_len += len(motif)

# Add text for the exon
context.set_source_rgb(0,0,0)
context.move_to(100*(2+redo_nummo)+name_len*4, 100)
context.show_text("Exon")
 
surface.finish()

 


