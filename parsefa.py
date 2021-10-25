#!/usr/bin/env python

################################################################
## This program will take fasta file (two- or multi-line) and ##
##              produce a two-line fasta file                 ##
################################################################

import argparse
def get_args():
    parser = argparse.ArgumentParser("a program to produce two-line fastas")
    parser.add_argument("-f", "--file", type=str, help="fasta file to pass through program", required=True)
    parser.add_argument("-o", "--output", type=str, help="name of output", required=True)

    return parser.parse_args()
args = get_args()

###assign arguments to variables inside of program
fasta = args.file
out = args.output

###open output file

outfa = open(out, "w")
sequence = ''

###read file
with open (fasta, "r") as fh:
    for line in fh:
        if line[0] == '>':
            if sequence != '':
                #print(len(sequence))
                outfa.write(sequence+'\n')
            header = line.strip()
            outfa.write(str(header)+'\n')
            sequence = ''
        else:
            seq = line.strip()
            sequence += seq
    #print(len(sequence))
    outfa.write(sequence+'\n')


outfa.close()

