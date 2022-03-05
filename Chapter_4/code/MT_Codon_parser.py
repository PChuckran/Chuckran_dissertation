#! /usr/bin/env python

import pandas as pd
import numpy as np
import csv
import re
import sys

# Read codons frequency and GC
##create codon dictionary
input_file = sys.argv[1]
output_file_name = input_file.strip(".fna")+"_codon_frequency.csv"

CodonsDict = {
    "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
    "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
    "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
    "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
    "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
    "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
    "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
    "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
    "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
    "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
    "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
    "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
    "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
    "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
    "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
    "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}

NuclDict = {"A": 0, "C": 0, "G": 0, "T": 0}

colnames = ["ID", "contig"]
colnames.extend(CodonsDict.keys())
colnames.extend(NuclDict.keys())
# Create output file, write headers
output_file = open(output_file_name, "w")
wr = csv.writer(output_file, quoting=csv.QUOTE_ALL)
wr.writerow(colnames)

output_line = list()
seq = str()
ID = str()
contig_name = int()
row = 0
with open(input_file, 'r') as list_paths:
    # For each line, identify if it is a label or string of AA and parse accordingly
    for line in list_paths:
        if line.startswith(">"):
            # If it's the very first line, go straight to parse
            if row == 0:
                next
            else:
                for i in range(0, len(seq), 3):
                    codon = seq[i:i+3].upper()
                    if len(codon) == 3:
                        CodonsDict[codon] = CodonsDict[codon] + 1
                    else:
                        print(ID)
                        print(codon+" is not a codon")
                    for nuc in codon:
                        NuclDict[nuc] = NuclDict[nuc] + 1
                output_line.append(ID)
                output_line.append(contig_name)
                output_line.extend(list(CodonsDict.values()))
                output_line.extend(list(NuclDict.values()))        
                # Marks the end of a sequeence and the start of a new one
                # Append the last sequence to the DF
                wr.writerow(output_line)
                # Re-initialize dictionary so everything is back to zero
                seq = str()
                CodonsDict = dict.fromkeys(CodonsDict, 0)
                NuclDict = dict.fromkeys(NuclDict, 0)
                output_line = list()
            row += 1
            ID = line.strip("\n").replace(" # ", ";").split(';')[0]
            ID = ID.replace('>','')
            contig_name = int(re.search('^Ga[0-9]+_([0-9]+)_*', ID).group(1))
        elif line.startswith(tuple(['G', 'T', 'C', 'A', 'a', 't', 'c', 'g'])):
            seq = seq+line.strip("\n")
        else:
            print("space")
            print(line)
    # To parse the last sequence 
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3].upper()
        if len(codon) == 3:
            CodonsDict[codon] = CodonsDict[codon] + 1
        else:
            print(codon+" is not a codon")
        for nuc in codon:
            NuclDict[nuc] = NuclDict[nuc] + 1
    output_line.append(ID)
    output_line.append(contig_name)
    output_line.extend(list(CodonsDict.values()))
    output_line.extend(list(NuclDict.values()))
    wr.writerow(output_line)
  
output_file.close()
