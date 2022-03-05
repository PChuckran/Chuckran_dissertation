#! /usr/bin/env python

import pandas as pd
import numpy as np
import math
import time
import glob
import json

# create a dictionary of each AA and corresponding codon
aa_dict = {'A':["GCT", "GCC", "GCA", "GCG"],
           'C':["TGT", "TGC"],
           'D':["GAT", "GAC"],
           'E':["GAA", "GAG"],
           'F':["TTT", "TTC"],
           'G':["GGT", "GGC", "GGA", "GGG"],
           'H':["CAT", "CAC"],
           'I':["ATT", "ATC", "ATA"],
           'K':["AAA", "AAG"],
           'L':["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
           'M':["ATG"],
           'N':["AAT", "AAC"],
           'P':["CCT", "CCC", "CCA", "CCG"],
           'Q':["CAA", "CAG"],
           'R':["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
           'S':["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
           'T':["ACT", "ACC", "ACA", "ACG"],
           'V':["GTT", "GTC", "GTA", "GTG"],
           'W':["TGG"],
           'Y':["TAT", "TAC"]}

aa_dict_sub = {'A':["GCT", "GCC", "GCA", "GCG"],
           'C':["TGT", "TGC"],
           'D':["GAT", "GAC"],
           'E':["GAA", "GAG"],
           'F':["TTT", "TTC"],
           'G':["GGT", "GGC", "GGA", "GGG"],
           'H':["CAT", "CAC"],
           'I':["ATT", "ATC", "ATA"],
           'K':["AAA", "AAG"],
           'L':["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
           'N':["AAT", "AAC"],
           'P':["CCT", "CCC", "CCA", "CCG"],
           'Q':["CAA", "CAG"],
           'R':["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
           'S':["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
           'T':["ACT", "ACC", "ACA", "ACG"],
           'V':["GTT", "GTC", "GTA", "GTG"],
           'Y':["TAT", "TAC"]}

codon_list = ['TTT','TTC','TTA','TTG','CTT','CTC',
                       'CTA','CTG','ATT','ATC','ATA','ATG',
                       'GTT','GTC','GTA','GTG','TAT','TAC',
                       'TAA','TAG','CAT','CAC','CAA','CAG',
                       'AAT','AAC','AAA','AAG','GAT','GAC',
                       'GAA','GAG','TCT','TCC','TCA','TCG',
                       'CCT','CCC','CCA','CCG','ACT','ACC',
                       'ACA','ACG','GCT','GCC','GCA','GCG',
                       'TGT','TGC','TGA','TGG','CGT','CGC',
                       'CGA','CGG','AGT','AGC','AGA','AGG',
                       'GGT','GGC','GGA','GGG']

codons_to_exclude = ["TGG", "ATG", "TAG", "TAA", "TGA"]


def calculate_MILC(input_freq, base_freq):
    Msum = 0
    aa_included = []
    L = 0
    for aa in aa_dict:
        Ma = 0
        codons_to_get = aa_dict[aa]
        sum_codons = float(sum(input_freq[codons_to_get]))
        if sum_codons == 0:
            next
        else:
            aa_included.append(aa)
            for codon in codons_to_get:
                Mcodon = 0
                Oc = float(input_freq[codon])
                L = L + Oc
                if Oc == 0:
                    next
                else:
                    fc = Oc/sum_codons
                    #print("fc "+str(fc))
                    gc = float(base_freq[codon])
                    #print("gc "+str(gc))
                    #print("Oc " +str(Oc))
                    Mcodon = Oc * np.log(fc/gc)
                    Ma = Ma + Mcodon
            Msum = Msum + Ma
    rsum = 0 
    for aaM in aa_included:
        ra = len(aa_dict[aaM])
        rsum = rsum + (ra-1)
    Cor = (rsum/L)-0.5
    MILC = (Msum/L)-Cor
    return MILC

def get_wi(codon_wi, rscu_freq, co_ct):
    RSCU_co = rscu_freq[codon_wi]
    wi = [RSCU_co] * co_ct
    return wi

def get_CAI_and_FOP(taxa_id, input_freq, rscu_for_taxa, op_codon_list):
    #calculate fop
    CAI = float()
    op_ct = float()
    codon_to_subtract = float()
    codon_total = float()
    wi_ct = list()
    for codon in codon_list:
        codon_abundance = input_freq[codon]
        codon_total = codon_total + codon_abundance 
        if codon in codons_to_exclude:
            codon_to_subtract = codon_to_subtract + codon_abundance
        elif codon in op_codon_list:
            op_ct = op_ct + codon_abundance
            wi_ct = wi_ct + get_wi(codon, rscu_for_taxa, codon_abundance)
        else:
            wi_ct = wi_ct + get_wi(codon, rscu_for_taxa, codon_abundance)
    fop = op_ct/(codon_total-codon_to_subtract)
    CAI = math.exp(np.log(wi_ct).mean())
    if CAI == 0:
        print("""
        Error: CAI equals zero, indicated underflow problem
         - Setting to NA
        """)
        print(taxa_id)
        #print(wi_ct)
        CAI = "NaN"
    else:
        next
    return CAI, fop

phy_freq = open("/scratch/pfc25/transcript_quality/DNA/co_by_phy/phy_freq.json", "r").read()
phy_aa_max = open("/scratch/pfc25/transcript_quality/DNA/co_by_phy/phy_aa_max.json", "r").read()
phy_OP = open("/scratch/pfc25/transcript_quality/DNA/co_by_phy/phy_OP.json", "r").read()
taxa_freq = open("/scratch/pfc25/transcript_quality/DNA/co_by_taxa/taxa_freq.json", "r").read()
taxa_aa_max = open("/scratch/pfc25/transcript_quality/DNA/co_by_taxa/taxa_aa_max.json", "r").read()
moo = open("/scratch/pfc25/transcript_quality/DNA/co_by_taxa/taxa_OP.json", "r")
taxa_OP = moo.read()

phy_freq = json.loads(phy_freq)
phy_aa_max = json.loads(phy_aa_max)
phy_OP = json.loads(phy_OP)
taxa_freq = json.loads(taxa_freq)
taxa_aa_max = json.loads(taxa_aa_max)
taxa_OP = json.loads(taxa_OP)

sample_freq = pd.read_csv(glob.glob("*codon_frequency.csv")[0])

taxa_names = ["gene_ID", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage"]
#taxa_file_name = sys.argv[4]
taxa_file_name = glob.glob("*.a.phylodist.txt")[0]
sample_tax = pd.read_csv(taxa_file_name, delimiter="\t", header = None, names = taxa_names)
sample_tax["homolog_taxon_oid"] = sample_tax["homolog_taxon_oid"].astype("str").str.strip(".0")
sample_tax = sample_tax[sample_tax["percent_identity"] > 50]

sample_freq = sample_freq.merge(sample_tax, how = "left", on = "gene_ID")

sample_freq["phylum"] = sample_freq["lineage"].str.extract("^\w+;(\w+)").astype('str')

phylum_in_dict = phy_aa_max.keys()
taxa_in_dict = taxa_aa_max.keys()
cai_list = list()
fop_list = list()
milc_list = list()
pcai_list = list()
pfop_list = list()
pmilc_list = list()
for line in range(len(sample_freq)):
    #NA until proven otherwise
    CAI_G = "NaN"
    FOP_G = "NaN"
    milc_out = "NaN"
    pCAI_G = "NaN"
    pFOP_G = "NaN"
    pmilc_out = "NaN"
    input_line = sample_freq.loc[line]
    taxa_id_in = input_line["homolog_taxon_oid"]
    phylum_id_in = input_line["phylum"]
    if taxa_id_in in taxa_in_dict:
        op_codon_list_in = taxa_OP[taxa_id_in]
        rscu_for_taxa_in = taxa_aa_max[taxa_id_in]
        CAI_G, FOP_G =get_CAI_and_FOP(taxa_id_in, input_line, rscu_for_taxa_in, op_codon_list_in)
        milc_freq = taxa_freq[taxa_id_in]
        milc_out = calculate_MILC(input_line, milc_freq)
    else:
        next
    if phylum_id_in in phylum_in_dict:
        phy_op_codon_list_in = phy_OP[phylum_id_in]
        rscu_for_phy_in = phy_aa_max[phylum_id_in]
        pCAI_G, pFOP_G =get_CAI_and_FOP(phylum_id_in, input_line, rscu_for_phy_in, phy_op_codon_list_in)
        phy_milc_freq = phy_freq[phylum_id_in]
        pmilc_out = calculate_MILC(input_line, phy_milc_freq)
    else:
        next
    cai_list.append(CAI_G)
    fop_list.append(FOP_G)
    milc_list.append(milc_out)
    pcai_list.append(pCAI_G)
    pfop_list.append(pFOP_G)
    pmilc_list.append(pmilc_out)
    
sample_freq["CAI_t"]= cai_list
sample_freq["FOP_t"]= fop_list
sample_freq["MILC_t"]= milc_list
sample_freq["CAI_p"]= pcai_list
sample_freq["FOP_p"]= pfop_list
sample_freq["MILC_p"]= pmilc_list

base_name = taxa_file_name.strip('.a.phylodist.txt')
sample_freq.to_csv(base_name+"_indicies.csv")


