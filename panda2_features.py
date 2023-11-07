import os
import sys
import glob
import numpy as np
import pandas as pd
from panda2_psiblast import df_psiblast
from panda2_psitopiden import *
from panda2_paac import df_paac
from panda2_diamond import df_diamond
from panda2_esm import df_esm
###############################################
#  
###############################################
fasta_dir = sys.argv[1]
if fasta_dir[-1] == '/': fasta_dir = fasta_dir[:-1]
feature_dir = fasta_dir+'.2'
proteins, sequences = [], []
for fasta_f in glob.glob(fasta_dir+'/*fasta'): # read in random order
    protein, sequence = "", ""
    for line in open(fasta_f):
        if '>' in line:
            protein = line[1:].rstrip()
        else:
            sequence+=line.rstrip()
    proteins.append(protein)
    sequences.append(sequence)
## df proteins + sequences
dictionary = {'proteins': proteins, 'sequences': sequences}
df = pd.DataFrame(data=dictionary)
## df run psiblast + fa3
df = df_psiblast(fasta_dir,df)
## df + psiblast top 10
## df + iden
df = df_psitopiden(df)
## df + diamondscore
df = df_diamond(fasta_dir+'.fasta',df)
## df + paac
df = df_paac(df)
## df + esmseq
df = df_esm(df)
## save pkl
#print(df)
df.to_pickle(f'{feature_dir}/panda2_features.pkl')
