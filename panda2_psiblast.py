import os,sys,time
import pandas as pd
import numpy as np

# Run psiblast and parser
def df_psiblast(fasta_dir, df):
    blast_path = "" # /blast-2.2.23/bin/blastpgp
    /home/czhao/tools/blast-2.2.23/bin/blastpgp
    if fasta_dir[-1] == '/': fasta_dir = fasta_dir[:-1]
    psiblast_dir = f'{fasta_dir}.2/psiblast_run/'
    results = []
    if not os.path.exists(psiblast_dir):
        os.system(f'mkdir -p {psiblast_dir}')
    with open(f'{fasta_dir}.2/psiblast_cmd.sh','w') as w_fh:
        for prot,seq in zip(df.proteins,df.sequences):
            fa_f = f'{psiblast_dir}/{prot}.fa'
            with open(fa_f,'w') as fa_fh:
                fa_fh.write(f'>{prot}\n{seq}\n')
            if not os.path.exists(fa_f+'2'):
                w_fh.write(f'{blast_path} -i {fa_f} -d data-cafa/train_data.fa -j 2 -b 1000 -o {fa_f}2\n')
            if not os.path.exists(fa_f+'3'):
                w_fh.write(f'perl parse_blast_result_2.pl {fa_f}2 {fa_f}3 2\n')
            results.append(f'{fa_f}3')
    os.system(f'sh {fasta_dir}.2/psiblast_cmd.sh')
    df['fa3'] = results
    return df
