import sys,os
import pandas as pd
import numpy as np
from utils import Ontology
def df_diamond(fasta_f, df):
    # run diamond bit score
    feature_dir = fasta_f.replace('.fasta','.2')
    os.system(f'./diamond blastp -d data-cafa/train_data_train.dmnd -q {fasta_f} --outfmt 6 qseqid sseqid bitscore > {feature_dir}/diamond.res')
    # calculate diamond score
    go_rels = Ontology('data-cafa/go.obo', with_rels=True)
    train_data_file = 'data-cafa/train_data.pkl'
    train_df = pd.read_pickle(train_data_file)
    annotations = train_df['annotations'].values
    annotations = list(map(lambda x: set(x), annotations))
    terms_df = pd.read_pickle('data-cafa/terms.pkl')
    terms = terms_df['terms'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}

    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        prot_index[row.proteins] = i
    
    # BLAST Similarity (Diamond)
    diamond_scores = {}
    with open(f'{feature_dir}/diamond.res') as f:
        for line in f:
            it = line.strip().split()
            if it[0] not in diamond_scores:
                diamond_scores[it[0]] = {}
            diamond_scores[it[0]][it[1]] = float(it[2])

    blast_preds = []
    for i, row in enumerate(df.itertuples()):
        annots = {}
        annots_dia = np.zeros(len(terms))
        prot_id = row.proteins
        # BlastKNN
        if prot_id in diamond_scores:
            sim_prots = diamond_scores[prot_id]
            allgos = set()
            total_score = 0.0
            for p_id, score in sim_prots.items():
                allgos |= annotations[prot_index[p_id]]
                total_score += score
            allgos = list(sorted(allgos))
            sim = np.zeros(len(allgos), dtype=np.float32)
            for j, go_id in enumerate(allgos):
                s = 0.0
                for p_id, score in sim_prots.items():
                    if go_id in annotations[prot_index[p_id]]:
                        s += score
                sim[j] = s / total_score
            for go_id, score in zip(allgos, sim):
                annots[go_id] = score
                if go_id in terms_dict:
                    annots_dia[terms_dict[go_id]] = score 
        annots_dia = np.around(annots_dia, decimals=3)
        #print(prot_id,np.sum(annots_dia))
        blast_preds.append(annots_dia)
    df['diamond'] = blast_preds
    return df

