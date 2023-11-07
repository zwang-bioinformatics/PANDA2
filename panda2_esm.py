import os,sys,time
import pandas as pd
import numpy as np
import torch
import esm
# load esm 
model, alphabet = esm.pretrained.esm1_t34_670M_UR50S()
batch_converter = alphabet.get_batch_converter()

def df_esm(df):
    esmseq = []
    for i, row in enumerate(df.itertuples()):
        data = [(row.proteins,row.sequences)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[34])
        token_embeddings = results["representations"][34]
        esmseq.append(token_embeddings[0].mean(0).numpy()) # best use [1:]
    df['esmseq'] = esmseq
    return df
