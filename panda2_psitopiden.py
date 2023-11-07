import pandas as pd
import pickle
import numpy as np
from scipy import sparse
import math
import sys
import networkx as nx
import numpy as np

sys.path.append("/home/czhao/tools/pygosemsim/")
import pygosemsim
from pygosemsim import similarity
from pygosemsim import graph
from pygosemsim import download
from pygosemsim import annotation
import functools
from pygosemsim import term_set

def df_psitopiden(df):
    # load materail
    train_f = 'data-cafa/train_data.pkl'
    df_train = pd.read_pickle(train_f)

    G3 = graph.from_resource("go.cafa3",'')
    terms_file = 'data-cafa/terms.pkl'
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    # create protein object
    prots = {}
    for prot, annotations in zip(df_train.proteins,df_train.annotations):
        Prot = Protein(prot,terms,annotations)
        prots[prot] = Prot
    for prot, fa3 in zip(df.proteins,df.fa3):
        Prot = Protein(prot,terms)
        Prot.psiblast_slim = fa3
        prots[prot] = Prot
    # collect features
    psitops, idens = [],[]
    for i, row in enumerate(df.itertuples()):
        psitops.append(sparse.coo_matrix(prots[row.proteins].blasttop(prots).T))
        idens.append(sparse.coo_matrix(np.expand_dims(prots[row.proteins].identity(prots,G3,terms_dict),axis=1)))
    df['psitops'] = psitops
    df['idens'] = idens
    return df

class Protein:
    def __init__(self,prot,classmap,annotations = []):
        self.prot = prot
        self.classmap = classmap
        self.psiblast_slim = ''
        self.prop = annotations
        self.target = [int(go in self.prop) for go in self.classmap]
        
    def psib_parser(self):
        evalue = {}
        identity = {}
        for line in open(self.psiblast_slim,'r'):
            terms = line.rstrip().split()
            if terms[0] == self.prot: continue
            if terms[1][0] == 'e': 
                 terms[1] = '1'+terms[1]
            if float(terms[1])>=1: continue
            evalue[terms[0]] = float(terms[1])
            identity[terms[0]] = float(terms[2])
        return evalue,identity
    
    def panda_score(self,evalue,base = 40):
        return min(-math.log10(float(evalue))/base,1)
    
    def blasttop(self,prots,top = 10,base = 40):
        top_targets = []
        evalue,identity = self.psib_parser()
        for counter,i in enumerate(evalue.keys()):
            if counter==10: break
            if evalue[i] == 0:
                score = 1
            else:
                score = self.panda_score(evalue[i],base)
            #print(i,evalue[i],score)
            top_targets.append(np.array(prots[i].target)*score)
        top_targets =  np.asarray(top_targets)
        top_shape = top_targets.shape
        #print(type(top_targets),np.sum(top_targets),top_targets,top_shape)
        if top_shape[0] == 0:
            top_targets = np.zeros((top,len(self.classmap)))
        else:
            top_padding = np.zeros((top-top_shape[0],top_shape[1]))
            top_targets = np.vstack((top_targets,top_padding))
        return top_targets

    def identity(self,prots,G,terms_dict):
        # from seq similarity
        evalue,identity = self.psib_parser()
        candi_gos = {}
        identity = {k: v for k, v in sorted(identity.items(), key=lambda item: item[1])}
        evalue = dict((k, v) for k, v in evalue.items() if v < 1)
        #print(identity)
        # similar seq proteins into candidates go terms
        for k in identity:
            if not k in evalue: continue
            if k == self.prot: continue
            if not k in prots: continue
            blast_prop = prots[k].prop
            for go in blast_prop:
                if not go in candi_gos:
                    candi_gos[go] = identity[k]
        # get propagation and occurances
        larger,occurs = self.dic_propagation_occur(G,candi_gos)
        # norm with 
        larger_nm = {}
        for go in larger:
            larger_nm[go] = larger[go]*(occurs[go]/(2*(1+occurs[go])) + 1/2)
        larger_nm_list = np.zeros(len(terms_dict))
        for go in larger_nm:
            if go in self.classmap:
                larger_nm_list[terms_dict[go]] = larger_nm[go]
        return larger_nm_list
    
    def identity2text(self,larger_nm):
        # dictionary to text
        text = ''
        for go in larger_nm:
            text+= f'{go}|{larger_nm[go]}\t'
        return text
    
    def dic_propagation_occur(self, G,goterms):
        larger = goterms.copy()
        occurs = {}
        for go in goterms:
            # check go term in network
            if not go in G:
                # check go term in alt ids
                if go in G.alt_ids:
                    go = G.alt_ids[go]
                    # check go alt term in network
                    if not go in G: continue
                else: continue
            # ancestor
            ancestor = set(annotation.propagation(G,go))
            for term in ancestor:
                if not term in larger:
                    larger[term] = 0
                # max of self score and children's score
                larger[term] = max(larger[term],goterms[go])
                # occurs counter
                if term in occurs:
                    occurs[term] += 1
                else:
                    occurs[term] = 1
            if go in occurs:
                occurs[go] += 1
            else:
                occurs[go] = 1
        return larger,occurs

