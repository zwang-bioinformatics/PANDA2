import numpy as np
import pandas as pd
import sys,os
import datetime
import pickle

from utils import Ontology, FUNC_DICT

import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.nn import Sequential as Seq, Linear as Lin, ReLU
from torch_geometric.nn import MetaLayer
import torch.nn.functional as F
from torch_scatter import scatter_mean
from math import sqrt

model_tag = 'gcn37onlyesm_cafa'
ont = 'bcm'
data_dir = '/home/czhao/Web_Server/PANDA2_server/data-cafa'
go_file = f'{data_dir}/go.obo'
terms_file = f'{data_dir}/terms.pkl'
edge_file = f'{data_dir}/train_edges_{ont}.pkl'

in_file = sys.argv[1]
out_file = in_file.replace('panda2_features.pkl','panda2_prediction.txt')

load_model = True


def main():
    # load GO
    go_rels = Ontology(f'{data_dir}/go.obo', with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    nb_classes = len(terms)
    #print('GO terms loaded',datetime.datetime.now())
    # load edge, edfe feat
    with open(edge_file,'rb') as fh:
        edges, edges_feat = pickle.load(fh)
    edges = np.asarray(edges)
    edges_feat = np.asarray(edges_feat)
    edges_feat = np.full(edges_feat.shape,1)
    edges = np.concatenate((np.array(edges[:,:2]), np.array([edges[:,1],edges[:,0]]).transpose()))
    edges_feat = np.concatenate((edges_feat, edges_feat))
    #print('GO Graph loaded',datetime.datetime.now())
    # load pickle
    pred_df = pd.read_pickle(in_file)
    pred_generator = DFGenerator(pred_df, terms_dict, nb_classes, 1,edges,edges_feat)
    #print('Train+Test generator ready',datetime.datetime.now())
    # load model
    device = torch.device("cpu")
    model_gcn = GNet(nb_classes,0.).to(device)
    model_cnn = Deepgo(nb_classes).to(device)

    if load_model:
        model_gcn.load_state_dict(torch.load(f"{model_tag}.gcn"))
        model_cnn.load_state_dict(torch.load(f"{model_tag}.cnn"))
        predictions = predicate(model_cnn, model_gcn, device, pred_generator, pred_df)
        predictions.to_pickle(out_file.replace('txt','pkl'))
        prediction2text(terms, predictions, out_file)
        print(datetime.datetime.now(),'Prediction saved to',out_file)
        
def predicate(model_cnn, model_gcn, device, pred_generator, pred_df):
    model_cnn.eval()
    model_gcn.eval()
    preds = []
    pred_generator.reset()
    with torch.no_grad():
        while pred_generator.start < pred_generator.size:
            data,u,gfeat = pred_generator.next()
            data = Variable(torch.from_numpy(data),requires_grad=False).to(device, dtype=torch.float)
            u = Variable(torch.from_numpy(u),requires_grad=False).to(device, dtype=torch.float)
            x, edges, edges_feat, gbatch = toOnegraph(gfeat)
            x = Variable(torch.from_numpy(x)).to(device, dtype=torch.float)
            edges = Variable(torch.from_numpy(edges.transpose())).to(device, dtype=torch.long)
            edges_feat = Variable(torch.from_numpy(edges_feat)).to(device, dtype=torch.float)
            gbatch = Variable(torch.from_numpy(gbatch)).to(device, dtype=torch.long)
            pred_cnn = model_cnn(data,u)
            pred = torch.sigmoid(model_gcn(pred_cnn,x,u, edges, edges_feat, gbatch)).cpu().detach().numpy()
            preds.append(pred[0,0])
    pred_df['preds'] = list(preds)
    return pred_df

def prediction2text(terms, predictions, prediction_format_f):
    prediction_format_fh = open(prediction_format_f, 'w')
    prediction_format_fh.write(f'AUTHOR PANDA2\nMODEL 1\nKEYWORDS graph network, sequence alignment.\n')
    for i, row in enumerate(predictions.itertuples()):
        target = row.proteins
        scores = row.preds
        # sort scores
        scores_dic = {} 
        for i,score in enumerate(scores):
            score = round(score,2)
            if score > 1 or score<= 0.09:
                continue
            scores_dic[str(terms[i])] = score
        scores_sorted = dict(sorted(scores_dic.items(), key=lambda item: item[1],reverse=True))
        # save scores
        for key in scores_sorted:
            prediction_format_fh.write(f'{target}\t{key}\t{"%.2f" % scores_sorted[key]}\n')
    prediction_format_fh.write('END\n')
    prediction_format_fh.close()

class DFGenerator(object):

    def __init__(self, df, terms_dict, nb_classes, batch_size,edges,edges_feat):
        self.start = 0
        self.size = len(df)
        self.df = df
        self.batch_size = batch_size
        self.nb_classes = nb_classes
        self.terms_dict = terms_dict
        self.edges = torch.from_numpy(edges)
        self.edges_feat = torch.from_numpy(edges_feat)

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def reset(self):
        self.start = 0

    def next(self):
        if self.start < self.size:
            batch_index = np.arange(
                self.start, min(self.size, self.start + self.batch_size))
            df = self.df.iloc[batch_index]
            data_onehot = np.zeros((len(df), 4, 2000), dtype=np.int32)
            u = np.zeros((len(df),1,21+1280))
            data_feat = []
            for i, row in enumerate(df.itertuples()):
                seq = row.sequences
                u[i,0,1:21] = row.paac/100
                u[i,0,21:] = row.esmseq
                nodefeat = row.psitops #torch.from_numpy(row.node_feat)
                nodefeat = nodefeat.todense()
                dia = np.expand_dims(np.asarray(row.diamond),axis=1)
                iden = (row.idens).todense()
                nodefeat = np.concatenate((nodefeat,dia,iden),axis=1)
                data_feat.append([nodefeat,self.edges,self.edges_feat])
            self.start += self.batch_size
            return (data_onehot,u,data_feat)
        else:
            self.reset()
            return self.next()

def toOnegraph(gfeat):
    x2, edge_index2, edge_feat2,batchs = [],[],[],[]
    for i,agraph in enumerate(gfeat):
        x, edge_index, edge_feat = agraph
        #print(x.shape,edge.shape,edge_f.shape,edge_feat.shape[0])
        x2.append(x)
        edge_feat2.append(edge_feat)
        #print(edge_feat.shape[0]*i)
        edge_index = edge_index + x.shape[0]*i
        edge_index2.append(edge_index)
        batchs.append(np.full((x.shape[0]),i))
        #cat
    x2 = np.concatenate(x2)
    edge_index2 = np.concatenate(edge_index2)
    edge_feat2 = np.concatenate(edge_feat2)
    batch = np.concatenate(batchs)

    return x2, edge_index2, edge_feat2, batch

class EdgeModel(torch.nn.Module):
    def __init__(self, F_x, F_e, F_u, hiddens, dr=0.):
        super(EdgeModel, self).__init__()
        in_channels = 2*F_x + F_e + F_u
        self.edge_mlp = Seq(Lin(in_channels, hiddens),
                        nn.Dropout(p=dr),
                        ReLU(), 
                        Lin(hiddens, hiddens+5),
                        nn.Dropout(p=dr))

    def forward(self, src, dest, edge_attr,u, batch):
        # source, target: [E, F_x], where E is the number of edges.
        # edge_attr: [E, F_e]
        # u: [B, F_u], where B is the number of graphs.
        # batch: [E] with max entry B - 1.
        out = torch.cat([src, dest, edge_attr,u[batch]], 1)
        return self.edge_mlp(out)

class NodeModel(torch.nn.Module):
    def __init__(self, F_x, F_e,F_u,h1,h2,h3, dr=0.):
        super(NodeModel, self).__init__()
        in_channels = F_x + F_e
        #print(in_channels,h1,h2,h3)
        self.node_mlp_1 = Seq(Lin(in_channels, h1), 
                          nn.Dropout(p=dr),
                          ReLU(), 
                          Lin(h1, h2),
                          nn.Dropout(p=dr))
        self.node_mlp_2 = Seq(Lin(h2+F_x+F_u, h2), 
                          nn.Dropout(p=dr),
                          ReLU(), 
                          Lin(h2, h3),
                          nn.Dropout(p=dr))

    def forward(self, x, edge_index, edge_attr,u, batch):
        # x: [N, F_x], where N is the number of nodes.
        # edge_index: [2, E] with max entry N - 1.
        # edge_attr: [E, F_e]
        # u: [B, F_u]
        # batch: [N] with max entry B - 1.
        row, col = edge_index
        out = torch.cat([x[row], edge_attr], dim=1)
        out = self.node_mlp_1(out)
        out = scatter_mean(out, col, dim=0, dim_size=x.size(0))
        out = torch.cat([x, out, u[batch]], dim=1)
        return self.node_mlp_2(out)

class GlobalModel(torch.nn.Module):
    def __init__(self, F_x, F_e, F_u, h1,h2, dr=0.):
        super(GlobalModel, self).__init__()
        in_channels = F_x + F_e + F_u
        #in_channels = F_x + F_u 
        self.global_mlp = Seq(Lin(in_channels, h1),
                          nn.Dropout(p=dr),
                          ReLU(),
                          Lin(h1, h2),
                          nn.Dropout(p=dr))

    def forward(self, x, edge_index, edge_attr, u, batch):
        row, col = edge_index
        edge_u = scatter_mean(edge_attr, col, dim=0, dim_size=x.size(0))
        out = torch.cat([u, scatter_mean(x, batch, dim=0), scatter_mean(edge_u, batch, dim=0)], dim=1)
        return self.global_mlp(out)


class GNet(torch.nn.Module):
    def __init__(self, nb_classes,dropout):
        super(GNet, self).__init__()
        # deepgoplus
        self.nb_classes = nb_classes
        ########## for network 1
        # F_x, F_e, F_u, hiddens
        # edge hidden+5
        # node 3*hidden
        # glob hidden+5
        self.metalayer1 = MetaLayer(EdgeModel(12,1,20,5,dropout), NodeModel(12,10,20,15,20,20,dropout),GlobalModel(20,10,20,25,30,dropout))
        self.metalayer2 = MetaLayer(EdgeModel(20,10,30,5,dropout), NodeModel(20,10,30,20,25,20,dropout), GlobalModel(20,10,30,25,20,dropout))
        self.metalayer5 = MetaLayer(EdgeModel(23,10,20,5,dropout), NodeModel(23,10,20,23,10,1,dropout), GlobalModel(10,1,20,10,1,dropout))

    def forward(self, deepgo, x, u,edge, edge_attr, batch):
        # deepgocnn
        deepgo = torch.squeeze(deepgo,1)
        node_dpg = deepgo.view(self.nb_classes*(torch.max(batch)+1),1)

        #print('x',x.size(),'edge',edge.size(),'edgeatt',edge_attr.size(),'b',batch.size(),'dpg size',dpg.size(),'node_dpg',node_dpg.size(),'u',u.size())
        u = torch.squeeze(u[:,:,1:21],1)

        x1, edge_attr1,u = self.metalayer1(x, edge, edge_attr=edge_attr, u=u, batch=batch)
        x1, edge_attr1,u = self.metalayer2(x1, edge, edge_attr=edge_attr1,u=u,  batch=batch)
        #x1, edge_attr1,u = self.metalayer3(x1, edge, edge_attr=edge_attr1,u=u,  batch=batch)
        #x1, edge_attr1,u = self.metalayer4(x1, edge, edge_attr=edge_attr1,u=u,  batch=batch)
        # append node_dpg to x1
        x1 = torch.cat([x1,node_dpg,x[:,10:]],1)
        x1, edge_attr1,_ = self.metalayer5(x1, edge, edge_attr=edge_attr1, u=u, batch=batch)
        #print('x1',x1.size(),'edge',edge.size(),'edgeatt',edge_attr.size(),'b',batch.size(),'dpg size',dpg.size(),'node_dpg',node_dpg.size())
        # last layer
        #
        out = x1.view(torch.max(batch)+1,1,-1)
        return out

class Deepgo(nn.Module):
    def __init__(self, nb_classes):
        super(Deepgo, self).__init__()
        self.seqlen = 2000
        ########## for network 1
        self.block1 = nn.ModuleList()
        for kernel_size in range(8,129,8):
            cv1d = nn.Conv1d(4, 512, kernel_size=kernel_size, stride=1)#, bias=False)
            nn.init.xavier_uniform_(cv1d.weight)
            self.block1.append(cv1d)

            maxpool = nn.MaxPool1d(self.seqlen - kernel_size + 1)
            self.block1.append(maxpool)

        ########## network2
        #self.fc1 = nn.Linear(8192+1280, nb_classes)
        self.fc1 = nn.Linear(1280, nb_classes)

    def forward(self, x, u):
        """
        # 1d residual
        xs = []
        for i in range(int(len(self.block1)/2)):
            layer = self.block1[i*2]
            maxp1d = self.block1[i*2+1]
            x_tmp = maxp1d(layer(x))
            x_tmp = torch.flatten(x_tmp, start_dim=1)
            x_tmp = torch.unsqueeze(x_tmp,1)
            xs.append(x_tmp)
        x2 = torch.cat(xs, 2)
        x2 = torch.cat([x2,u[:,:,21:]],2)
        x2 = self.fc1(x2)
        """
        x2  = self.fc1(u[:,:,21:])
        return x2

if __name__ == '__main__':
    main()
