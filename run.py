#!/usr/bin/env python
# coding: utf-8

# In[4]:


from enrichment_analysis import *
from enrichment_visualization import *
from computing_similarity_modified import *


# In[7]:


import sys
import argparse


# In[10]:


parser=argparse.ArgumentParser(description='TEAL')
parser.add_argument('type1',type=str,help='ea,vis,tss',default='ea') ## vis need 4 parameters,seta,setb,timea,timeb
                                                        ## enrichment analysis need 2, time and inputgenelist
                                                        ## tss need 6,seta, setb, endtime, directed_only, ignore_duplicates and aspect
parser.add_argument('--save_path',type=str,help='save_path')
parser.add_argument('--timea',type=str,help='time or endtime')
parser.add_argument('--lista',type=str,help='seta_path or gene list')
parser.add_argument('--listb',type=str,help='setb_path')
parser.add_argument('--timeb',type=str,help='timeb')
parser.add_argument('--directed',type=str,default=True,help='directed only')
parser.add_argument('--duplicates',type=str,default=False,help='ignore duplicates or not')
parser.add_argument('--aspect',type=str,default='biological_process',help='aspect of gene ontology')
args=parser.parse_args()


# In[9]:


if __name__=='__main__':
    conduction=args.type1
    print(conduction)
    lista=[]
    with open(args.lista) as f:
        for x in f.readlines():
            lista.append(x.rstrip())     
    if conduction=='ea':##conducting enrichment analysis
        f1=enerich(args.timea)
        f2=f1.enrichmentanalysis(lista)
        f2.to_csv(args.save_path+'/ea.csv')
    elif conduction=='vis':##visualization
        listb=[]
        with open(args.listb) as f:
            for x in f.readlines():
                listb.append(x.rstrip())  
        f1=viz_to_root(lista,listb,args.timea,args.timeb)
        f1.draw(args.save_path)
    elif conduction=='tss':
        listb=[]
        with open(args.listb) as f:
            for x in f.readlines():
                listb.append(x.rstrip()) 
        f1=similarity(args.timea,directed_only=args.directed,ignore_duplicates=args.duplicates,aspect=args.aspect)
        f2=f1.set_similarity(lista,listb)
        print(f2)

