#!/usr/bin/env python
# coding: utf-8

# In[33]:


from computing_similarity_modified import *
import pandas as pd
import obonet
import networkx as nx
#####generate ontology file for bingo,below example generate file of biological process terms
from scipy.stats import hypergeom


# In[44]:


class bingo(networks):
    def __init__(self,time):
        file_url=corresponding_time(time)
        super().__init__(file_url[0],file_url[1])

    def export(self):
        ontology=self.ontology
        with open('./ontology.txt','w') as f:
            f.write('(curator=GO)(type=all)' + '\n')

            for x in ontology.nodes:
                name = ontology.nodes[x]['name']
                x1 = x.split(':')[1]
                if x in ['GO:0008150', 'GO:0005575', 'GO:0003674']:
                    f.write(x1 + '=' + name + '\n')
                else:
                    name = name + ' [isa: '
                    parents = ontology.nodes[x]['is_a']
                    for parent in parents:
                        parent_name = parent.split(':')[1]
                        name = name + str(parent_name) + ' '
                    name += ']'
                f.write(x1 + ' = ' + name + '\n')

        annotation=self.annotation
        df=nx.to_pandas_edgelist(annotation)
        df['target']=df['target'].apply(lambda x:x.split(':')[1]).astype('str')
        df['output'] = df['source'] + ' = ' + df['target']
        output = df['output'].to_list()
        with open('./annotation.txt', 'w') as f1:
            f1.write('(species=Homo Sapien)(type=all)(curator=GO)' + '\n')
            for line in output:
                f1.write(line + '\n')


# In[37]:


def fdr(p_vals):
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)

    fdr1 = [x * len(p_vals) / y for x, y in zip(p_vals, ranked_p_values)]
    fdr1 = [x if x <= 1 else 0 for x in fdr1]
    return fdr1

class enerich(networks):
    def __init__(self,time):
        file_url=corresponding_time(time)
        super().__init__(file_url[0],file_url[1])
        
        
        
    def enrichmentanalysis(self, inputgene):
        ontology = self.ontology
        annotation = self.annotation_genes
        annotation = nx.to_pandas_edgelist(annotation)
        associated = annotation.loc[annotation['source'].isin(inputgene)]
        related_go = list(set(associated['target'].to_list()))
        associated_go_list = []
        M = len(set(annotation.source.to_list()))  # M
        N = len(inputgene)  # N
        df = pd.DataFrame(columns=['GO term', 'p_value','k','M','N','n'])
        num = 0
        for x in related_go:
            ancestors = list(nx.descendants(ontology, x))
            associated_go_list.extend(ancestors)
            associated_go_list.append(x)
            associated_go_list = list(set(associated_go_list))
        for term in associated_go_list:

            descendants = list(nx.ancestors(ontology, term))
            descendants.append(term)
            k = len(associated.loc[associated['target'].isin(descendants), 'source'].drop_duplicates())  # k
            #k = len(set([x[0] for x in associated if x[1] in descendants]))
            n = len(annotation.loc[annotation['target'].isin(descendants), 'source'].drop_duplicates())  # n
            #n = len(set([x[0] for x in annotation.edges() if x[1] in descendants]))
            P = hypergeom.sf(k - 1, M, n, N)
            df.loc[num] = [term, P,k,M,N,n]
            num += 1
        df.sort_values(by='p_value',inplace=True)
        df['adjusted']=fdr(df['p_value'])
        return df

