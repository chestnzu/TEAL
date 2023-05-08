#!/usr/bin/env python
# coding: utf-8

# In[1]:


import graphviz
import networkx as nx


# In[2]:


from computing_similarity_modified import *
from enrichment_analysis import *


# In[74]:


class viz_to_root:
    def __init__(self, seta, setb, onta, ontb=None):
        self.set1 = seta
        self.set2 = setb

        if ontb == None:
            ontb = onta
        self.ont1 = onta
        self.ont2 = ontb

    def draw(self,savepath=None):
        ont1 = self.ont1
        ont2 = self.ont2
        seta=self.set1
        setb=self.set2
        seta_edges = []
        for x in seta:
            namespace = ont1.nodes[x]['namespace']
            if namespace == 'biological_process':
                edges = nx.all_simple_edge_paths(ont1, x, 'GO:0008150')
            elif namespace == 'molecular_function':
                edges = nx.all_simple_edge_paths(ont1, x, 'GO:0003674')
            else:
                edges = nx.all_simple_edge_paths(ont1, x, 'GO:0005575')
            for x in list(edges):                
                    [seta_edges.append(elements) for elements in x if elements not in seta_edges]


        setb_edges = []
        for x in setb:
            namespace = ont2.nodes[x]['namespace']
            if namespace == 'biological_process':
                edges = nx.all_simple_edge_paths(ont2, x, 'GO:0008150')
            elif namespace == 'molecular_function':
                edges = nx.all_simple_edge_paths(ont2, x, 'GO:0003674')
            else:
                edges = nx.all_simple_edge_paths(ont2, x, 'GO:0005575')
            for x in list(edges):                
                    [setb_edges.append(elements) for elements in x if elements not in setb_edges]

        g_edge=seta_edges+setb_edges
        g_edge=list(set(g_edge))
        g=graphviz.Digraph(comment='ontology_visualization')
        g.graph_attr['rankdir'] = 'BT'
        g.format='png'
        for x in g_edge:
            source=x[0].replace(':','_')
            target=x[1].replace(':','_')
            g.edge(source,target)

        common=set(seta)&set(setb)
        onlya = set(seta) - set(setb)
        onlyb = set(setb) - set(seta)

        for x in common:
            x=x.replace(':','_')
            g.node(x,color='yellow',style='filled')
        for x in onlya:
            x=x.replace(':','_')
            g.node(x,color='red',style='filled')
        for x in onlyb:
            x=x.replace(':','_')
            g.node(x,color='lightblue',style='filled')
        
        g.attr(layout='fdp')
        g.graph_attr['rankdir']='TB'
        if savepath==None:           
            g.render('./testfile{0}.gv',format='jpg',view=True)
        else:
            g.render(savepath+'goviz.gv',format='jpg',view=True)



