#!/usr/bin/env python
# coding: utf-8

# In[2]:


import obonet,re,gzip,time,os
import pandas as pd
import networkx as nx
import numpy as np
from scipy.spatial import distance
import requests
import networkx
import matplotlib.pyplot as plt
from collections import Counter


# In[3]:


def obonet_modif(link):
    file=obonet.read_obo(link,ignore_obsolete=False)
    added_dict={}
    for term in file.nodes:
        attributes=file.nodes[term]##a dict
        is_obsolete=attributes.get('is_obsolete','false')=='true'
        if is_obsolete:##如果term已经obsolete，则直接给obsolete标签，否则给valid标签。因为在nodes中的点不会存在除这两种状态外的其他状态
            attributes['status']='obsolete'
        else:
            attributes['status']='valid'
        alter=attributes.get('alt_id',False)
        if alter:
            added_dict.update({x:term for x in alter})
    file.add_nodes_from([(node,{'current':id,'status':'merged'}) for (node,id) in added_dict.items()])##给altered merged标签
    
    return file
################find the url to corresponding version of GO and GOA#########################
############################################################################################
r=requests.get('https://go-data-product-release.s3.amazonaws.com/?list-type=2&delimiter=/&prefix=')
pattern=re.compile('>[0-9]{4}\-[0-9]{2}\-[0-9]{2}/')
all_date=pattern.findall(r.text)
for x in all_date:
    all_date.remove(x)
    x=x.replace('>','')
    x=x.replace('/','')
    all_date.insert(0,x)
all_date.sort()
def corresponding_time(time):
    annotation_online_path='http://release.geneontology.org/0000-00-00/annotations/goa_human.gaf.gz'
    ontology_online_path='http://release.geneontology.org/0000-00-00/ontology/go.obo'
    for x in all_date:
        if time<x:
            index=all_date.index(x)
            break
        else:
            index=-1
    if index==0:
        time1=all_date[0]
    else:
        time1=all_date[index-1]
    annotation_online_path = annotation_online_path.replace('0000-00-00', time1)
    ontology_online_path = ontology_online_path.replace('0000-00-00', time1)
    return [annotation_online_path, ontology_online_path]


# In[9]:


class network:
    def __init__(self,annotation='http://current.geneontology.org/annotations/goa_human.gaf.gz',ontology='http://purl.obolibrary.org/obo/go.obo',directed_only=True):
        if re.match(r'^http:/{2}\w.+$',annotation):
#             file=requests.get(annotation)
#             with open('annotation.gz','wb') as f:
#                 f.write(file.content)
#             f1 = gzip.open('annotation.gz', 'rb')
#             f_out = open('text.txt', 'w')
#             file_content = f1.read()
#             f_out.write(file_content.decode('utf-8'))
#             f.close()
#             f_out.close()
#             with open('text.txt', 'r') as out:
#                 num = 0
#                 for line1 in out.readlines():
#                     if line1.startswith('!'):
#                         num += 1
#                         continue
#                     else:
#                         break
#             annotation_f = pd.read_csv('text.txt', sep='\t', skiprows=num, header=None,usecols=[1,3,4],names=['gene product','relationship','GO term'])
            output=requests.get(annotation)
            with open('annotation.gz','wb') as f:
                f.write(output.content)
            f1 = gzip.open('annotation.gz', 'rb')

            items=[]
            name=['gene product','relationship','GO term','evidence code']
            for line in f1.readlines():
                line=line.decode('utf-8')
                if line.startswith('!')==False:
                    elements=line.split('\t')
                    items.append((elements[1],elements[3],elements[4],elements[6]))
            os.remove('annotation.gz')
            annotation_f=pd.DataFrame.from_records(items,columns=name)
        else:
            with open(annotation, 'r') as out:
                num = 0
                for line1 in out.readlines():
                    if line1.startswith('!'):
                        num += 1
                        continue
                    else:
                        break
            annotation_f = pd.read_csv(annotation, sep='\t', skiprows=num, header=None,usecols=[1,3,4,6],names=['gene product','relationship','GO term','evidence code'])
        t1=time.time()
        G_a = nx.from_pandas_edgelist(annotation_f, source='gene product', target='GO term', edge_key='relationship', create_using=nx.MultiDiGraph())
        G_o = obonet.read_obo(ontology)
        G_o_obsolete=obonet_modif(ontology)
        annotation_gene=annotation_f[['gene product','GO term']].copy()
        annotation_gene.drop_duplicates(inplace=True)
        G_a_gene=nx.from_pandas_edgelist(annotation_gene,source='gene product',target='GO term',create_using=nx.MultiDiGraph())

        if directed_only:
            G_o.remove_edges_from([(x,y) for x,y,z in G_o.edges(keys=True) if z!='is_a'])
        self.ontology=G_o
        self.annotation=G_a
        self.annotation_genes=G_a_gene
        self.ontology_with_obsolete=G_o_obsolete

    def categories(self):
        G_o = self.ontology
        bp = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'biological_process']
        mf = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'molecular_function']
        cc = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'cellular_component']
        return [bp, mf, cc]

    def information_content_value(self,ignore_duplicates=False):
        bp = self.categories()[0]
        mf = self.categories()[1]
        cc = self.categories()[2]
        ontology = self.ontology
        if ignore_duplicates==False:
            annotation = self.annotation
        else:
            annotation=self.annotation_genes
        IC_record_bp = {}
        IC_record_mf = {}
        IC_record_cc = {}

        aspects = ['bp', 'mf', 'cc']
        for aspect in aspects:
            IC_file = eval('IC_record_' + aspect)
            for term in eval(aspect):
                try:
                    child = [subterm for subterm in nx.ancestors(ontology, term)]
                except networkx.exception.NetworkXError:
                    IC_file[term] = len(annotation.in_edges(term))

                else:
                    child.append(term)
                    ICs = np.sum([len(annotation.in_edges(child_term)) for child_term in child])
                    IC_file[term]=ICs  # 这一步用于计算information content
        return [IC_record_bp, IC_record_mf, IC_record_cc]
    
class networks(network):
    def __init__(self,annotation,ontology,directed_only=True,ignore_duplicates=False,aspect='biological_process'):
        super().__init__(annotation,ontology,directed_only)
        self.IC_values=super().information_content_value(ignore_duplicates)
        self.aspect=aspect
####here, the MICA is based on IC values.         
    def MICA(self,a,b,return_type='ID',subterms=None):
        root_list=['GO:0008150','GO:0003674','GO:0005575']
        namespace=['biological_process','molecular_function','cellular_component']
        aspect=self.aspect
        IC_values=self.IC_values
                ##对于正常情况
        if subterms:
            ontology=self.ontology.subgraph(subterms)
        else:
            ontology=self.ontology
        idx_a=namespace.index(ontology.nodes[a]['namespace'])
        idx_b=namespace.index(ontology.nodes[b]['namespace'])
        ##解决特殊情况
        if a==b:
            if return_type=='IC':
                return IC_values[idx_a][a]
            else:
                return a
        if a in root_list:
            if return_type=='IC':
                return IC_values[idx_a][a]
            else:
                return a
        if b in root_list:
            if return_type=='IC':
                return IC_values[idx_b][b]
            else:
                return b
        

        
        parentA = [x for x in nx.descendants(ontology, a)]
        parentA.append(a)
        parentB = [x for x in nx.descendants(ontology, b)]
        parentB.append(b)
        common_ancestors = set(parentA) & set(parentB)
        
        mica_value=np.inf
        mica=''
        for x in common_ancestors:
            idx1=namespace.index(aspect)
            x_score=IC_values[idx1][x]
            if x_score<mica_value:
                mica=x
                mica_value=x_score
            else:
                continue
        if return_type=='IC':
            return mica_value
        else:
            return mica
    def similarity_score_position(self,image=False):
        aspects=['biological_process','molecular_function','cellular_component']
        aspect=self.aspect
        idx=aspects.index(aspect)
        IC=self.IC_values[idx]
        terms=list(IC.keys())
        denominator=max(IC.values())
        maxp=-np.log(1/(denominator+0.001))
        network1=self.ontology.subgraph(terms)
        IC_non_zero_id=[x for x,y in IC.items() if y!=0]
        have_d=[x for x in IC_non_zero_id if len(nx.ancestors(network1,x))>2]
        for term in have_d:
            if IC[term]>50:
                continue
            n=0
            descendants=nx.ancestors(network1,term)
            for de in descendants:
                if IC[de]==0:
                    continue
                else:
                    n+=1
            if n<2:
                have_d.remove(term)
        IC_non_zero=[y for x,y in IC.items() if y!=0 and x in have_d]
        x0=[(-np.log(x/(denominator+0.001)))/maxp for x in IC_non_zero]
        x1=Counter(x0)
        if image==True:
            fx=dict(sorted(x1.items()))
            fig,ax=plt.subplots()
            f1=ax.scatter(fx.keys(),fx.values())
            plt.show()
        return x0


# In[1]:


class similarity:
    def __init__(self,endtime,starttime=None,directed_only=True,ignore_duplicates=False,aspect='biological_process'):#time should be in yyyy-mm-dd format

        if starttime==None:
            info2=info1=corresponding_time(endtime)
        else:
            info1=corresponding_time(starttime)
            info2=corresponding_time(endtime)
        #self.start_IC=networks(info1[0],info1[1],directed_only=True)
        self.end_IC=networks(info2[0],info2[1],directed_only,ignore_duplicates,aspect)
        self.score_evaluate=self.end_IC.similarity_score_position()
        self.aspect=aspect
    def set_similarity(self, seta, setb): #,only_interested=None,):
        namespace=['biological_process','molecular_function','cellular_component']
        if type(seta).__name__ != 'list':
            if seta.endswith('.txt'):
                sa = pd.read_csv(seta, sep='\t', header=None, names=['GO'])
                seta = sa['GO'].to_list()
            else:
                return 'wrong input'
        if type(setb).__name__ != 'list':
            if setb.endswith('.txt'):
                sb = pd.read_csv(setb, sep='\t', header=None, names=['GO'])
                setb = sb['GO'].to_list()
            else:
                return 'wrong input'
        aspect=self.aspect
        e_ont=self.end_IC.ontology
        e_net=self.end_IC
        IC_values=e_net.IC_values#information_content_value()   
        e_obs_nodes=e_net.ontology_with_obsolete.nodes
        if seta==[] or setb==[]:
            return 'empty input'
        setA=[]
        for term in list(seta): 
            attribute=e_obs_nodes[term]
            status=attribute['status']
            if status=='valid':
                setA.append(term)
            elif status=='merged': ##将merged term换成其alternative terms
                sub=attribute['current']
                setA.append(sub)
            else:                  ##对于obsolete term,查看其是否有replace_by 或者 consider，如果都没有，那么直接删除
                replaceby=attribute.get('replaced_by',[])
                consider=attribute.get('consider',[])
                a1=[x for x in replaceby]
                a2=[x for x in consider]
                setA.extend(a1)
                setA.extend(a2)

        setB=[]
        for term in list(setb): 
            attribute=e_obs_nodes[term]
            status=attribute['status']
            if status=='valid':
                setB.append(term)
            elif status=='merged': ##将merged term换成其alternative terms
                sub=attribute['current']
                setB.append(sub)
            else:                  ##对于obsolete term,查看其是否有replace_by 或者 consider，如果都没有，那么直接删除
                replaceby=attribute.get('replaced_by',[])
                consider=attribute.get('consider',[])
                a1=[x for x in replaceby]
                a2=[x for x in consider]
                setB.extend(a1)
                setB.extend(a2)
        setA=[x for x in setA if e_ont.nodes[x]['namespace']==aspect]
        setb=[x for x in setB if e_ont.nodes[x]['namespace']==aspect]
        if setb==[] or setA==[]:
            return 'Unrelated'
        onlyA=set(setA)-set(setb)
        onlyB=set(setb)-set(setA)
        common_set=set(setA)&set(setb)
        IC_value=IC_values[namespace.index(aspect)]
        max_value=max(IC_value.values())
        max_score=-np.log(1/max_value)
        #scores=sum([-np.log(IC_value[x]/max_value)/max_score for x in common_set])*2
        scores=len(common_set)*2#*(len(common_set)/len(set(setA)|set(setb)))
        for term1 in onlyA:
            value=min([e_net.MICA(term1,x,'IC') for x in setb])
            value=-np.log(value/max_value)
            value_norm=value/max_score
            scores+=value_norm
        for term2 in onlyB:
            value=min([e_net.MICA(term2,x,'IC') for x in setA])
            value=-np.log(value/max_value)
            value_norm=value/max_score
            scores+=value_norm
        final_score=scores/(len(setA)+len(setb))
        print(final_score)
        score_quantile=self.score_evaluate
        #score_quantile=list(set(score_quantile))
        intervals=np.percentile(score_quantile,[30,60,80])
        if final_score==1:
            return 'Identical'
        if final_score<=intervals[0]:
            return 'Unrelated'
        if final_score>intervals[0] and final_score<=intervals[1]:
            return 'Related'
        if final_score>intervals[1] and final_score<=intervals[2]:
            return 'Similar'
        else:
            return 'Compatible'


# In[ ]:




