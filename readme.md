# Welcome to the homepage of Temporal Enrichmentnt Analysis tooL(TEAL)

This repository contains the codes for the python package of the paper: 'Understanding the Impact of Gene Ontology Evolution with Temporal Semantic Similarity‘

## About 

Although Gene Ontology Knowledge Base (GOKB) is widely used in biological studies, the evolution of GOKB and its impact on studies has not attracted enough attention, espcially in downstream analysis, in which GOKB is one of the commonly used background knowledge base. 

By taking the temporal component of GOKB into consideration, TEAL enables users to compute the semantic similarity across versions. For studies that shares research targets, TEAL can help researchers build the consensus based on enrichment analysis outcomes. Also, TEAL can also be used to evaluate the reusability of outcomes of previously published studies, in which out-of-date versions of GOKB were used. 

## how to use TEAL

### enrichment analysis 
TEAL can be used to conduct enrichment analysis using the version of GOKB at a given timepoint. To run a enrichment analysis, users need to provide:
   1) a gene list
   2) a timepoint ('yyyy-mm-dd')
   3) saving path

The output file includes ['GO terms','p value','k','M','N','n','adjust p value']. 


### enrichment analysis results visualization
Besides conducting enrichment analysis, TEAL is also able to visualize the enriched terms, from two sets. This makes it easier for users to know during the evolution of GOKB, how the enrichment analysis results was influenced by GOKB evolution, even the input gene set stays the same. Also, users can use the visualization module to see how the enriched terms of two studies are connected, regardless of GO versions. 

To visualize the results, users need to provide:
1) GO term set A
2) GO term set B
3) the corresponding version of GO term set A
4) (optional) the corresponding version of GO term set B
5) saving path
![Nsp8](https://github.com/chestnzu/TEAL/assets/40864288/52eabf25-5d85-4f0d-bcf7-539d8b6c8d56)
Figure 1. Example visualization outcomes. 


For visualization
step 1) provide the location of two term list.
     2) provide the time-point coressponding to the set 
