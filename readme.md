# Welcome to the homepage of Temporal Enrichmentnt Analysis tooL(TEAL)

This repository contains the codes for the python package of the paper: 'Temporal Semantic Similarity: Exploring the Evolution of Knowledge in Biomedical Ontologies‘

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
            Figure 1. Example visualization of Nsp8 and its interacted human protein. The list of interacted human proteins is retrieved from the paper 'https://www.nature.com/articles/s41586-020-2286-9'. We conducted enrichment analysis by using TEAL with two different versions of GOKB: 2017 Oct and 2021 Dec and visualized using TEAL. Red nodes are uniquely enriched when applying the former release of the GOKB, while blue nodes represent uniquely enriched terms in Dec 2021. Yellow nodes are terms enriched in both 2017 and 2021.


### Temporal Semantic Similarity 

As GOKB updates almost monthly, both the Gene Ontology and Gene Ontology Annotation keep changing. To capture the temporal component of GOKB and reduce its influence on enrichment analysis, we proposed a method named 'Temporal Semantic Similarity' (TSS). For more details, you can read our paper. 

To compute the TSS score, users must provide:
1) a time point. e.g., '2017-10-25'
2) aspects. default:'biological_process'
3) relationships. If only 'is_a' relationship is considered. 
4) GO term set A and GO term set B
