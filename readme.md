# Welcome to the homepage of Temporal Enrichmentnt Analysis tooL(TEAL)

This repository contains the codes for the python package of the paper: 'Understanding the Impact of Gene Ontology Evolution with Temporal Semantic Similarity‘

## About 

Although Gene Ontology Knowledge Base (GOKB) is widely used in biological studies, the evolution of GOKB and its impact on studies has not attracted enough attention, espcially in downstream analysis, in which GOKB is one of the commonly used background knowledge base. 

By taking the temporal component of GOKB into consideration, TEAL enables users to compute the semantic similarity across versions. For studies that shares research targets, TEAL can help researchers build the consensus based on enrichment analysis outcomes. Also, TEAL can also be used to evaluate the reusability of outcomes of previously published studies, in which out-of-date versions of GOKB were used. 



TEAL can help you to evaluate if the results from an enrichment analysis results are still comparable or not. 

For enrichment analysis. 
step 1) provide a timepoint ('yyyy-mm-dd')
     2) provide the input gene list
     3) a dataframe will be returned

For visualization
step 1) provide the location of two term list.
     2) provide the time-point coressponding to the set 
