# TBRU_RIFcost
==============================

## Introduction

This document illustrates the structure of the data deposited on Zenodo. You can download the data from here:
[![DOI](https://www.zenodo.org/badge/DOI/10.5281/zenodo.4903635.svg)](https://doi.org/10.5281/zenodo.4903635)


```
.\
├───Cost_data
├───Growth_Curves
├───Proteome
└───RNAseq
```

## RNAseq
Processed data used for the paper - raw data available at the ArrayExpress repository of the European Bioinformatics Institute under the E-MTAB-7359 accession.

### Processed transcriptomics data
Files containing read counts per gene obtained by htseq-count. The sample list as well as a collation of all the data.

- *.htc
- EVOLUTION_SET_compensation_genesets.csv
- EVOLUTION_SET_RNAseq_dataframe.tsv [**SUMMARY TEXT FILE**]
- EVOLUTION_SET_rnaseq_samples.txt
- GENETIC_DIVERSITY_RNAseq_dataframe.tsv [**SUMMARY TEXT FILE**]
- GENETIC_DIVERSITY_rnaseq_samples.txt

### Differential expression summary
Starting with the htc files as input, we used DESeq2 to perform differential expression analysis on the data. The analysis pertinent to each file is outlined next to it.

- DESeq_EVOLUTION_SET_compensation.csv [RifR vs rest of dataset]
- DESeq_EVOLUTION_SET_comp_rpoBrpoBE.csv [RifR vs RifRevo]
- DESeq_EVOLUTION_SET_rif_RS.csv [DS vs RifR]
- DESeq_EVOLUTION_SET_rif_wtE_rpoBE.csv [DSevo vs RifRevo]
- DESeq_EVOLUTION_SET_rif_wt_rpoBE.csv [DS vs RifRevo]
- DESeq_GENETIC_DIVERSITY_rif_RS.csv [all RifR vs all DS]
- DESeq_GENETIC_DIVERSITY_all.csv [all RifR vs all DS adjusted for genetic background]
- DESeq_GENETIC_DIVERSITY_pair_N0052.csv [N0052 RifR vs DS]
- DESeq_GENETIC_DIVERSITY_pair_N0072.csv [N0072 RifR vs DS]
- DESeq_GENETIC_DIVERSITY_pair_N0145.csv [N0145 RifR vs DS]
- DESeq_GENETIC_DIVERSITY_pair_N0155.csv [N0155 RifR vs DS]
- DESeq_GENETIC_DIVERSITY_pair_N0157.csv [N0157 RifR vs DS]


### Support files for analysis
These files contain the mappers for the genesets used for the analysis as well as annotation information for *Mycobacterium tuberculosis* H37Rv.

- Geneset_mapper.txt
- H37Rv.json
- H37Rv_extra_knowledge.txt
- mycobacterium_tuberculosis_h37rv_2_genome_summary_per_gene.txt
- Peterson_module_genes.tsv
- Peterson_module_TFs.csv

### Support files for the network analysis
Mostly contains the positional information for the randomly generated graphs for the sake of reproducibilty.

GRAPH_ANALYSIS_node_positions.json
GRAPH_ANALYSIS_RpoB_modules_nodes.json
GRAPH_ANALYSIS_RpoB_modules_positions.json
GRAPH_ANALYSIS_RpoB_nodes_expanded.json
GRAPH_ANALYSIS_RpoB_node_positions_expanded.json


## Proteome

Processed data used for the analyses are found in the files listed below. The mass spectrometry proteomics data have been deposited  to  the  ProteomeXchange  Consortium  via  the  PRIDE  partner  repository  with  the dataset  identifier  PXD011568.

### Proteomics data
Label-free quantification for each sample for the two strains sets described in the paper are collated in the following two files.

- EVOLUTION_SET_SWATH_MS_dataframe.txt [**SUMMARY TEXT FILE**]
- GENETIC_DIVERSITY_SWATH_MS_dataframe.txt [**SUMMARY TEXT FILE**]


### Proteomic differential expression
We used MSstat package to derive differential expression analyses as outlined in the short description bracket.

- MSstats_EVOLUTION_SET_compensation.tsv [RifR vs rest of dataset]
- MSstats_EVOLUTION_SET_rpoBCrpoB.tsv [RifRevo vs RifR]
- MSstats_EVOLUTION_SET_rpoBCwt.tsv [RifRevo vs DS]
- MSstats_EVOLUTION_SET_rpoBCwtE.tsv [RifRevo vs DSevo]
- MSstats_EVOLUTION_SET_rpoBwt.tsv [RifR vs DS]
- MSstats_GENETIC_DIVERSITY_L1.txt [all RifR vs all DS for Lineage 1 - N0072, N0157]
- MSstats_GENETIC_DIVERSITY_L2.txt [all RifR vs all DS for Lineage 2 - N0052, N0145, N0155]
- MSstats_GENETIC_DIVERSITY_RpoB.txt [all RifR vs all DS]
- MSstats_GENETIC_DIVERSITY_RpoB_All.txt [all RifR vs all DS adjusted for genetic background]
- MSstats_GENETIC_DIVERSITY_RpoB_N0052.txt [N0052 RifR vs DS]
- MSstats_GENETIC_DIVERSITY_RpoB_N0072.txt [N0072 RifR vs DS]
- MSstats_GENETIC_DIVERSITY_RpoB_N0145.txt [N0145 RifR vs DS]
- MSstats_GENETIC_DIVERSITY_RpoB_N0155.txt [N0155 RifR vs DS]
- MSstats_GENETIC_DIVERSITY_RpoB_N0157.txt [N0157 RifR vs DS]

## Growth_Curves
Data pertinent to the measurement of growth rates of all the strains analysed in this study.

- EVOLUTION_SET_GROWTH_CURVES.csv
- EVOLUTION_SET_GROWTH_CURVES_HEMIN.csv
- GENETIC_DIVERSITY_GROWTH_CURVES.csv

## Cost_data
Datafiles pertinent to the analysis of the fitness cost of the RpoB Ser450Leu substitution as analysed through the lens of global changes in the expression profile.

### External data
External data containing different metrics for the estimation of the metabolic costs of protein sysnthesis based on the constituent amino acids. And annotation and sequence data for each protein coding sequence in the H37Rv genome.
- Barton_2010_supp_1.csv
- Barton_2010_supp_2.csv
- GCF_000195955.2_ASM19595v2_feature_table.txt
- GCF_000195955.2_ASM19595v2_protein.faa

### Summary proteomic/transcriptomic data
The data below were used to estimate the impact of RpoB Ser450Leu on overall gene expression in different strains of *Mycobacterium tuberculosis*. These files are almost identical to the SUMMARY proteomic and transcriptomic data above.

- GENETIC_DIVERSITY_SUMMARY_PROTEOME.csv
- GENETIC_DIVERSITY_SUMMARY_RNASEQ.csv
