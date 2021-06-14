# TBRU_RIFcost
==============================


[![DOI](https://zenodo.org/badge/XXXXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXXX)

The data pertinent to the analyses can be on zenodo, with DOI-references:
[![DOI](https://www.zenodo.org/badge/DOI/10.5281/zenodo.4903635.svg)](https://doi.org/10.5281/zenodo.4903635)

### This study was performed by:
Andrej Trauner, Amir Banaei-Esfahani, Sebastian M. Gygli, Philipp Warmer, Julia Feldmann, Mattia Zampieri, Sonia Borrell, Ben C. Collins, Christian Beisel, Ruedi Aebersold and Sebastien Gagneux.

## Abstract for Expression dysregulation as a mediator of fitness costs in antibiotic resistance
Antimicrobial resistance (AMR) poses a threat to global health and the economy. Rifampicin-resistant Mycobacterium tuberculosis accounts for a third of the global AMR burden. Gaining the upper hand on AMR requires a deeper understanding of the physiology of resistance. AMR often results in a fitness cost in absence of drug. Identifying the molecular mechanisms underpinning this cost could help strengthen future treatment regimens. Here, we used a collection of M. tuberculosis strains providing an evolutionary and phylogenetic snapshot of rifampicin resistance, and subjected them to genome-wide transcriptomic and proteomic profiling to identify key perturbations of normal physiology. We found that the clinically most common rifampicin resistance-conferring mutation RpoB Ser450Leu imparts considerable gene expression changes, many of which are mitigated by the compensatory mutation in RpoC Leu516Pro. However, our data also provide evidence for pervasive epistasis: the same resistance mutation imposed a different fitness cost and functionally distinct changes to gene expression in genetically unrelated clinical strains. Finally, we report a likely post-transcriptional modulation of gene expression that is shared in most of the tested strains carrying RpoB Ser450Leu, resulting in an increased abundance of proteins involved in central carbon metabolism. These changes contribute to a more general trend, in which the disruption of the composition of the proteome correlates with the fitness cost of the RpoB Ser450Leu mutation in different strains.

### Importance
Antimicrobial resistance poses a threat to global health and the economy. It is widely accepted that, in the absence of antibiotics, drug resistance mutations carry a fitness cost. In the case of rifampicin resistance in fast-growing bacteria, this cost stems from a reduced transcription rate of the RNA polymerase resulting in slower ribosome biosynthesis. However, this relationship does not apply in the slow-growing Mycobacterium tuberculosis, where the true mechanism of fitness cost of rifampicin resistance as well as the impact of compensatory evolution remain unknown. Here we show, using global transcriptomic and proteomic profiling of selected M. tuberculosis mutants and clinical strains, that rifampicin resistance due to RpoB Ser450Leu in M. tuberculosis imparts complex gene expression changes that are a target for compensation, and might therefore lie at the root of the fitness cost of resistance. Understanding these dependencies might result in better combination treatments in the future.

Here we provide a record of the analyses and offer a structure for the 

## Project Organization
------------
```
    ├── LICENSE
    ├── README.md          <- The top-level README to serve as a guide for future analyses.
    ├── data
    │   ├── README.md      <- NB, these data are not here. To be downloaded from Zenodo [10.5281/zenodo.4903635]
    │   ├── RNAseq         <- Data, external, raw and processed pertinent to transcriptomic profiling of our samples
    │   ├── Proteome       <- Data, external, raw and processed pertinent to proteomic profiling of our samples
    │   ├── Growth_Cruves  <- Data for the growth curves of bacterial cultures reported in this study.
    │   └── Cost_data      <- Data for the determination of the overall cost of expression.
    │
    ├── notebooks          <- Jupyter notebooks with the analysis.
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         use for creating a virtual environment for the analysis `pip freeze > requirements.txt`
    │
    └── src                <- Source code for use in this project.
        │
        ├── data           <- Scripts used to process raw sequencing data
        │
        ├── DESeq          <- R script for differential expression analysis
        │
        └── TbX_module_ver1.py         <- Custom script for data manipulation - needed for the Jupyter notebook.
```

## SETTING UP THE VIRTUAL PYTHON ANALYSIS ENVIRONMENT [cost-venv]
----------

```
#With virtualenv on Mac/Linux
python -m venv cost-venv #create
source cost-venv/bin/activate #activate
python -m pip install -r requirements.txt #install dependencies
```
```
#With virtualenv on Windows
python -m venv cost-venv #create
cost-venv\Scripts\activate.bat #activate
python -m pip install -r requirements.txt #install dependencies
```
```
#With conda
conda create --name cost-venv --file requirements.txt
```