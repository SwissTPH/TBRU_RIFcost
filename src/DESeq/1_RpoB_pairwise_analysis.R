#Analysis of TbX011 and TbX008 pairwise comparisons. 
#The main question is: can we identify the size fo the S450L perturbation.

#Set working directory
#raw/RNAseq/
#setwd("C:/Users/trauner/Documents/git/cost-of-resistance/data/raw/RNAseq/")
#setwd("~/Documents/Work/analysis_logs/cost-of-resistance/data/raw/RNAseq/")

#setwd("~/Documents/Work/BC2_home/cost-of-resistance/data/raw/RNAseq/")
#setwd("/Users/Zuco/Documents/Work/BC2_home/TBRU_git/cost-of-resistance/data/raw/RNAseq")
#setwd("/Volumes/TB_Research_Seq/Users/Andrej/analysis_logs/cost-of-resistance/data/raw/RNAseq/")

#Load libraries
library('DESeq2')
library('FactoMineR')

#Load TbX008 data and data structure
datafile8 = 'TbX008_dataframe.tsv'
design_file8 = 'TbX008_rnaseq.txt'
TbX_DataTable = read.table(datafile8, header=TRUE, row.names=1)
TbX_Design = read.table(design_file8, header=TRUE, row.names=1)

#Add evolved status
TbX_Design$evolved <- factor(c('no', 'no', 'no', 
                               'no', 'no', 'no', 
                               'yes', 'yes', 'yes',
                               'yes', 'yes', 'yes'))

#Extract only N0155-N1981 data
TbX_Design_evo = subset(TbX_Design, TbX_Design$evolved=="yes")
TbX_DataTable_evo = TbX_DataTable[,rownames(TbX_Design_evo)]

#Run DESeq
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable_evo,
                              colData = TbX_Design_evo,
                              design = ~ rif)

dds <- DESeq(dds)
res <- results(dds)	

evo_sig = table(res$padj<0.05)["TRUE"]

write.csv(res, file='../../interim/1_DESeq_TbX008_pair_N1588.csv')

#Extract only N0155-N1981 data
TbX_Design_N0155 = subset(TbX_Design, TbX_Design$evolved=="no")
TbX_DataTable_N0155 = TbX_DataTable[,rownames(TbX_Design_N0155)]

#Run DESeq
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable_N0155,
                              colData = TbX_Design_N0155,
                              design = ~ rif)

dds <- DESeq(dds)
res <- results(dds)	

N0155_sig = table(res$padj<0.05)["TRUE"]


write.csv(res, file='../../interim/1_DESeq_TbX008_pair_N0155.csv')


#Load data and data structure
datafile11 = 'TbX011_dataframe.tsv'
design_file11 = 'TbX011_rnaseq.txt'
TbX_DataTable = read.table(datafile11, header=TRUE, row.names=1)
TbX_Design = read.table(design_file11, header=TRUE, row.names=1)

#Do overall analysis
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable,
                              colData = TbX_Design,
                              design = ~ rif)

dds <- DESeq(dds)
res <- results(dds)	

write.csv(res, file='../../interim/1_DESeq_TbX011_all.csv')

#stratify by background
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable,
                              colData = TbX_Design,
                              design = ~ rif/background)

dds <- DESeq(dds)
res <- results(dds)	

write.csv(res, file='../../interim/1_DESeq_TbX011_all_slash_background.csv')

#Extract only N0157-N2030 data
TbX_Design_N0157 = subset(TbX_Design, TbX_Design$background=="N0157")
TbX_DataTable_N0157 = TbX_DataTable[,rownames(TbX_Design_N0157)]

#Run DESeq
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable_N0157,
                              colData = TbX_Design_N0157,
                              design = ~ rif)

dds <- DESeq(dds)
res <- results(dds)	

N0157_sig = table(res$padj<0.05)["TRUE"]

write.csv(res, file='../../interim/1_DESeq_TbX011_pair_N0157.csv')

#Extract only N0072-N2027 data
TbX_Design_N0072 = subset(TbX_Design, TbX_Design$background=="N0072")
TbX_DataTable_N0072 = TbX_DataTable[,rownames(TbX_Design_N0072)]

#Run DESeq
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable_N0072,
                              colData = TbX_Design_N0072,
                              design = ~ rif)

dds <- DESeq(dds)
res <- results(dds)	

N0072_sig = table(res$padj<0.05)["TRUE"]

write.csv(res, file='../../interim/1_DESeq_TbX011_pair_N0072.csv')

## SAMPLES WITH ONLY 2 REPLICATES PER STRAIN

#Extract only N0145-N1888 data from TbX011
TbX_Design_N0052 = subset(TbX_Design, TbX_Design$background=="N0052")
TbX_DataTable_N0052 = TbX_DataTable[,rownames(TbX_Design_N0052)]

#Run DESeq
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable_N0052,
                              colData = TbX_Design_N0052,
                              design = ~ rif)

dds <- DESeq(dds)
res <- results(dds)	

N0052_sig11 = table(res$padj<0.05)["TRUE"]

write.csv(res, file='../../interim/1_DESeq_TbX011_pair_N0052.csv')

#Extract only N0145-N1888 data from TbX011
TbX_Design_N0145 = subset(TbX_Design, TbX_Design$background=="N0145")
TbX_DataTable_N0145 = TbX_DataTable[,rownames(TbX_Design_N0145)]

#Run DESeq
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable_N0145,
                              colData = TbX_Design_N0145,
                              design = ~ rif)

dds <- DESeq(dds)
res <- results(dds)	

N0145_sig11 = table(res$padj<0.05)["TRUE"]

write.csv(res, file='../../interim/1_DESeq_TbX011_pair_N0145.csv')


#Extract only N0155-N1981 data from TbX011
TbX_Design_N0155 = subset(TbX_Design, TbX_Design$background=="N0155")
TbX_DataTable_N0155 = TbX_DataTable[,rownames(TbX_Design_N0155)]

#Run DESeq
dds <- DESeqDataSetFromMatrix(countData = TbX_DataTable_N0155,
                              colData = TbX_Design_N0155,
                              design = ~ rif)

dds <- DESeq(dds)
res <- results(dds)	

N0155_sig11 = table(res$padj<0.05)["TRUE"]

write.csv(res, file='../../interim/1_DESeq_TbX011_pair_N0155.csv')