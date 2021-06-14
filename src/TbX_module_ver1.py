##Analysis of RNAseq samples

import scipy as sp
from scipy import stats as ss
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
import json
from sklearn import utils as sku

mpl.rcParams['font.sans-serif']='Helvetica'
mpl.rcParams['legend.numpoints'] = 1

import warnings
warnings.filterwarnings("ignore")

class knowledge(dict):
    """Class containing a summary of genome information based on Broad annotation"""
    def __init__(self, filename, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self

        _all = []
        for line in open(filename):
           line=line.strip()
           split=line.split('\t')
           _all.append(split)

        _genes = []

        self.index = np.array(_all[0][1:])
        for value in _all[1:]:
           _genes.append(value[0])
           self[value[0]] = value[1:]
        self.genes = np.array(sorted(_genes))
        self.populated_genesets=[]

        self['mapper']={} #dictionary to hold geneset categories and their descriptions.

        try:
            for line in open('Geneset_mapper.txt'):
                split=line.split('\t')
                self.mapper[split[0].strip()]=split[1].strip()
        except:
            print('No mapper file, make sure "Geneset_mapper.txt" is in the root directory')

    def add_extras(self, filename='H37Rv_extra_knowledge.txt'):
        """Include extra information not stored in the basic knowledge file.
        """
        new_index = []
        new_knowledge = {}

        for line in open(filename): #this part parses the extra_knowledge file, mind the structure, and updates the existing knowledge to match it.
            line = line.strip()
            split = line.split('\t')
            new_fields = 0
            if split[0] == '>ADD':
                new_index=split[1:]
                new_fields+=len(split[1:])
                for existing_gene in self.genes:
                    new_knowledge[existing_gene] = ['',]*new_fields

            if split[0][:2] =='Rv' and split[0] in new_knowledge.keys(): #may need to make it if split[0]=!'>ADD' in the future
                ind = new_index.index(split[1])
                if new_knowledge[split[0]][ind]!='': new_knowledge[split[0]][ind]+=',' #makes sure that multiple entries in the same category are added separated as a comma, consistent with normal.
                new_knowledge[split[0]][ind]+=split[2]

        updated_index = list(self.index)+new_index #make new index list
        self.index=np.array(updated_index) #update knowledge index list
        for k,v in new_knowledge.items(): #update data stored in knowledge
            self[k]+=v

    def fetch(self, gene, attribute='GENOME ONTOLOGY'):
       """Fetch information for a gene.

       Arguments:
       ---------
       gene: gene_code of interest
       attribute: attribute of interest
       _examples_ |'SYMBOL'|'NAME'|'GENOME ONTOLOGY'|'ENZYME CODE'|'KEGG'|'COG'|
       full list at self.index

       Output:
       -------
       Ordered np.array of sorted individual entries
       """

       attribute_index = np.where(self.index==attribute)[0]
       _items = []
       if len(self[gene]) > attribute_index[0]: _items += sorted(self[gene][attribute_index[0]].split(','))

       return np.array(_items)


    def populate_list(self, attribute='GENOME ONTOLOGY', label='GO'):
       """Generate a non-redundant sorted list of all categories in attribute

       Arguments:
       ---------
       attribute: attribute of interest
       _examples_ |'SYMBOL'|'NAME'|'GENOME ONTOLOGY'|'ENZYME CODE'|'KEGG'|'COG'|
       full list at self.index
       label: used to name the resulting category i.e. self.label
       _suggested_ |'GO'|'KEGG'|'EC'|'PWAY'|

       Output:
       -------
       Relevant category added to the current knowledge instance.
       """

       _all = []
       for gene in self.genes:
           _all+=list(self.fetch(gene, attribute))

       self[label] = np.array(sorted(set(_all)))[1:]

    def populate_matrix(self, attribute='GENOME ONTOLOGY', label='GO'):
       """Generate a geneset matrix using populate_list output as base.

       Arguments:
       ---------
       attribute: attribute of interest
       _examples_ |'SYMBOL'|'NAME'|'GENOME ONTOLOGY'|'ENZYME CODE'|'KEGG'|'COG'|
       full list at self.index
       label: used to name the resulting category i.e. self.label
       _suggested_ |'GO'|'KEGG'|'EC'|'PWAY'|

       Output:
       -------
       Relevant category matrix added to the current knowledge instance.
       """

       self.populate_list(attribute, label)
       matrix_label = '%s_matrix' %label

       self[matrix_label]=np.zeros((len(self.genes),len(self[label]))) #make recipient matrix

       for index, gene in enumerate(self.genes):
           gene_ind = np.where(np.in1d(self[label], self.fetch(gene,attribute)))
           self[matrix_label][index][gene_ind]+=1

    def populate_geneset(self, geneset=['GENOME ONTOLOGY'], label=['GO']):
        """Populate selected curated genesets.

        Arguments:
        ----------
        geneset: any of the genesets specified in the Broad file
        _example_ 'GENOME ONTOLOGY'|'KEGG'|'COG'|PATHWAY|TB_PATHWAY|TB_REACTION|PFAM
        label: labels to be used for output matrices.

        Output:
        -------
        Addition of category lists (self.label) and category matrices (self.label_matrix)
        to the current instance of knowledge.

        Knowledge matrix notes: rows correspond to genes and columns to categories.
        sum on axis=0 gives the overall number of genes put in that category
        sum on axis=1 gives the number of categories per gene.

        """
        for index, _geneset in enumerate(geneset):
            self.populate_matrix(_geneset, label[index])
            self.populated_genesets.append(label[index])

class Bunch(dict):
    """Container object for datasets: dictionary-like object that
       exposes its keys as attributes."""

    def __init__(self, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self


    #Data import/export functions
    def import_genome_annotation(self, genome_info = 'H37Rv.json'):
        """Imports the genome annotation of an organism of interest"""

        try:
            details = json.load(open(genome_info,'r'))

            self['gene_info'] = details
            self['all_genes'] = sorted([x for x in details.keys() if x!='info'])

        except IOError:
            print('Cannot find and open {}.'.format(genome_info))

    def import_dataset(self, filename, N=10, format='htseq'):
        """Import experimental data formated as an 'experiment_index.txt' for 'htseq'
           or M table for 'microarray', or aFLQ.tsv for proteome and stored in the working folder"""
        #make N calculated internally.
        pos=0

        data_dictionary_temp = {}

        if format not in ['htseq','microarray', 'proteome']:
            print('Only HTSeq, microarray and proteome implemented so far. Reformat the data.')


        if format=='htseq':
            for line in open(filename):
                parameters = line.strip().split('\t')
                if parameters[0]=='strain':
                    self['categories'] = parameters
                    for parameter in self.categories:
                        self[parameter.strip()]=[]
                if parameters[0]!='strain':
                    _filename=parameters[-1].strip()
                    for pind, pval in enumerate(self.categories):
                        self[pval].append(parameters[pind])
                    #self['lineage'].append(parameters[1])
                    #self['isoniazid'].append(parameters[2])
                    #self['sequencer'].append(parameters[3])
                    #self['filename'].append(_filename)
                    for information in open(_filename):
                        measurement = information.split('\t')
                        if measurement[0][0:2]=='__': measurement[0]=measurement[0][2:]
                        if measurement[0] in data_dictionary_temp:
                            data_dictionary_temp[measurement[0]][pos]+=int(measurement[1])
                        if measurement[0] not in data_dictionary_temp:
                            data_dictionary_temp[measurement[0]]=np.zeros(N)
                            data_dictionary_temp[measurement[0]][pos]+=int(measurement[1])
                    pos+=1

            htseq_extras=['no_feature', 'ambiguous', 'too_low_aQual', 'not_aligned', 'alignment_not_unique']

            for extra in htseq_extras:
                self[extra]=data_dictionary_temp[extra]
                data_dictionary_temp.pop(extra)

            self['genes'] = np.array(sorted(data_dictionary_temp.keys()))

            _data=[]
            for gene in self.genes:
                _data.append(data_dictionary_temp[gene])

            self['data'] = np.array(_data)

            self.import_genome_annotation(genome_info = 'H37Rv.json')

            self['samples'] = [sam[:7] for sam in self.filename]

        if format == 'microarray':
            print('Checking for knowledge...')
            try:
                if self.knowledge.genes: print('Ready to import')
            except:
                print('Need to import knowledge...trying now')
            try:
                self.import_knowledge()
            except:
                print('Cannot import knowledge, make sure default file in folder. - Crapping out.')

            print('Done importing knowledge, moving on to the good stuff.')
            temp_data = np.zeros((len(self.knowledge.genes), int(N)+1))
            _genes = list(self.knowledge.genes)

            for line in open(filename):
                split = line.split('\r')
                for cols in split:
                    col_split = cols.split(',')
                    if col_split[0] == 'SystematicName':
                        self['strain'] = col_split[1:]
                    if col_split[0] in self.knowledge.genes:
                        values = np.array(col_split[1:], dtype=float)
                        values = np.append(values, 1)
                        index = _genes.index(col_split[0])
                        temp_data[index]+=values
            fix = np.where(temp_data[:,-1] == 0)[0] #This is here to avoid division by 0 for genes without probes.
            temp_data[:,-1][fix]+=1
            _data = temp_data.T/temp_data[:,-1]
            self.data = _data.T[:,:-1]
            self['genes']=_genes

        if format == 'proteome':

            _data = []
            _proteins_quantified = []

            #preliminary load of the knowledge.
            print('Trying to load the knowledge')
            try:
                self['knowledge'] = knowledge('mycobacterium_tuberculosis_h37rv_2_genome_summary_per_gene.txt')
            except:
                print('Cannot find the knowledge file, make sure it is in the directory')

            _genes = list(self.knowledge.genes)


            for line in open(filename):
                split = line.split('\r')
                for cols in split:
                    col_split = cols.split('\t')
                    if col_split[0] == 'Protein':
                        self['sample'] = col_split[1:]
                    if col_split[0] in self.knowledge.genes:
                        _data.append(col_split[1:])
                        _proteins_quantified.append(_genes.index(col_split[0]))

            self['data'] = np.array(_data, dtype=float)
            self['genes']= np.array(_genes)[_proteins_quantified]
            self['strain'] = [''.join(['N',sam[8:12]]) for sam in self.sample]

            self.import_knowledge(proteome=True)

    def import_knowledge(self, filename='mycobacterium_tuberculosis_h37rv_2_genome_summary_per_gene.txt', extras=True, proteome=False):
        """Import knowledge file."""

        self['knowledge']=knowledge(filename)

        if proteome: #trim the genelist to include only proteins for which I have data.
            _all_genes = self.knowledge.genes #get list of all proteins
            self.knowledge['genes'] = self.genes #re-define the genelist in knowledge before populating the matrices
            self.knowledge['all_genes'] = _all_genes #store the complete genelist in a new recipient.

        if extras==True:
            self.knowledge.add_extras()
            self.knowledge.populate_geneset(geneset=['GENOME ONTOLOGY', 'KEGG', 'COG', 'PATHWAY','TB_PATHWAY','TB_REACTION','PFAM','ENZYME CODE', 'REGION OF DIFFERENCE', 'REGULATOR', 'GENE CLASS', 'COEXPRESSION'], label=['GO', 'KEGG', 'COG', 'PWY','TBPWY','TBRXN','PF','EC', 'RD','REG','CL', 'CAT'])
        if extras!=True:
            self.knowledge.populate_geneset(geneset=['GENOME ONTOLOGY', 'KEGG', 'COG', 'PATHWAY','TB_PATHWAY','TB_REACTION','PFAM','ENZYME CODE'], label=['GO', 'KEGG', 'COG', 'PWY','TBPWY','TBRXN','PF','EC'])


    def export_data(self, output_format='variable', trim=False):
        """export loaded dataframe as a .tsv dataframe.

        Args:
        -----
        output_format: 'variable' will return the matrix to be stored locally,
        a filename will store the file as a file.
        trim: A list (or np.array) of indices to be excluded. E.g. np.arange(0,74) excludes all RNAs.

        Output:
        -----
        The matrix of measurements or filename.tsv stored in the IDE's root.
        """

        fg = np.vstack((self.genes,self.genes)).T

        output_dataframe = np.hstack((fg,self.data)).T[1:].T
        if trim is not False:
            fg = np.vstack((self.genes[trim],self.genes[trim])).T
            output_dataframe = np.hstack((fg,self.data[trim])).T[1:].T

        _temp_header = []
        for _category in self.categories[1:]: #so that it excludes the strains cateogry.
            _temp_header.append([_category]+list(self[_category]))

        header = np.vstack((_temp_header, ['genes']+list(self.samples)))

        output_dataframe = np.vstack((header,output_dataframe))

        if output_format == 'variable': return output_dataframe
        else: np.savetxt(output_format, output_dataframe, fmt='%s', delimiter=u'\t')


    def import_DEcsv(self, filename, p_cutoff=0.05, label='DESeq', DEtype='DESeq'):
        """Import DE processed sample.
        """

        index_name = '%s_index' %label
        genes_name = '%s_genes' %label
        data_name = '%s_data' %label
        sig_name = '%s_significant' %label
        non_name = '%s_nonsignificant' %label

        _DE_genes = []
        _DE_data=[]
        for line in open(filename):
            if DEtype=='MSstats':
                split2 = line.strip().split('\t') #not a csv made with excel
                if split2[0]=='Protein':
                    self[index_name]=np.array(split2[2:])
                if split2[0]!='Protein':
                    _DE_genes.append(split2[0])
                    _DE_data.append(split2[2:])
            if DEtype!='MSstats':
                split=line.split('\r')
                for cols in split:
                    cols_split=cols.split(',')
                    #Import from DESeq
                    if cols_split[1]=='id' and DEtype=='DESeq':
                        self[index_name]=np.array(cols_split[2:])
                    if cols_split[1]!='id' and DEtype=='DESeq':
                        _DE_genes.append(cols_split[1])
                        _DE_data.append(cols_split[2:])
                    #Import from edgeR
                    if cols_split[1]=='logFC' and DEtype=='edgeR':
                        self[index_name]=np.array(cols_split[1:])
                    if cols_split[1]!='logFC' and DEtype=='edgeR':
                        _DE_genes.append(cols_split[0])
                        _DE_data.append(cols_split[1:])
                    #Import from limma
                    if cols_split[1]=='logFC' and DEtype=='limma':
                        self[index_name]=np.array(cols_split[1:])
                        _DE_genes=list(self.knowledge.genes)
                        temp_data = np.zeros((len(self.knowledge.genes), 5)) #4 columns plus one counting the number of instances
                    if cols_split[1]!='logFC' and DEtype=='limma':
                        if cols_split[0] in _DE_genes:
                            values = np.array(cols_split[1:], dtype=float)
                            values = np.append(values, 1)
                            index = _DE_genes.index(cols_split[0])
                            temp_data[index]+=values
                    #Import from DESeq2
                    if cols_split[1]=='baseMean' and DEtype=='DESeq2':
                        self[index_name]=np.array(cols_split[1:])
                    if cols_split[1]!='baseMean' and DEtype=='DESeq2':
                        _DE_genes.append(cols_split[0])
                        _DE_data.append(cols_split[1:])


        _DE_data = np.array(_DE_data)
        if DEtype=='limma':
            fix = np.where(temp_data[:,-1] == 0)[0]
            temp_data[fix]+=np.array([0,1,1,1,1])
            _DE_data = (temp_data.T/temp_data[:,-1])[:-1].T
        NA = np.where(_DE_data=='NA')
        _DE_data[NA]=1
        _DE_data = np.array(_DE_data, dtype=float)
        self[genes_name] = np.array(_DE_genes)
        self[data_name] = _DE_data
        if DEtype=='DESeq': #Define significant and non-significant genes for DESeq
            self[sig_name] = np.where(_DE_data[:,6]<p_cutoff)[0]
            self[non_name] = np.where(_DE_data[:,6]>p_cutoff)[0]
        if DEtype=='edgeR': #Define significant and non-significant genes for edgeR
            self[sig_name] = np.where(_DE_data[:,3]<p_cutoff)[0]
            self[non_name] = np.where(_DE_data[:,3]>p_cutoff)[0]
        if DEtype=='DESeq2': #Define significant and non-significant genes for DESeq2
            self[sig_name] = np.where(_DE_data[:,5]<p_cutoff)[0]
            self[non_name] = np.where(_DE_data[:,5]>p_cutoff)[0]
        if DEtype=='MSstats': #Define significant and non-significant genes for MSstats
            self[sig_name] = np.where(_DE_data[:,5]<p_cutoff)[0]
            self[non_name] = np.where(_DE_data[:,5]>p_cutoff)[0]

    def descriptor(self, genes):
        """Find genesets that encompass most of the genes in a list"""

        if type(genes) is not np.ndarray: genes=np.array(genes)
        _query_indices = np.where(np.in1d(self.knowledge.genes, genes))

        N=float(len(genes))

        overlap90 = []
        any_overlap = []

        for geneset in self.knowledge.populated_genesets:
            for category in self.knowledge[geneset]:
                _indices = self.fetch_category(category, mapped=False)[0]
                _gene_indices = np.where(np.in1d(_query_indices, _indices))
                _genes = list(genes[_gene_indices])
                overlap = np.sum(np.in1d(_query_indices, _indices))/N
                if overlap != 0: any_overlap.append((overlap, category, _genes))
                if overlap >= 0.9: overlap90.append((overlap, category, _genes))

        overlap90 = sorted(overlap90, reverse=True)
        any_overlap = sorted(any_overlap, reverse=True)

        if len(overlap90) < 0:
                print('----------\nOverlapping genesets (90%+)\n----------\n%\tGeneset\tDescription\Genes covered')
                for (overlap, category, _genes) in overlap90:
                    if category[:2]=='PF': _category=category.split('.')[0]
                    else: _category=category
                    print('%.1f%%\t%s\t%s' %(100*overlap, category, self.knowledge.mapper.get(_category,''), _genes))

        if len(overlap90) >= 0:
                print('----------\nTop Overlapping genesets\n----------\n%\tGeneset\tDescription\tGenes covered')
                if len(any_overlap) < 15: K = len(any_overlap)
                else: K=15
                for (overlap, category, _genes) in any_overlap[:K]:
                    if category[:2]=='PF': _category=category.split('.')[0]
                    else: _category=category
                    print('%.1f%%\t%s\t%s\t%s' %(100*overlap, category, self.knowledge.mapper.get(_category,''), _genes))

    #Analysis
    def geneset_enrichment(self, method='ORA_fisher', dataset='DESeq', DEtype='DESeq', geneset='GO', category='all'):
        """Calculate gene set enrichment p-values.

        Arguments:
        ----------
        method: Statistical method for p-value calculation. Implemented methods:
            'ORA_fisher' perform an over-representation analysis based on Fisher's exact
            'ORA_chi2' perform an over-representation analysis based on Chi2 test
            'CERNO' Coincident Extreme Ranks in Numerical Observations adapted from (Kunnath-Velayudhan et al, 2010, PNAS 107(33):14703:14708).

        dataset: Data set to use for the analysis. It is important that the data
                 set is imported and that the label mathces that added to the Bunch().

        DEtype: Source of DE data, expects either 'DESeq' or 'edgeR'

        geneset: Gene set to use, it is important that the gene set has been
                 imported and added to the class prior to running the function.

        category: Calculate enrichment for particular category within the defined
                  geneset. Default set to 'all'.

        Output:
        -------
        np.array(p_values), the order reflect that of self.knowledge.geneset

        Notes:
        ------
        The Bunch() instance should have a .knowledge class attached and its
        genesets should be populated. A DE experiment should also be imported
        and its labels match that specified on top.
        """

        _dataset_genes = '%s_genes' %dataset
        _dataset_significant = '%s_significant' %dataset
        _dataset_matrix = '%s_data' %dataset

        _pvalindex, _foldindex = 5, 3
        if DEtype == 'DESeq': _pvalindex, _foldindex = 5, 3
        if DEtype == 'DESeq2': _pvalindex, _foldindex = 4, 1
        if DEtype == 'DESeqFold': _pvalindex, _foldindex = 3, 3
        if DEtype == 'edgeR': _pvalindex, _foldindex = 2, 0
        if DEtype == 'limma': _pvalindex, _foldindex = 2, 0

        #Map the significant genes from self.DE_significant onto the self.knowledge.genes. This is done so that the gene_specific indices
        #match those specified in the self.knowledge.geneset_matrix, as they may not always match.

        _mapped_significant_genes = np.where(np.in1d(self.knowledge.genes, self[_dataset_genes][self[_dataset_significant]]))

        _geneset_flavour = geneset
        _geneset_matrix = '%s_matrix' %geneset

        _n = len(_mapped_significant_genes[0]) #number of significant genes
        _ks = np.sum(self.knowledge[_geneset_matrix][_mapped_significant_genes],axis=0) # calculation of significant genes falling into a geneset
        _ms = np.sum(self.knowledge[_geneset_matrix], axis=0) # number of all the genes in a group
        _N = len(self.knowledge.genes) # number of genes in the genome

        _pval_ranks = np.argsort(self[_dataset_matrix][:,_pvalindex])*1.
        #_fold_change = self[_dataset_matrix][:,_foldindex]

        _geneset_pvals=[]
        #_geneset_meanchange=[]
        #_geneset_mapped_indices=[]

        if category=='all':
            for index,_category in enumerate(self.knowledge[_geneset_flavour]):
                _category_gene_ranks = np.where(np.in1d(self[_dataset_genes], np.array(self.knowledge.genes)[np.where(self.knowledge[_geneset_matrix][:,index]==1)])) #get the indices of genes with a certain category attached, note these are mapped to the self.genes not self.knowledge.genes
                # I think this is wrong: contingency_table=np.array([[_ks[index],_ms[index]-_ks[index]],[_n-_ks[index],_N+_ks[index]-_n-_ms[index]]])
                contingency_table=np.array([[_ks[index],_ms[index]],[_n-_ks[index],_N-_ms[index]]])

                #_mean_change = np.mean(_fold_change[_category_gene_ranks])
                #_geneset_meanchange.append(_mean_change)
                #_geneset_mapped_indices.append(_category_gene_ranks)

                if method=='ORA_fisher':
                   odds, pval = ss.fisher_exact(contingency_table)
                   _geneset_pvals.append(pval)

                if method=='CERNO':
                    _S = -2*sum(np.log(_pval_ranks[_category_gene_ranks]/_N))
                    df = len(_category_gene_ranks[0])
                    pval = ss.chi2.sf(_S, df)
                    _geneset_pvals.append(pval)

        if category!='all':
            index = np.where(self.knowledge[_geneset_flavour]==category)[0][0]

            _category_gene_ranks = np.where(np.in1d(self[_dataset_genes], np.array(self.knowledge.genes)[np.where(self.knowledge[_geneset_matrix][:,index]==1)]))

            contingency_table=np.array([[_ks[index],_ms[index]],[_n-_ks[index],_N-_ms[index]]]) #changed as above.

            if method=='ORA_fisher':
               odds, pval = ss.fisher_exact(contingency_table)
               _geneset_pvals.append(pval)

            if method=='CERNO':
                _S = -2*sum(np.log(_pval_ranks[_category_gene_ranks]/_N))
                df = len(_category_gene_ranks[0])
                pval = ss.chi2.sf(_S, df)
                _geneset_pvals.append(pval)

        return np.array(_geneset_pvals)

    def FDR_cutoff(self, pvals, q_cutoff=0.05):
        """Determine position of FDR cutoff

        Arguments:
        ----------
        pvals: list/np.array of pvals

        q_cutoff: FDR you are willing to tolerate.

        Output:
        -------
        index at which H0 is accpeted.

        """
        try:
            _pvals = np.array(pvals, dtype=float)
        except:
            print('Cannot convert pvals input into np.array')

        qvals=[]
        click = False

        ranks = np.argsort(_pvals)

        for index, value in enumerate(_pvals[ranks]):
            q = len(_pvals)*value/(index+1)
            qvals.append(q)
            if q>=q_cutoff and click==False: click=index

        return click

    def knn_search_rv(self, x, D, K):
        """Find K nearest neighbours of x among D
        x: Rv in self.genes
        D: data
        K: Number of neighbours"""
        _gene_index = list(self.genes).index(x)

        #Euclidean distances from the other points
        sqd = np.sqrt(np.sum((D-D[_gene_index])**2, axis=1))
        idx = np.argsort(sqd) # sorting
        #return the indexes of K nearest neighbours
        return idx[:K]

    def collate_genesets(self, method='CERNO', dataset='DESeq', DEtype='DESeq', FDR_cutoff=0.05):
        """Generate a sorted list of significant geneset categories.

        Arguments:
        ----------
        geneset: By default happens for all genesets in self.knowledge.populated_genesets.
        method: 'CERNO'|'ORA_fisher'|'ORA_chi2' method for geneset_enrichment
        dataset: dataset of interest for geneset_enrichment
        DEtype: 'DESeq'|'edgeR' dataset type for geneset_enrichment
        FDR_cutoff: False positive rate you are willing to tolerate.

        Output:
        -------
        sorted list of tuples [(pval, categoy),...]

        Notes:
        ------
        Make sure that all the relevant data are imported and mounted.
        """

        _significant_genesets = []

        for _geneset in self.knowledge.populated_genesets:

            _geneset_pvals = self.geneset_enrichment(method, dataset, DEtype, _geneset) #calculate pvals for all categories in a geneset group
            _FDR_index = self.FDR_cutoff(_geneset_pvals, q_cutoff = FDR_cutoff) #determine the FDR cutoff using the Benjamini-Hochberg approach.

            if _FDR_index>0: #make sure that only those genestes for which FDR makes the cut are included.
                _FDR_ranked = np.argsort(_geneset_pvals)[:_FDR_index]
                current = zip(_geneset_pvals[_FDR_ranked], self.knowledge[_geneset][_FDR_ranked])
                _significant_genesets+=current

        self['significant_genesets']=sorted(_significant_genesets)

    def fetch_category(self, geneset_category, dataset='DESeq', mapped=True):
        """Generate a np.array of all the gene indices pertinent to a geneset category.

        Arguments:
        ----------
        geneset: specific category of a geneset in any of the genesets: e.g. KEGG'|'GO'|'COG'|'PWY'
                e.g. PWY: ALADEG-PWY. The first part of the input will be used to identify the geneset.
        mapped: True|False self.genes and self.knowledge.genes are overallping but not identical
                this argument specifies whether the outputed indices should refer to self.genes
                (mapped = True) or self.knowledge.genes (mapped = False).

        Output:
        -------
        np.array(gene_indices)

        Notes:
        ------
        Requires the geneset matrix to be populated in knowledge.

        """

        geneset = geneset_category.split(':')[0] #identify geneset from category name
        _geneset_column = np.where(self.knowledge[geneset]==geneset_category)[0][0] #identify the column number in the populated geneset index file.
        _geneset_matrix = '%s_matrix' %geneset
        _dataset_genes = '%s_genes' %dataset

        _unmapped_geneset_indices = np.where(self.knowledge[_geneset_matrix][:,_geneset_column])
        _mapped_geneset_indices = np.where(np.in1d(self[_dataset_genes], self.knowledge.genes[_unmapped_geneset_indices]))

        if mapped==True: return _mapped_geneset_indices
        if mapped==False: return _unmapped_geneset_indices


    #Plots
    def correlation_plot(self, replicates=3, normalised=None, output_data=False, cmap='seismic'):
        """Plots a chequerboard plot of sample or gene Spearman correlations

        replicates - number of replicates per sample
        library - barcode library to be used, if available
        normalised - normalise the data and fit to specific range. Correct usage: normalised=Normalize(vmin=M, vmax=X) where M and X are the upper and lower bound.
        output_data - set to True if you want to save the plotted matrix
        cmap - name of colormap as defined by matplotlib.pyplot.colormaps()
        control - control sample
        treatment - if unspecified all the instances will be plotted, if not only the one of interest will be added.
        filtered, proportion - filtering parameters.
        """

        _sample_names = []
        _allnames = []
        for index, strain in enumerate(self.strain):
            sample = '%s_%s' %(self.strain[index], index+1)
            _allnames.append(sample)
            if sample not in _sample_names: _sample_names.append(sample)


        r_correlation, pvals = ss.spearmanr(self.data)

        plt.imshow(r_correlation, cmap, norm=normalised, interpolation='none')
        if replicates!=1:
            plt.yticks(np.arange(replicates/2,len(self.strain),replicates),_sample_names)
            plt.xticks(np.arange(replicates/2,len(self.strain),replicates),_sample_names)
        if replicates==1:
            plt.yticks(np.arange(replicates/2,len(self.strain),replicates),_allnames)
            plt.xticks(np.arange(replicates/2,len(self.strain),replicates),_allnames, rotation=45)
        plt.xlabel('Strain',size=18)
        plt.ylabel('Strain', size=18)
        plt.suptitle('Sample correlation matrix', size=24)
        if replicates!=1:
            plt.vlines(np.arange(replicates, len(self.strain), replicates)-0.5,-0.5,len(self.strain)-.5)
            plt.hlines(np.arange(replicates, len(self.strain), replicates)-0.5,-0.5,len(self.strain)-.5)
        plt.colorbar()
        plt.show()

    def scatterplots(self, x=[0], y=[1], mode='sum', scale='log'):
        """Plot a scatterplot of data.

        Arguments:
        ----------
        x: a list of indices to be considered, len of list must be 1 or same as y, unless mode is sum in which case 1+
        y: a list of indices to be considered, len of list must be 1 or same as x, unless mode is sum in which case 1+
        mode ['sum'|'onex_ally'|'oney_allx'|'paired']: determines what is plotted. 'sum' will sum across samples, one_all will plot all samples against first of the other, 'paired' requires x and y to have the same length plot paired
        scale ['log', 'linear']: determines the scale of the axis. It will be scaled as min-max
        Output:
        -------
        Scatterplot
        """
        xy=x+y
        x=np.array(x)
        y=np.array(y)
        xy=np.array(xy)

        xlen = len(x)
        ylen = len(y)

        color=['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black']

        if mode=='sum':
            sample_label = '%s vs %s, sum' % (self.strain[x[0]], self.strain[y[0]])
            plt.scatter(np.sum(self.data[:,x],axis=1), np.sum(self.data[:,y],axis=1), s=18, marker='o',label=sample_label, facecolors='none', edgecolors='black', alpha=0.2)
            plt.xlabel('%s counts' % self.strain[x[0]], size=18)
            plt.ylabel('%s counts' % self.strain[y[0]], size=18)
            plt.suptitle(sample_label, size=24)


        if mode=='onex_ally':
            x1 = np.repeat(x[0], ylen)

            for i,v in enumerate(y):
                sample_label = '%s vs %s_%s' % (self.strain[x1[i]], self.strain[y[i]], i)
                plt.scatter(self.data[:,x1[i]], self.data[:,y[i]], s=18, marker='o',label=sample_label, facecolors='none', edgecolors=color[i], alpha=0.2)

            plt.xlabel('%s counts' % self.strain[x1[0]], size=18)
            plt.ylabel('%s counts' % self.strain[y[0]], size=18)
            plt.suptitle('%s vs %s' % (self.strain[x1[0]], self.strain[y[0]]), size=24)
            plt.legend()


        if mode=='oney_allx':
            y1 = np.repeat(y[0], xlen)

            for i,v in enumerate(x):
                sample_label = '%s_%s vs %s' % (self.strain[x[i]],i ,self.strain[y1[i]])
                plt.scatter(self.data[:,x[i]], self.data[:,y1[i]], s=18, marker='o',label=sample_label, facecolors='none', edgecolors=color[i], alpha=0.2)

            plt.xlabel('%s counts' % self.strain[x[0]], size=18)
            plt.ylabel('%s counts' % self.strain[y1[0]], size=18)
            plt.suptitle('%s vs %s' % (self.strain[x[0]], self.strain[y1[0]]), size=24)
            plt.legend()


        if mode=='paired':
            for i,v in enumerate(y):
                sample_label = '%s_%s vs %s_%s' % (self.strain[x[i]], i, self.strain[y[i]], i)
                plt.scatter(self.data[:,x[i]], self.data[:,y[i]], s=18, marker='o',label=sample_label, facecolors='none', edgecolors=color[i], alpha=0.2)

            plt.xlabel('%s counts' % self.strain[x[0]], size=18)
            plt.ylabel('%s counts' % self.strain[y[0]], size=18)
            plt.suptitle('%s vs %s' % (self.strain[x[0]], self.strain[y[0]]), size=24)
            plt.legend()

        if scale=='log':
            sample_max=np.max(self.data[:,xy])
            plt.xscale('log')
            plt.yscale('log')
            plt.xlim(1,sample_max)
            plt.ylim(1,sample_max)

        plt.show()
    def boxplot(self, dataset='data'):
        """Plot boxplot for all samples from selected data table.

        dataset: usually self.data, but can be anything calculated and attached to the Bunch.
        """

        plt.boxplot(self[dataset])
        plt.yscale('log')
        plt.xticks(np.arange(1,len(self.strain)+1,1), self.strain, rotation=45)
        plt.xlabel('Sample', fontsize=18)
        plt.ylabel('Raw counts', fontsize=18)
        plt.suptitle('Mapped reads', fontsize=24)
        plt.show()

    def plotMA(self, dataset='DESeq'):
        """Generate MA plot marking significant and non-significant genes

        Arguments:
        ----------
        dataset: source of data to be plotted.

        Ouput:
        ------
        MA plot

        Notes:
        ------
        A DE analysis output must be attached to the Bunch() instance.
        """

        try:
            _data = '%s_data' %dataset
            _nonsignificant = '%s_nonsignificant' %dataset
            _significant = '%s_significant' %dataset

            plt.plot(self[_data][self[_nonsignificant]][:,0],self[_data][self[_nonsignificant]][:,3], 'ko', markerfacecolor='white', alpha=0.6)
            plt.plot(self[_data][self[_significant]][:,0],self[_data][self[_significant]][:,3], 'ko', markerfacecolor='red', alpha=0.9)
            plt.xscale('log')
            plt.xlim(1,10000000)
            plt.yscale('log')
            plt.ylim(0.0001, 10000)
            plt.hlines(1,0.01, max(self[_data][:,0]), color='red', linewidth=4, alpha=0.4)
            plt.xlabel('Normalised transcript counts', fontsize=18)
            plt.ylabel('Abundance ratio', fontsize=18)
            plt.suptitle('%s-analysed data - MA plot' %dataset, fontsize=24)

            plt.show()

        except:
            print('Could not generate MA plot, make sure the relevant dataset is correctly imported: self.import_DEcsv .')

    def GSplot(self, N=20, dataset = 'DESeq',DEtype='DESeq', geneset='all', return_genes=False):
        """Plot Geneset enrichments, covering N most significant geneset categories.

        Arguments:
        ----------
        N: number of geneset categories to plot
        dataset: DE dataset to be considered.
        DEtype: 'DESeq'|'edgeR'
        geneset: 'all'|'GO'|'GO:0000104' specify geneset group/category to plot
        return_genes: False|True printthe relative abundance per group, fold difference and adjusted p-values for significant genes in a set.

        Output:
        -------
        Geneset enrichment plot, tab delimited list of significan genes if desired.

        Notes:
        ------

        """

        _foldindex = 4
        if DEtype == 'DESeq': _foldindex = 4
        if DEtype == 'DESeq2': _foldindex = 1
        if DEtype == 'edgeR': _foldindex = 0
        if DEtype == 'limma': _foldindex = 0
        _dataset_significant = '%s_significant' %dataset
        _dataset_matrix = '%s_data' %dataset
        _dataset_genes = '%s_genes' %dataset

        if geneset=='all':
            geneset_list = self.significant_genesets[:N]
            if len(geneset_list)<N: N=len(geneset_list)

        if geneset!='all':
            geneset_list=[]
            for (p,category) in self.significant_genesets:
                if geneset in category: geneset_list.append((p,category))
            if N > len(geneset_list) and len(geneset_list)>0:
                N=len(geneset_list) #makes sure we don't go out of range
            if len(geneset_list)==0: #get data if a particular geneset is not in significant list.
                _gs = geneset.split(':')[0]
                cats = [x for x in self.knowledge[_gs] if geneset in x]
                N=len(cats)
                for cat in cats: geneset_list.append((1.0, cat))

        fig, ax = plt.subplots(N, sharex=True) #makes the backbone - N subplots in a single column
        if N>1: axt = np.array([x.twinx() for x in ax]) #calls all the second axes in order to split the label.
        if N==1: axt=ax.twinx()


        for index in range(N):
            (p,category)=geneset_list[index]

            if category[:2]=='PF': _category=category.split('.')[0] #hack to get around the fact that the Broad made up their categories from PFAM.
            else: _category=category

            description = self.knowledge.mapper.get(_category,'') #This will get the description for a geneset or nothing if not in the list.
            _indices = self.fetch_category(category, dataset, mapped=True)[0]
            all_genes = len(_indices)
            _indices1 = np.array(np.intersect1d(self[_dataset_significant], _indices))
            sig_genes = len(_indices1)
            bins = np.arange(-5,5.5,0.5)
            _hist, ban = sp.histogram (self[_dataset_matrix][:,_foldindex][_indices], bins)
            _hist_s, ban_s = sp.histogram(self[_dataset_matrix][:,_foldindex][_indices1], bins)
            sample_legendL = '%s/%s, pval = %.2E' %(sig_genes, all_genes, p)
            sample_legendR = '%s, %s' %(category, description)

            if N > 1:
                ax[index].bar(bins[:-1], _hist*1./max(_hist), width=0.5, color='black', alpha=0.2)
                axt[index].bar(bins[:-1], _hist_s*1./max(_hist), width=0.5, color='red', alpha=0.8)

                ax[index].set_yticks([])
                ax[index].set_ylim(0,1.)
                ax[index].set_ylabel(sample_legendL, rotation='horizontal', horizontalalignment='right')
                axt[index].set_yticks([])
                axt[index].set_ylim(0,1.)
                axt[index].set_ylabel(sample_legendR, rotation='horizontal', horizontalalignment='left')
                axt[index].vlines(0,0,1, 'r', linewidth=2)
                title = 'GSE-plot, %s, geneset group - %s' %(dataset, geneset)


            if N == 1:
                ax.bar(bins[:-1], _hist*1./max(_hist), width=0.5, color='black', alpha=0.2)
                axt.bar(bins[:-1], _hist_s*1./max(_hist), width=0.5, color='red', alpha=0.8)

                ax.set_yticks([])
                ax.set_ylim(0,1.)
                ax.set_ylabel(sample_legendL, rotation='horizontal',horizontalalignment='right')
                axt.set_yticks([])
                axt.set_ylim(0,1.)
                axt.set_ylabel(sample_legendR, rotation='horizontal',horizontalalignment='left')
                axt.vlines(0,0,1, 'r', linewidth=2)
                title = 'GSE-plot, %s, geneset - %s' %(dataset, geneset)



        plt.subplots_adjust(left=0.2, right=0.45, hspace=0)
        plt.xticks(np.arange(-10,10.5,5.))
        plt.xlim(-10,10)

        plt.suptitle(title, size=24)
        plt.show()

        if return_genes==True:
            _dataset_index='%s_index' %dataset
            if DEtype == 'DESeq':
                print('----------\n%s\n----------\nGene\t%s\t%s\t%s\t%s' %(geneset,self[_dataset_index][1], self[_dataset_index][2], self[_dataset_index][4], self[_dataset_index][6]))
                for ind in _indices:
                    print('%s\t%s\t%s\t%s\t%s' %(self[_dataset_genes][ind], self[_dataset_matrix][ind,1], self[_dataset_matrix][ind,2], self[_dataset_matrix][ind,4], self[_dataset_matrix][ind,6]))
            if DEtype == 'edgeR' or DEtype == 'limma':
                print('----------\n%s\n----------\nGene\t%s\t%s\t%s' %(geneset,self[_dataset_index][0], self[_dataset_index][2], self[_dataset_index][3]))
                for ind in _indices:
                    print('%s\t%s\t%s\t%s' %(self[_dataset_genes][ind], self[_dataset_matrix][ind,0], self[_dataset_matrix][ind,2], self[_dataset_matrix][ind,3]))


    def count_plot(self, genes, adjust=True, return_genes=False, category=None, cmap='seismic', normalised=None, figsize=(12,12)):
        """Plot counts for genes of interest.

        Arguments:
        ----------
        genes: one or more genes
        return_genes: False|True printthe data as well using show_DEsummary()
        cmap - name of colormap as defined by matplotlib.pyplot.colormaps()
        normalised - normalise the data and fit to specific range. Correct usage: normalised=Normalize(vmin=M, vmax=X) where M and X are the upper and lower bound.

        Output:
        -------
        Heatmap of selected genes

        Notes:
        ------

        """

        dats = self.data

        if adjust==True:
            dats = self.data.T/np.sum(self.data,axis=1)
            dats = dats.T

        if type(genes) is not np.ndarray: genes=np.array(genes)

        if return_genes==True: self.show_DEsummary(genes)

        _indices = np.where(np.in1d(self.genes,genes))
        if category != None:
            _indices2 = self.fetch_category(category, dataset, mapped=True)[0]+73
            genes2 = self.genes[_indices2]

        plt.figure('Genes', figsize=figsize)
        plt.imshow(dats[_indices], cmap, norm=normalised, interpolation='none')
        plt.yticks(np.arange(0,len(genes)),self.genes[_indices])
        plt.xticks(np.arange(0,len(self.strain)),self.strain, rotation=45)
        plt.xlabel('Sample', size=18)
        plt.ylabel('Gene', size=18)
        plt.suptitle('Gene counts', size=24)
        plt.colorbar()
        plt.subplots_adjust(left=0.3, right=0.7, hspace=0)
        plt.show()

        if category!=None:
            plt.figure('Category', figsize=figsize)
            plt.imshow(dats[_indices2], cmap, norm=normalised, interpolation='none')
            plt.yticks(np.arange(0,len(genes2)),genes2)
            plt.xticks(np.arange(0,len(self.strain)),self.strain, rotation=45)
            plt.xlabel('Sample', size=18)
            plt.ylabel('Gene',size=18)
            plt.suptitle('Gene counts for %s' %category, size=24)
            plt.colorbar()
            plt.subplots_adjust(left=0.3, right=0.7, hspace=0)
            plt.show()

    def count_plot2(self, genes, return_genes=False, category=None, cmap='seismic', normalised=None, figsize=(12,12)):
        """Plot counts for genes of interest.

        Arguments:
        ----------
        genes: one or more genes
        return_genes: False|True printthe data as well using show_DEsummary()
        cmap - name of colormap as defined by matplotlib.pyplot.colormaps()
        normalised - normalise the data and fit to specific range. Correct usage: normalised=Normalize(vmin=M, vmax=X) where M and X are the upper and lower bound.

        Output:
        -------
        Heatmap of selected genes

        Notes:
        ------

        """

        if type(genes) is not np.ndarray: genes=np.array(genes)

        if return_genes==True: self.show_DEsummary(genes)

        _indices = np.where(np.in1d(np.array(self.genes),genes))
        #print_indices
        if category != None:
            _indices2 = self.fetch_category(category, mapped=True)[0]
            genes2 = np.array(self.genes)[_indices2]

        plt.figure('Genes', figsize=figsize)
        plt.imshow(self.data[_indices], cmap, norm=normalised, interpolation='none')
        plt.yticks(np.arange(0,len(genes)),np.array(self.genes)[_indices])
        plt.xticks(np.arange(0,len(self.strain)),self.strain, rotation=45)
        plt.xlabel('Sample', size=18)
        plt.ylabel('Gene', size=18)
        plt.suptitle('Gene counts', size=24)
        plt.colorbar()
        plt.subplots_adjust(left=0.3, right=0.7, hspace=0)
        plt.show()

        if category!=None:
            plt.figure('Category', figsize=figsize)
            plt.imshow(self.data[_indices2], cmap, norm=normalised, interpolation='none')
            plt.yticks(np.arange(0,len(genes2)),genes2)
            plt.xticks(np.arange(0,len(self.strain)),self.strain, rotation=45)
            plt.xlabel('Sample', size=18)
            plt.ylabel('Gene', size=18)
            plt.suptitle('Gene counts for %s' %category, size=24)
            plt.colorbar()
            plt.subplots_adjust(left=0.3, right=0.7, hspace=0)
            plt.show()

    #Misc
    def show_DEsummary(self, genes, dataset='DESeq', DEtype='DESeq'):
        """Shows summary for genes of interest

        Arguments:
        ----------
        genes: one or more genes
        dataset: specifies dataset to use
        DEtype: 'DESeq'|'edgeR'

        Output:
        -------
        printout of adjusted counts, fold change and adjusted p-value (DESeq).
        printout fold change, p-value and FDR-adjuted p-value (edgeR).

        Notes:
        ------
        Make sure to use quotation marks and if using more than one gene use square brackets
        """

        if type(genes) is not np.ndarray: genes=np.array(genes)
        _dataset_matrix = '%s_data' %dataset
        _dataset_genes = '%s_genes' %dataset
        _dataset_index = '%s_index' %dataset

        _indices = [gen_ind for (gen_ind, value) in enumerate(self[_dataset_genes]) if value in genes]

        if DEtype=='DESeq':
            print('----------\nGenes of interest\n----------\nGene\t%s\t%s\t%s\t%s' %(self[_dataset_index][1], self[_dataset_index][2], self[_dataset_index][4], self[_dataset_index][6]))
            for ind in _indices:
                print('%s\t%s\t%s\t%s\t%s' %(self[_dataset_genes][ind], self[_dataset_matrix][ind,1], self[_dataset_matrix][ind,2], self[_dataset_matrix][ind,4], self[_dataset_matrix][ind,6]))

        if DEtype=='edgeR':
            print('----------\nGenes of interest\n----------\nGene\t%s\t%s\t%s' %(self[_dataset_index][0], self[_dataset_index][2], self[_dataset_index][3]))
            for ind in _indices:
                print('%s\t%s\t%s\t%s' %(self[_dataset_genes][ind], self[_dataset_matrix][ind,0], self[_dataset_matrix][ind,2], self[_dataset_matrix][ind,3]))


## TO DO
# import_data: add a sample_name category, fix plots to reflect samples.
# GSEA using Chi, Bootstrap maybe also. Find a way to add directionality to the data. Import mapping of categories as part of knowledge. Think about visualisations. Set up a single command pipeline.
# Fix scatterplot.
# export_data: add trim function, to delete all RD for given lineages and all structural RNAs.
# knowledge: the following gene entries were modified by placing them all into one gene operons. These facts were not expermentally verified, merely added to patch a gap in the reference file. Genes: ['Rv0157A', 'Rv0724A', 'Rv2023A', 'Rv3769', 'Rv3324A', 'Rv0469', 'Rv0491', 'Rv1638', 'Rv2143', 'Rv2529', 'Rv2917', 'Rv3219']
# add DEtype to plotMA and count_plot
