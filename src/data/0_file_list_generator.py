#!/bin/python

########
# snipett to generate input file for pipeline
########

import sys

#identify experiemnt
experiment = str(sys.argv[1])

#load rnaseq.sorted.list -> sorted list of paths to rnaseq fastqs
with open(sys.argv[2]) as _input:
    f = _input.read()

#identify queries: RS-numbers pertinent to an experiment
#and the number of input fastqs.
if experiment in ['TbX007', 'tbx7', 'TbX7']:
    queries = ['RS010{}'.format(x) for x in range(18,36)]
    inputs = 1

if experiment in ['TbX008', 'tbx8', 'TbX8']:
    queries = ['RS010{}'.format(x) for x in range(36,51)]
    inputs = 1

if experiment in ['TbX011', 'tbx11', 'TbX11']:
    queries = ['RS010{}'.format(x) for x in range(51,75)]
    inputs = 2

else:
    print('Invalid experiment: TbX007|TbX008|TbX011')
    sys.exit(1)

#generate an info dictionary

info = {x: {0:'',
            1:'',
            'QUERY':x,
            'GFF': '/scicore/home/gagneux/trauner/ref/gff/MTB_ancestor_sRNA_collated.gff' ,
            'GFF_AS':'/scicore/home/gagneux/trauner/ref/gff/MTB_ancestor_reference_exonAS.gff'}
        for x in queries}

#populate the info dictionary
for x in queries:
    pos=0
    for line in f.split('\n'):
        if x in line and 'fastq.gz' in line:
            info[x][pos]=line.strip()
            pos=1

#write info dictionary to file
output_file = open('{}_rnaseq_sRNA.csv'.format(experiment),'w')

if inputs==2:
    for x in queries:
        line = '{fastq1},{fastq2},{query},{gff},{gff2}\n'.format(
            fastq1=info[x][0],
            fastq2=info[x][1],
            query=x,
            gff=info[x]['GFF'],
            gff2=info[x]['GFF_AS'])
        output_file.write(line)

if inputs==1:
    for x in queries:
        line = '{fastq1},{query},{gff},{gff2}\n'.format(
            fastq1=info[x][0],
            query=x,
            gff=info[x]['GFF'],
            gff2=info[x]['GFF_AS'])
        output_file.write(line)

output_file.close()
