#!/bin/bash

################################################################################
#
# Script was written for mapping of fastq RNAseq reads, and counting features using htseq.
#
# (Primitive) usage: ./script path/to/RNAseq.fastq sample.name reference.gff
#
# Fastq-files are expected to be single end, have a .fastq extension. This iteration will map to the ancestor genome
#
# Outputs count tables for fastq files based on specified gff.
#
# Andrej Trauner, 06. November 2014
#                 24. November 2017
#
################################################################################


########################### Definitions ########################################

######## Define filenames and reference paths:

module load Python/2.7.11-goolf-1.7.20
module load SAMtools/1.3.1-goolf-1.7.20
module load BWA/0.7.13-goolf-1.7.20
module load HTSeq/0.6.1p1-goolf-1.7.20-Python-2.7.11

SAMPLENAME=$3

#concatenate two gzipped fastq files.
cat $1 $2 > $TMPDIR/$SAMPLENAME.fastq.gz

FILENAME=$(basename $TMPDIR/$SAMPLENAME.fastq.gz) ### take everything from the script argument that comes after the slash
BASENAME=$(basename $1 .fastq) ### take everything from the script argument that comes after the slash and remove '.fastq'
PATHNAME=$(dirname $1) ### take everything from the script argument that comes before the slash
REFERENCE=/scicore/home/gagneux/trauner/ref/MTB_ancestor_reference.fasta
REFERENCEGFF=$4
REFERENCEGFFAS=$5

echo "processing $SAMPLENAME.fastq.gz"
echo "mapping $SAMPLENAME to $REFERENCE"


#bwa index -a bwtsw $REFERENCE


### The following command takes an SE fastq and aligns it to the index using BWA

bwa aln $REFERENCE "$TMPDIR/$SAMPLENAME.fastq.gz" > $TMPDIR/$SAMPLENAME.sai
bwa samse $REFERENCE $TMPDIR/$SAMPLENAME.sai "$TMPDIR/$SAMPLENAME.fastq.gz" > $TMPDIR/$SAMPLENAME.sam

echo "alignment of $SAMPLENAME done"

######### The following command take sam file and generate a sorted bam used for feature counting
echo "sorting $SAMPLENAME"

samtools view -bS $TMPDIR/$SAMPLENAME.sam > $TMPDIR/$SAMPLENAME.bam
samtools sort -n $TMPDIR/$SAMPLENAME.bam -o $TMPDIR/$SAMPLENAME.nsorted.bam
echo "sorting of $SAMPLENAME by name done"
#samtools sort $TMPDIR/$BASENAME.bam $TMPDIR/$BASENAME.csorted
#echo "sorting of $BASENAME by coordinate done"

######### The following command counts the features based on the gff file.
echo "counting features in $SAMPLENAME for sRNA"
#samtools view -h -F 4 $TMPDIR/$BASENAME.nsorted.bam | htseq-count -m intersection-nonempty -f sam $REFERENCEGFF > $SAMPLENAME.ancestor.htc
htseq-count -m intersection-nonempty -f bam $TMPDIR/$SAMPLENAME.nsorted.bam $REFERENCEGFF > "$SAMPLENAME"_sRNA.ancestor.htc
##NB THIS NEEDS NAME SORTED BAM

echo "counting features in $SAMPLENAME for antisense transcripts"
#samtools view -h -F 4 $TMPDIR/$BASENAME.nsorted.bam | htseq-count -m intersection-nonempty -f sam $REFERENCEGFF > $SAMPLENAME.ancestor.htc
htseq-count -m intersection-nonempty -f bam $TMPDIR/$SAMPLENAME.nsorted.bam $REFERENCEGFFAS > "$SAMPLENAME"_AS.ancestor.htc
##NB THIS NEEDS NAME SORTED BAM

echo "counting of $SAMPLENAME done"

######### The following generates a trace file.
#echo "generating pileup for $BASENAME"
#samtools mpileup $TMPDIR/$BASENAME.csorted.bam | awk 'BEGIN {FS="\t"}; {print $2,FS,$4}' - > $SAMPLENAME.ancestor.trace
##NB THIS NEEDS COORDINATE SORTED BAM
#echo "piling up of $BASENAME done
