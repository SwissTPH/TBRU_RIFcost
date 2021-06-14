#!/bin/sh

#SBATCH --job-name=AS_mapping
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --qos=1day
#SBATCH --output=/scicore/home/gagneux/trauner/logs/%j.%a.o
#SBATCH --error=/scicore/home/gagneux/trauner/logs/%j.%a.e
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=andrej.trauner@unibas.ch
#SBATCH --open-mode=append

module load Python/2.7.11-goolf-1.7.20
module load SAMtools/1.3.1-goolf-1.7.20
module load BWA/0.7.13-goolf-1.7.20
module load HTSeq/0.6.1p1-goolf-1.7.20-Python-2.7.11

FILENAME=$(basename $1) ### take everything from the script argument that comes after the slash
BASENAME=$(basename $1 .bam) ### take everything from the script argument that comes after the slash and remove '.bam'
PATHNAME=$(dirname $1) ### take everything from the script argument that comes before the slash
REFERENCE=/scicore/home/gagneux/trauner/ref/MTB_ancestor_reference.fasta
SAMPLENAME=$2
REFERENCEGFFAS=$3


echo "sorting $SAMPLENAME"
samtools sort -n $1 -o $TMPDIR/$SAMPLENAME.nsorted.bam
echo "sorting of $SAMPLENAME by name done"

echo "counting features in $SAMPLENAME for antisense transcripts"
htseq-count -m intersection-nonempty -f bam $TMPDIR/$SAMPLENAME.nsorted.bam $REFERENCEGFFAS > "$SAMPLENAME"_AS.ancestor.htc