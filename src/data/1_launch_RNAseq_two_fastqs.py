#!/bin/python

#SBATCH --job-name=mapping
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

################################################################################
#
# Script for batch submission to SLURM using python
#
# Primitive usage: sbatch --array 1-18 launch__RNAseq.py INPUT_FILE
#
# INPUT_FILE is a csv, with 3 columns:
# <PATH_FASTQ>,<ID>,<GFF>
#
# Andrej Trauner, 24. November 2017
#                 08. October 2019
#
################################################################################


import sys
import os

#get the environment variable "SLURM_ARRAY_TASK_ID" to
#allow for array job submission to work properly
slurm_task_id = os.environ.get('SLURM_ARRAY_TASK_ID', None)

#define script - note the space at the end.
#_script = '/scicore/home/gagneux/trauner/cost_of_dr/1_BWA_rnaseq_SE_ancestor_sRNA.sh '
_script = '/scicore/home/gagneux/trauner/cost-of-resistance/src/data/1_Stranded_Position_count_two_fastqs.sh '

#Read the input file.
with open(sys.argv[1]) as _input_file:
    f = _input_file.read()

#Submit the array job
if slurm_task_id:
    slurm_task_id = int(slurm_task_id)
    #Get the line of interest from the input list
    _input = f.split('\n')[slurm_task_id-1].strip()

    #Assemble command - notice the space as the string
    #separator
    _command = _script + ' '.join(_input.split(','))

    print('<INPUT COMMAND>:   {}'.format(_command))
    #run the command.
    os.system(_command)

#This might be uneccesary...
if not slurm_task_id:
    for _input in f:
        _input = _input.strip()
        _path, _id, _gff = _input.split(',')

        _command = _script + ' '.join(_input.split(','))

        os.system(_command)

#Clean up and exit gracefully.
sys.exit()
