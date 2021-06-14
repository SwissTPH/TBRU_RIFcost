#!/bin/python

#SBATCH --job-name=FtCount
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2:00:00
#SBATCH --output=/scicore/home/gagneux/trauner/logs/%j.%a.o
#SBATCH --error=/scicore/home/gagneux/trauner/logs/%j.%a.e
#SBATCH --qos=6hours

################################################################################
#
# Script for batch submission to SLURM using python
#
# Primitive usage: qbatch --array 1-18 launch__FtCount.py INPUT_FILE
#
# INPUT_FILE is a csv, with 3 columns:
# <PATH_FASTQ>,<ID>,<GFF>
#
# Andrej Trauner, 02. January 2017
#
################################################################################


import sys
import os

#get the environment variable "SLURM_ARRAY_TASK_ID" to
#allow for array job submission to work properly
slurm_task_id = os.environ.get('SLURM_ARRAY_TASK_ID', None)

#define script - note the space at the end.
_script = '/scicore/home/gagneux/trauner/git/cost-of-resistance/src/data/1_Stranded_Position_count.sh '

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
