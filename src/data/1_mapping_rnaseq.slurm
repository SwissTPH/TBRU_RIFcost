#!/bin/python

#SBATCH --job-name=rpoB_pathways
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

import sys
import os

with open(sys.argv[1]) as _f:
    _input = _f.read()

for line in _input.strip():
    info = line.strip().split(',')
    cmd = "./scicore/home/gagneux/trauner/cost-of-resistance/src/data/1_Stranded_Position_count.sh {} {} {}".format(info[0],info[1],info[2])
    os.system(cmd)