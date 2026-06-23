#!/bin/bash

#SBATCH -J run_purge_dups
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 24

#source activate 2024LabRetreat #this environment contains minimap2

set -euo pipefail

execs=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/software/purge_dups

export PATH=$execs/scripts:$PATH

python $execs/scripts/run_purge_dups.py -p bash h_spumosa_inbred.config.json $execs/bin h_spumosa_inbred 
