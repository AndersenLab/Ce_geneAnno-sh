#!/bin/bash

#SBATCH -J syri
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 48

source activate syri

REF="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/genomes/N2.PRJNA13758.WS283.genome.fa"
QUERY="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/genomes/CB4856.PRJNA275000.WS272.chromIDsFixed.fa"
ALN="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/genomes/out_i90_l100.CB4856.coords"
DELTA="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/genomes/out_i90_l100.CB4856.delta"

cd /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/output

syri -c $ALN -F T -d $DELTA -r $REF -q $QUERY --prefix CB4856_syri
