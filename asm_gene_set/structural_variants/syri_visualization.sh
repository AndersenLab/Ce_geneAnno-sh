#!/bin/bash

#SBATCH -J syri_vis
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -c 4

source activate plotsr

GENOMES="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/genomes/genomes_paths.tsv"
SYRI_OUT="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/output/CB4856_syrisyri.out"

cd /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/output
plotsr --sr $SYRI_OUT --genomes $GENOMES -o syri_vis_CB4856_nuclear.png -d 600 --chr I --chr II --chr III --chr IV --chr V --chr X
