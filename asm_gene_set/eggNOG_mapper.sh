#!/bin/bash

#SBATCH -J eggnog                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -n 24                            # Number of cores
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/scripts/SLURM_output/eggNOG_euk.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/scripts/SLURM_output/eggNOG_euk.rr

source activate eggnog

emapper.py -m diamond -i /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein.fa \
	--data_dir /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/software/eggNOG_mapper/databases \
	--seed_ortholog_evalue 0.001 --tax_scope eukaryota \
	--output_dir /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/output \
	--temp_dir /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/temp \
	-o N2_longestIso_background
