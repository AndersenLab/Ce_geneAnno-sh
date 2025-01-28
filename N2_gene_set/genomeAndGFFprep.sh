#!/bin/bash

#SBATCH -J geneAnno                # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -n 12                           # Number of cores
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/initialRun.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/initialRun.rr 


#=================================

# Calling core, soft-core, rare, and private genes in Caenorhabditis elegans

#=================================
source activate gene_annotation

raw_data="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data"
processed_data="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data"

ref_genome="c_elegans.PRJNA13758.WS283.genome.fa.fai"
GFF="c_elegans.PRJNA13758.WS283.csq.gff3"

# Make 1kb windows of the reference genome
bedtools makewindows -g $raw_data/$ref_genome -w 1000 > $processed_data/N2.WS238.1kb_windows.bed

# Extract gene coordinates from the gff3
awk '$3 == "gene" {print $1, $4, $5, $9}' OFS="\t" $raw_data/$GFF | sed 's/ID=gene://;s/;.*//' > $processed_data/refGenesCoordinates.bed