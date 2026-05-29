#!/bin/bash

#SBATCH -J SV_merging                    # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 48                           # Number of cores

source activate jasmine

mkdir -p tmp/$SLURM_JOB_ID
mkdir -p output

TMP="tmp/$SLURM_JOB_ID"
input="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/elegans/strain_dirs/all_vcfs/filtered_vcfs/all_vcfs.tsv"

jasmine --output_genotypes threads=48 outdir=$TMP file_list=$input out_file=output/141_SVs_merged.vcf
