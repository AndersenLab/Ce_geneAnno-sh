#!/bin/bash

#SBATCH -J pav
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 48

module load singularity

container_image="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/container_images/pav-becklab-latest-20250209.sif"
analysis_dir="${PWD}"
pav_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/software/pav"

singularity run --bind $analysis_dir:/analysis_dir \
    --bind $pav_dir:/pav_dir \
    $container_image -c 48 

    