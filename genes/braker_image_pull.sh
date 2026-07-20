#!/bin/bash

#SBATCH -J sing_pull
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=braker_image_pull.oe  # Output log file
#SBATCH --error=braker_image_pull.rr

echo "SCRIPT IS RUNNING! HERE IS A LOG MESSAGE"

module load singularity

cd /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/container_images

singularity pull --name loconn13999-braker3_20260720.sif docker://teambraker/braker3:v3.1.1
