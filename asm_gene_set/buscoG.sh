#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 12

source activate /data/eande106/software/conda_envs/busco

mkdir -p busco

busco -i $file -c 12 -m genome -l /vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/busco_downloads/lineages/nematoda_odb10/ --out_path busco -c 12
