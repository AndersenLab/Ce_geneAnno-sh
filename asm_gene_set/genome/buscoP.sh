#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --output=busco.oe
#SBATCH --job-name="busco"

source activate busco

busco -i $file -c 12 -m prot -l /vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/busco_downloads/lineages/nematoda_odb10/ -o /busco/proteome/$(basename $file .braker.protein.fa).busco -c 12
