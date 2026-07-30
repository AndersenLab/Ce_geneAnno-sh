#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --output=busco.oe
#SBATCH --job-name="busco"

source activate busco

busco -i $file -c 12 -m prot -l /vast/eande106/data/DBs/BUSCO/nematoda_odb10/ --out_path <> -o $(basename $file .braker.protein.fa).busco
