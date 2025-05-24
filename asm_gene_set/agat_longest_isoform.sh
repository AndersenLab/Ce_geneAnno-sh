#!/bin/bash

#SBATCH -J isoform                    # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p express                     # Partition/Queue name
#SBATCH -t 00:30:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 2                           # Number of cores

source activate agat

out_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform"

agat_sp_keep_longest_isoform.pl  -f $file -o $out_dir/${file%.*}.longest.gff3
