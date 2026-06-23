#!/bin/bash

#SBATCH -J Prot_filt
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -c 12

source activate MSA

# for file in /vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/CYTO_P450_SC_OGS/MSA/*.fa; do sbatch --export=file=$file <script>; done

OUT_DIR="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/CYTO_P450_SC_OGS/MSA/identity_matrices"

mkdir -p $OUT_DIR

og_name=$(basename $file .fa)

# Run clustalo to get percent identity matrix
clustalo -i $file --distmat-out=$OUT_DIR/${og_name}_identity.txt --percent-id --full 

# Convert to clean TSV
sed -E 's/[[:space:]]+/\t/g' $OUT_DIR/${og_name}_identity.txt > $OUT_DIR/${og_name}_identity.tsv
