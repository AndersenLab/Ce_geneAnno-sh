#!/bin/bash

#SBATCH -J pavPREP
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pavprep.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pavprep.rr 

pav="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav"
assemblies="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/assemblies"
SLURM="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output"


# tail -n +2 $pav/assemblies.tsv | while IFS=$'\t' read -r strain path; do
#     mkdir -p $pav/strain_dirs/$strain
#     cp $pav/config.json  $pav/strain_dirs/$strain

#     cd  $pav/strain_dirs/$strain
#     echo -e "NAME\tHAP_unphased\n$strain\t$path" > $pav/strain_dirs/$strain/assemblies.tsv
    
#     mkdir -p $pav/strain_dirs/$strain/assemblies
#     cp $assemblies/${strain}.hifi.inbred.fa $pav/strain_dirs/$strain/assemblies
# done

for strain_dir in $pav/strain_dirs/*/; do
    cd $strain_dir
    if ! ls *.vcf.gz &>/dev/null; then 
        rm -r .snakemake data log results temp

        strain_name=$(basename $strain_dir)
        sbatch --job-name="pav${strain_name}" \
        --output="$SLURM/pav_${strain_name}.oe" \
        --error="$SLURM/pav_${strain_name}.rr" \
        /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/pav.sh
    fi 
done
