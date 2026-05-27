#!/bin/bash

pav="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/elegans"
assemblies="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/elegans/assemblies" 
SLURM="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output"

#tail -n +2 $pav/141_strains.tsv | while IFS=$'\t' read -r strain path; do
#    mkdir -p $pav/strain_dirs/$strain
#    cp $pav/config_Ce.json  $pav/strain_dirs/$strain/config.json
#
#    cd  $pav/strain_dirs/$strain
#    echo -e "NAME\tHAP_unphased\n$strain\t$path" > $pav/strain_dirs/$strain/assemblies.tsv
#    
#    mkdir -p $pav/strain_dirs/$strain/assemblies
#    cp $assemblies/$strain.* $pav/strain_dirs/$strain/assemblies
#done

for strain_dir in $pav/strain_dirs/sixth25/*/; do
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
