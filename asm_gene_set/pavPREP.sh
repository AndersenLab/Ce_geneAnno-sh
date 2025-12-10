#!/bin/bash

pav="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/elegans"
# strain_dirs=($pav/strain_dirs/*/)
# strain_dir=${strain_dirs[$SLURM_ARRAY_TASK_ID-1]}
assemblies="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/elegans/assemblies" 
SLURM="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output"

tail -n +2 $pav/141_strains.tsv | while IFS=$'\t' read -r strain path; do
    mkdir -p $pav/strain_dirs/$strain
    cp $pav/config_Cb.json  $pav/strain_dirs/$strain/config.json

    cd  $pav/strain_dirs/$strain
    echo -e "NAME\tHAP_unphased\n$strain\t$path" > $pav/strain_dirs/$strain/assemblies.tsv
    
    mkdir -p $pav/strain_dirs/$strain/assemblies
    cp $assemblies/$strain.* $pav/strain_dirs/$strain/assemblies
done

for strain_dir in $pav/strain_dirs/*/; do
    cd $strain_dir
    if ! ls *.vcf.gz &>/dev/null; then 
        rm -r .snakemake data log results temp

        strain_name=$(basename $strain_dir)
        sbatch --job-name="pav${strain_name}" \
        --output="$SLURM/pav_${strain_name}briggsae.oe" \
        --error="$SLURM/pav_${strain_name}briggsae.rr" \
        /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/pav.sh
    fi 
done




# for strain_dir in /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/strain_dirs/*/; do
#     cd $strain_dir
#     if ! ls *.vcf.gz &>/dev/null; then 
#         rm -r .snakemake data log results temp

#         strain_name=$(basename $strain_dir)
#         sbatch --job-name="pav${strain_name}" \
#         --output="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pav_${strain_name}LATENCY900.oe" \
#         --error="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pav_${strain_name}LATENCY900.rr" \
#         /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/pav.sh
#     fi 
# done



# sbatch --job-name="pavXECA2581" --output="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pav_XECA2581_LATENCY900.oe" --error="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pav_XECA2581_LATENCY900.rr" /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/pav.sh
# sbatch --job-name="pavECA36" --output="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pav_ECA36_LATENCY900.oe" --error="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pav_ECA36_LATENCY900.rr" /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/pav.sh
# sbatch --job-name="pavECA594" --output="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pav_ECA594_LATENCY900.oe" --error="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/pav_ECA594_LATENCY900.rr" /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/pav.sh

# cd $strain_dir
# if ! ls *.vcf.gz &>/dev/null; then 
#     rm -r .snakemake data log results temp

#     strain_name=$(basename $strain_dir)
#     sbatch --job-name="${strain_name}LATENCY900" \
#     --output="$SLURM/${strain_name}LATENCY900.oe" \
#     --error="$SLURM/${strain_name}LATENCY900.rr" \
#     /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/pav.sh
# fi
