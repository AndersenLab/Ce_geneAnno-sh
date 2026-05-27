#!/bin/bash

#SBATCH -J paftoolsSV
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/sniffles.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/sniffles.rr 

# for file in *.hifi_reads.sorted.bam; do sbatch --export=alignment=$file /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/sniffles2_SVs.sh individual; done

source activate sniffles

ref_genome="/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
hifi_alignments="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/lr_hifi/alignments"
out_snf="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/sniffles2/snf"
out_vcf="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/sniffles2/vcf"

# Running on indiviudal strains
if [[ $1 == "individual" ]]; then
    sniffles --threads 24 --input $alignment --vcf $out_vcf/$(basename $alignment .sorted.bam).vcf --snf $out_snf/$(basename $alignment .sorted.bam).snf
fi 

# Running in "combine" mode for all strains
if [[ $1 == "combined" ]]; then
    # To create strains.tsv I ran: 
    # for file in *.snf; do echo $PWD/$file >> FILES.tsv; done
    # awk -F'/' '{filename=$NF; sub(/\..*$/, "", filename); print $0 "\t" filename}' FILES.tsv > strains.tsv
    sniffles --input $out_snf/strains.tsv --vcf $out_snf/combined35strain.sniffles2.vcf
fi 