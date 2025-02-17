#!/bin/bash

#SBATCH -J paftoolsSV
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/paftools_SVs.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/paftools_SVs.rr 
 
source activate 2024LabRetreat

ref_genome="/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
asm_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/asssemblies"
out_paf="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/paf"
out_vcf="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf"

for hifi_asm in $asm_dir/*.hifi.inbred.fa; do
    paf="$out_paf/$(basename $hifi_asm .inbred.fa)_asm20.paf"
    vcf="$out_vcf/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    minimap2 -cx asm20 --cs -t 24 $ref_genome $hifi_asm > $paf # c flag - output paf
    
    sort -k6,6 -k8,8n $paf \
    | paftools.js call -l 1000 -L 5000 -f $ref_genome - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
done