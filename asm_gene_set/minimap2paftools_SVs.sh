#!/bin/bash

#SBATCH -J paftoolsSV
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/june_paftoolsRun.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/june_paftoolsRun.rr 
 
source activate 2024LabRetreat

# ref_genome_elegans="/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
# asm_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies"
# out_paf="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/paf/updated_genomes"
# out_vcf="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/updated_genomes"

# for hifi_asm in $asm_dir/*.hifi.inbred.fa; do
#     paf="$out_paf/$(basename $hifi_asm .inbred.fa)_asm20.paf"
#     vcf="$out_vcf/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
#     if [[ ! -f $vcf ]]; then
#         minimap2 -cx asm20 --cs -t 24 $ref_genome $hifi_asm > $paf # c flag - output paf
        
#         sort -k6,6 -k8,8n $paf \
#         | paftools.js call -l 1000 -L 5000 -f $ref_genome - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
#     fi
# done



# For 06/09_10 briggsae
ref_genome_briggsae="/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome.fa"
asm_dir_CB="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/briggsae"
out_paf0609_0610_CB="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/briggsae/paf/0609_0610"
out_vcf0609_0610_CB="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/briggsae/vcf/0609_0610/"

while read -r CB path_CB; do 
    paf="$out_paf0609_0610_CB/$CB.hifi_asm20.paf"
    vcf="$out_vcf0609_0610_CB/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    if [[ ! -f $vcf ]]; then
        minimap2 -cx asm20 --cs -t 24 $ref_genome_briggsae $asm_dir_CB/$CB.hifi.inbred.fa > $paf # c flag - output paf
        
        sort -k6,6 -k8,8n $paf \
        | paftools.js call -l 1000 -L 5000 -f $ref_genome_briggsae - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
    fi
done < "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/juneCB.tsv"



# For 06/09_10 tropicalis
ref_genome_tropicalis="/vast/eande106/data/c_tropicalis/genomes/NIC58_nanopore/June2021/c_tropicalis.NIC58_nanopore.June2021.genome.fa"
asm_dir_CT="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/tropicalis"
out_paf0609_0610_CT="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/tropicalis/paf/0609_0610"
out_vcf0609_0610_CT="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/tropicalis/vcf/0609_0610/"

while read -r CT path_CT; do 
    paf="$out_paf0609_0610_CT/$CT.hifi_asm20.paf"
    vcf="$out_vcf0609_0610_CT/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    if [[ ! -f $vcf ]]; then
        minimap2 -cx asm20 --cs -t 24 $ref_genome_tropicalis $asm_dir_CT/$CT.hifi.inbred.fa > $paf # c flag - output paf
        
        sort -k6,6 -k8,8n $paf \
        | paftools.js call -l 1000 -L 5000 -f $ref_genome_tropicalis - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
    fi
done < "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/juneCT.tsv"


# For 06/09_10 elegans
ref_genome_elegans="/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
asm_dir_CE="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans"
out_paf0609_0610_CE="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/elegans/paf/0609_0610"
out_vcf0609_0610_CE="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/elegans/vcf/0609_0610/"

while read -r CE path_CE; do 
    paf="$out_paf0609_0610_CE/$CE.hifi_asm20.paf"
    vcf="$out_vcf0609_0610_CE/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    if [[ ! -f $vcf ]]; then
        minimap2 -cx asm20 --cs -t 24 $ref_genome_elegans $asm_dir_CE/$CE.hifi.inbred.fa > $paf # c flag - output paf
        
        sort -k6,6 -k8,8n $paf \
        | paftools.js call -l 1000 -L 5000 -f $ref_genome_elegans - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
    fi
done < "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/juneCE.tsv"


# For 06/09_10 nigoni
ref_genome_nigoni="/vast/eande106/data/c_nigoni/genomes/PRJNA384657/WS276/c_nigoni.PRJNA384657.WS276.genomic.fa.gz"
asm_dir_CN="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/nigoni"
out_paf0609_0610_CN="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/nigoni/paf/0609_0610"
out_vcf0609_0610_CN="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/nigoni/vcf/0609_0610/"

while read -r CN path_CN; do 
    paf="$out_paf0609_0610_CN/$CN.hifi_asm20.paf"
    vcf="$out_vcf0609_0610_CN/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    if [[ ! -f $vcf ]]; then
        minimap2 -cx asm20 --cs -t 24 <(zcat $ref_genome_nigoni) $asm_dir_CN/$CN.hifi.inbred.fa > $paf # c flag - output paf
        
        sort -k6,6 -k8,8n $paf \
        | paftools.js call -l 1000 -L 5000 -f $ref_genome_nigoni - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
    fi
done < "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/juneCN.tsv"